#!/usr/bin/env python3
"""Outward-certify the refined three-band occupations 31 through 256.

The discovery witnesses are stored as one signed byte per composition.  This
checker reconstructs the canonical composition order, recomputes every
positive reference transfer, and bounds each occupation by its largest term
times its exact number of compositions.  The latter deliberately spends at
most log2(C(Q+2,2)) bits to keep the certificate compact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import arb, ctx

import certify_rm2sub_rm49_q1_outward as q1
import certify_rm2sub_rm49_q30_low_outward as q30


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "RM2SUB_RM49_Q31_Q256_INPUT_MANIFEST.json"
CHUNKS = (
    (31, 64),
    (65, 96),
    (97, 128),
    (129, 160),
    (161, 192),
    (193, 208),
    (209, 224),
    (225, 232),
    (233, 240),
    (241, 248),
    (249, 256),
)
GROUPS = (
    ("low", 1, 95, Fraction(1, 4)),
    ("central", 96, 416, Fraction(1, 2)),
    ("high", 417, 512, Fraction(8677722630069483, 10**16)),
)
UNIT_ROUNDOFF = 2.0**-53
MIN_SUBNORMAL = math.nextafter(0.0, math.inf)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest() -> tuple[dict[str, object], dict[str, Path]]:
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    paths: dict[str, Path] = {}
    for row in payload["files"]:
        path = (HERE / row["path"]).resolve()
        if sha256(path) != row["sha256"]:
            raise AssertionError(f"hash mismatch for {row['role']}")
        paths[str(row["role"])] = path
    return payload, paths


def parse_interval(value: str) -> tuple[int, int]:
    try:
        lower_text, upper_text = value.split(":", maxsplit=1)
        result = int(lower_text), int(upper_text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("interval must have form lower:upper") from error
    if result not in CHUNKS:
        raise argparse.ArgumentTypeError(f"interval must be one of {CHUNKS}")
    return result


def compositions(lower: int, upper: int) -> np.ndarray:
    rows = [
        (low, central, occupation - low - central)
        for occupation in range(lower, upper + 1)
        for low in range(occupation + 1)
        for central in range(occupation - low + 1)
    ]
    return np.asarray(rows, dtype=np.int16)


def rounding_inflation(operation_count: int) -> float:
    gamma = operation_count * UNIT_ROUNDOFF
    if gamma >= 1.0:
        raise ArithmeticError("rounding bound is vacuous")
    gamma /= 1.0 - gamma
    return math.nextafter(1.0 / (1.0 - gamma), math.inf)


def upper_convolve(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Upper-bound one positive convolution in binary64."""
    terms = min(len(left), len(right))
    raw = np.convolve(left, right)
    inflated = np.nextafter(
        raw * rounding_inflation(2 * terms + 2), math.inf
    )
    underflow_budget = (2 * terms + 2) * MIN_SUBNORMAL
    return np.nextafter(inflated + underflow_budget, math.inf)


def upper_binomial_transforms(maximum: int) -> list[np.ndarray]:
    transforms: list[np.ndarray] = []
    for _, _, _, probability in GROUPS:
        p = q1.rational(probability)
        one_minus_p = arb(1) - p
        table = np.zeros((maximum + 1, maximum + 1), dtype=np.float64)
        for count in range(maximum + 1):
            value = one_minus_p**count
            table[count, 0] = q1.upper_float(value)
            for live in range(count):
                value *= arb(count - live) * p
                value /= arb(live + 1) * one_minus_p
                table[count, live + 1] = q1.upper_float(value)
        transforms.append(table)
    return transforms


def upper_distributions(
    composition_rows: np.ndarray, maximum: int, transforms: list[np.ndarray]
) -> np.ndarray:
    result = np.zeros((len(composition_rows), maximum + 1), dtype=np.float64)
    for index, composition in enumerate(composition_rows):
        distribution = np.asarray([1.0])
        for count, transform in zip(composition, transforms, strict=True):
            count_int = int(count)
            distribution = upper_convolve(
                distribution, transform[count_int, : count_int + 1]
            )
        result[index, : len(distribution)] = distribution
    if not np.all(np.isfinite(result)) or np.any(result < 0.0):
        raise ArithmeticError("invalid upper distribution")
    return result


def upper_dot(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    terms = left.shape[1]
    raw = left @ right
    inflated = np.nextafter(
        raw * rounding_inflation(2 * terms + 2), math.inf
    )
    underflow_budget = (2 * terms + 2) * MIN_SUBNORMAL
    return np.nextafter(inflated + underflow_budget, math.inf)


def upper_matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    output = np.empty_like(left)
    for row in range(2):
        for column in range(2):
            first = np.nextafter(
                left[:, row, 0] * right[:, 0, column], math.inf
            )
            second = np.nextafter(
                left[:, row, 1] * right[:, 1, column], math.inf
            )
            output[:, row, column] = np.nextafter(first + second, math.inf)
    return output


def normalize_power_of_two(
    matrices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    maximum = np.max(matrices, axis=(1, 2))
    if np.any(maximum <= 0.0) or not np.all(np.isfinite(maximum)):
        raise ArithmeticError("cannot normalize a nonpositive matrix")
    _, exponents = np.frexp(maximum)
    normalized = np.ldexp(matrices, -exponents[:, None, None])
    normalized = np.nextafter(normalized, math.inf)
    return normalized, exponents.astype(np.int64)


def upper_scaled_power_moments(
    matrices: np.ndarray, power: int
) -> tuple[np.ndarray, np.ndarray]:
    count = len(matrices)
    result = np.zeros((count, 2, 2), dtype=np.float64)
    result[:, 0, 0] = 1.0
    result[:, 1, 1] = 1.0
    result_exponents = np.zeros(count, dtype=np.int64)
    base, base_exponents = normalize_power_of_two(matrices)
    exponent = power
    while exponent:
        if exponent & 1:
            result = upper_matmul(result, base)
            result_exponents += base_exponents
            result, normalization = normalize_power_of_two(result)
            result_exponents += normalization
        exponent >>= 1
        if exponent:
            base = upper_matmul(base, base)
            base_exponents *= 2
            base, normalization = normalize_power_of_two(base)
            base_exponents += normalization
    row_sum = np.nextafter(result[:, 0, 0] + result[:, 0, 1], math.inf)
    if np.any(row_sum <= 0.0) or not np.all(np.isfinite(row_sum)):
        raise ArithmeticError("invalid scaled matrix moment")
    return row_sum, result_exponents


def band_envelopes(
    outer_spectrum: dict[int, int],
) -> tuple[list[Fraction], list[int]]:
    envelopes: list[Fraction] = []
    maximizing_weights: list[int] = []
    for _, lower, upper, probability in GROUPS:
        best: Fraction | None = None
        best_weight = -1
        for weight, count in outer_spectrum.items():
            if not lower <= weight <= upper:
                continue
            reference = (
                Fraction(math.comb(q1.OUTER_BITS, weight), 1)
                * probability**weight
                * (1 - probability) ** (q1.OUTER_BITS - weight)
            )
            candidate = Fraction(count, 1) / reference
            if best is None or candidate > best:
                best = candidate
                best_weight = weight
        if best is None:
            raise AssertionError("empty spectrum band")
        envelopes.append(best)
        maximizing_weights.append(best_weight)
    if maximizing_weights != [32, 96, 420]:
        raise AssertionError(f"a density maximum changed: {maximizing_weights}")
    return envelopes, maximizing_weights


def load_chunk(
    lower: int, upper: int, paths: dict[str, Path]
) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    key = f"{lower:03d}_{upper:03d}"
    compact = json.loads(paths[f"compact_{key}"].read_text(encoding="utf-8"))
    if compact["schema"] != "rm2sub-fixed-rm-three-group-common-witness-v1":
        raise AssertionError("compact witness schema changed")
    if compact["parameters"]["occupation_interval"] != [lower, upper]:
        raise AssertionError("compact witness interval changed")
    rows = compositions(lower, upper)
    raw = paths[f"witness_{key}"].read_bytes()
    if hashlib.sha256(raw).hexdigest() != compact["packed_witnesses"]["sha256"]:
        raise AssertionError("packed witness hash disagrees with compact receipt")
    witnesses = np.frombuffer(raw, dtype=np.int8)
    if len(witnesses) != len(rows):
        raise AssertionError("packed witness count is wrong")
    if np.any(witnesses < -40) or np.any(witnesses > 15):
        raise AssertionError("packed witness lies outside the authenticated grid")
    return compact, rows, witnesses


def log_factorial_endpoints() -> tuple[np.ndarray, np.ndarray]:
    lower = np.zeros(q1.OUTER_ROWS + 1, dtype=np.float64)
    upper = np.zeros(q1.OUTER_ROWS + 1, dtype=np.float64)
    factorial = 1
    for value in range(1, q1.OUTER_ROWS + 1):
        factorial *= value
        interval = q1.log2_interval(arb(factorial))
        lower[value] = q1.lower_float(interval)
        upper[value] = q1.upper_float(interval)
    return lower, upper


def upper_location_logs(
    rows: np.ndarray, factorial_lower: np.ndarray, factorial_upper: np.ndarray
) -> np.ndarray:
    occupations = np.sum(rows, axis=1).astype(np.int16)
    result = np.full(len(rows), factorial_upper[q1.OUTER_ROWS])
    result = np.nextafter(
        result - factorial_lower[q1.OUTER_ROWS - occupations], math.inf
    )
    for column in range(3):
        result = np.nextafter(result - factorial_lower[rows[:, column]], math.inf)
    return result


def upper_density_logs(rows: np.ndarray, envelopes: list[Fraction]) -> np.ndarray:
    envelope_logs = np.asarray(
        [q1.upper_float(q1.log2_interval(q1.rational(value))) for value in envelopes]
    )
    result = np.zeros(len(rows), dtype=np.float64)
    for column in range(3):
        term = np.nextafter(rows[:, column] * envelope_logs[column], math.inf)
        result = np.nextafter(result + term, math.inf)
    return result


def upper_log2_exact_float(value: float) -> float:
    return q1.upper_float(q1.log2_interval(arb(value)))


def region_upper_matrices(
    witness: int,
    maximum: int,
    live_spectrum: dict[int, int],
    nonactivation: list[Fraction],
) -> tuple[np.ndarray, float]:
    log_surprisal = q1.rational(Fraction(int(witness), 10))
    surprisal = log_surprisal.exp()
    z = (-surprisal).exp()
    impulses = q30.impulse_matrices(z, live_spectrum, nonactivation, q1.STEP_BITS)
    regions = q30.region_matrices(impulses, maximum)
    upper = np.asarray(
        [
            [
                [q1.upper_float(matrix[row][column]) for column in range(2)]
                for row in range(2)
            ]
            for matrix in regions
        ],
        dtype=np.float64,
    )
    correction = q1.upper_float(
        arb(q1.BAD_WEIGHT) * surprisal / arb(2).log()
    )
    return upper, correction


def process_chunk(
    lower: int,
    upper: int,
    paths: dict[str, Path],
    envelopes: list[Fraction],
    live_spectrum: dict[int, int],
    nonactivation: list[Fraction],
    factorial_lower: np.ndarray,
    factorial_upper: np.ndarray,
) -> dict[str, object]:
    compact, rows, witnesses = load_chunk(lower, upper, paths)
    maximum = upper
    transforms = upper_binomial_transforms(maximum)
    distributions = upper_distributions(rows, maximum, transforms)
    location_logs = upper_location_logs(rows, factorial_lower, factorial_upper)
    density_logs = upper_density_logs(rows, envelopes)
    moment_logs = np.empty(len(rows), dtype=np.float64)
    correction_logs: dict[int, float] = {}
    minimum_reference_entry = math.inf

    for witness_value in np.unique(witnesses):
        witness = int(witness_value)
        selected = np.flatnonzero(witnesses == witness_value)
        region_upper, correction = region_upper_matrices(
            witness, maximum, live_spectrum, nonactivation
        )
        reference_flat = upper_dot(
            distributions[selected], region_upper.reshape((maximum + 1, 4))
        )
        reference = reference_flat.reshape((-1, 2, 2))
        minimum_reference_entry = min(
            minimum_reference_entry, float(np.min(reference))
        )
        row_sum, binary_exponent = upper_scaled_power_moments(
            reference, q1.OUTER_BITS
        )
        local_logs = np.fromiter(
            (upper_log2_exact_float(float(value)) for value in row_sum),
            dtype=np.float64,
            count=len(row_sum),
        )
        moment_logs[selected] = np.nextafter(
            local_logs + binary_exponent.astype(np.float64), math.inf
        )
        correction_logs[witness] = correction

    corrections = np.asarray(
        [correction_logs[int(witness)] for witness in witnesses], dtype=np.float64
    )
    inner_logs = np.minimum(
        0.0, np.nextafter(moment_logs + corrections, math.inf)
    )
    contribution_logs = np.nextafter(location_logs + density_logs, math.inf)
    contribution_logs = np.nextafter(contribution_logs + inner_logs, math.inf)

    occupations = np.sum(rows, axis=1).astype(np.int16)
    occupation_rows: list[dict[str, object]] = []
    for occupation in range(lower, upper + 1):
        selected = occupations == occupation
        maximum_term = float(np.max(contribution_logs[selected]))
        count = math.comb(occupation + 2, 2)
        if int(np.sum(selected)) != count:
            raise AssertionError("composition count changed")
        count_log = q1.upper_float(q1.log2_interval(arb(count)))
        log2_upper = math.nextafter(maximum_term + count_log, math.inf)
        occupation_rows.append(
            {
                "occupation": occupation,
                "composition_count": count,
                "maximum_composition_log2_upper": maximum_term,
                "log2_failure_probability_upper": log2_upper,
                "margin_bits_lower": -log2_upper,
            }
        )

    maximum_occupation = max(
        float(row["log2_failure_probability_upper"]) for row in occupation_rows
    )
    interval_count_log = q1.upper_float(
        q1.log2_interval(arb(upper - lower + 1))
    )
    interval_log2_upper = math.nextafter(
        maximum_occupation + interval_count_log, math.inf
    )
    if interval_log2_upper >= -q1.TARGET_MARGIN_BITS:
        raise AssertionError(f"Q={lower},...,{upper} does not clear 40 bits")
    return {
        "occupation_interval": [lower, upper],
        "composition_count": len(rows),
        "distinct_witness_count": int(len(np.unique(witnesses))),
        "minimum_positive_reference_entry": minimum_reference_entry,
        "log2_failure_probability_upper": interval_log2_upper,
        "margin_bits_lower": -interval_log2_upper,
        "binary64_per_composition_witness_margin_bits": compact[
            "per_composition_witness_interval_margin_bits_diagnostic"
        ],
        "occupation_rows": occupation_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--interval",
        type=parse_interval,
        help="certify one frozen chunk instead of the full Q=31,...,256 interval",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    ctx.prec = 256
    manifest, paths = load_manifest()
    _, q1_paths = q1.load_manifest()
    outer_spectrum = q1.load_outer_spectrum(q1_paths)
    live_spectrum = q1.load_live_spectrum(q1_paths)
    nonactivation = q30.load_nonactivation(q1_paths)
    envelopes, maximizing_weights = band_envelopes(outer_spectrum)
    factorial_lower, factorial_upper = log_factorial_endpoints()

    selected_chunks = (args.interval,) if args.interval is not None else CHUNKS
    chunk_rows = []
    for lower, upper in selected_chunks:
        row = process_chunk(
            lower,
            upper,
            paths,
            envelopes,
            live_spectrum,
            nonactivation,
            factorial_lower,
            factorial_upper,
        )
        chunk_rows.append(row)
        print(
            f"chunk,{lower},{upper},margin_bits_lower,"
            f"{row['margin_bits_lower']:.12f}",
            flush=True,
        )

    maximum_chunk = max(
        float(row["log2_failure_probability_upper"]) for row in chunk_rows
    )
    union_count_log = q1.upper_float(q1.log2_interval(arb(len(chunk_rows))))
    union_log2_upper = math.nextafter(maximum_chunk + union_count_log, math.inf)
    if union_log2_upper >= -q1.TARGET_MARGIN_BITS:
        raise AssertionError("selected chunk union does not clear 40 bits")

    selected_lower = selected_chunks[0][0]
    selected_upper = selected_chunks[-1][1]
    output = args.output or HERE / (
        f"rm2sub_rm49_t64_s14_q{selected_lower:03d}_q{selected_upper:03d}_outward.json"
    )
    payload = {
        "schema": "rm2sub-rm49-refined-band-outward-certificate-v1",
        "status": "OUTWARD_OCCUPATION_CERTIFICATE",
        "claim": {
            "event": (
                f"some occupation-{selected_lower} through occupation-{selected_upper} "
                "message encodes to output weight at most 13107"
            ),
            "log2_failure_probability_upper": union_log2_upper,
            "margin_bits_lower": -union_log2_upper,
            "clears_40_bits": True,
        },
        "coverage": {
            "occupation_interval": [selected_lower, selected_upper],
            "chunk_count": len(chunk_rows),
            "composition_count": sum(int(row["composition_count"]) for row in chunk_rows),
            "complete": True,
        },
        "groups": [
            {
                "name": group[0],
                "weights": [group[1], group[2]],
                "reference_probability_exact": f"{group[3].numerator}/{group[3].denominator}",
                "density_envelope_log2_interval": str(
                    q1.log2_interval(q1.rational(envelope))
                ),
                "maximizing_weight": maximizing_weight,
            }
            for group, envelope, maximizing_weight in zip(
                GROUPS, envelopes, maximizing_weights, strict=True
            )
        ],
        "chunk_rows": chunk_rows,
        "arithmetic": {
            "arb_precision_bits": ctx.prec,
            "binary64_unit_roundoff": UNIT_ROUNDOFF,
            "positive_convolution_rule": "inflate by 1/(1-gamma_(2n+2)) and add an explicit subnormal budget",
            "matrix_power_rule": "round every nonnegative product and sum toward +infinity; normalize only by exact powers of two",
            "aggregation_rule": "maximum term times exact number of terms within each occupation, interval, and chunk union",
        },
        "inputs": {
            "manifest": str(MANIFEST),
            "manifest_sha256": sha256(MANIFEST),
            "checker_sha256": sha256(Path(__file__).resolve()),
            "files": manifest["files"],
        },
        "scope": (
            "This certificate covers only the displayed occupation interval. "
            "The full distance theorem also requires the separate Q=1,...,30 certificates."
        ),
    }
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status,{payload['status']}")
    print(f"margin_bits_lower,{-union_log2_upper:.12f}")
    print(f"wrote,{output}")


if __name__ == "__main__":
    main()

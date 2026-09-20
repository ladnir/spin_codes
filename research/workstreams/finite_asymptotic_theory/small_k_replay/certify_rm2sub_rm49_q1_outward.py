#!/usr/bin/env python3
"""Certify the fixed-RM(4,9), RM2Sub-t64,s14 occupation-one bound.

The checker reads a finite witness table selected by a binary64 discovery
run. It canonicalizes each log-surprisal to an exact tenth. All probability
calculations then use python-flint Arb interval arithmetic.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx, fmpq


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "RM2SUB_RM49_Q1_INPUT_MANIFEST.json"
OUTPUT = HERE / "rm2sub_rm49_t64_s14_q1_outward.json"

OUTER_BITS = 512
OUTER_DIMENSION = 256
OUTER_ROWS = 256
MESSAGE_BITS = OUTER_DIMENSION * OUTER_ROWS
OUTPUT_BITS = OUTER_BITS * OUTER_ROWS
BAD_WEIGHT = 13107
STEP_BITS = 64
STATE_BITS = 14
EPOCHS_PER_REGION = OUTER_ROWS // STEP_BITS
TARGET_MARGIN_BITS = 40


Matrix = tuple[tuple[arb, arb], tuple[arb, arb]]
Vector = tuple[arb, arb]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rational(value: Fraction | int) -> arb:
    if isinstance(value, int):
        return arb(value)
    return arb(fmpq(value.numerator, value.denominator))


def upper_float(value: arb) -> float:
    return math.nextafter(float(value.upper()), math.inf)


def lower_float(value: arb) -> float:
    return math.nextafter(float(value.lower()), -math.inf)


def log2_interval(value: arb) -> arb:
    if not value > 0:
        raise ArithmeticError("logarithm input is not definitely positive")
    return value.log() / arb(2).log()


def zero_matrix() -> Matrix:
    return ((arb(0), arb(0)), (arb(0), arb(0)))


def identity_matrix() -> Matrix:
    return ((arb(1), arb(0)), (arb(0), arb(1)))


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    return (
        (left[0][0] + right[0][0], left[0][1] + right[0][1]),
        (left[1][0] + right[1][0], left[1][1] + right[1][1]),
    )


def matrix_scale(value: Matrix, scalar: arb) -> Matrix:
    return (
        (value[0][0] * scalar, value[0][1] * scalar),
        (value[1][0] * scalar, value[1][1] * scalar),
    )


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def matrix_power(value: Matrix, exponent: int) -> Matrix:
    result = identity_matrix()
    factor = value
    power = exponent
    while power:
        if power & 1:
            result = matrix_multiply(result, factor)
        power >>= 1
        if power:
            factor = matrix_multiply(factor, factor)
    return result


def vector_matrix(value: Vector, matrix: Matrix) -> Vector:
    return (
        value[0] * matrix[0][0] + value[1] * matrix[1][0],
        value[0] * matrix[0][1] + value[1] * matrix[1][1],
    )


def vector_add(left: Vector, right: Vector) -> Vector:
    return (left[0] + right[0], left[1] + right[1])


def load_manifest() -> tuple[dict[str, object], dict[str, Path]]:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    resolved: dict[str, Path] = {}
    for row in manifest["files"]:
        path = (HERE / row["path"]).resolve()
        actual = sha256(path)
        expected = str(row["sha256"])
        if actual != expected:
            raise AssertionError(
                f"hash mismatch for {row['role']}: {actual} != {expected}"
            )
        resolved[str(row["role"])] = path
    return manifest, resolved


def load_outer_spectrum(paths: dict[str, Path]) -> dict[int, int]:
    payload = json.loads(
        paths["outer_spectrum_certificate_input"].read_text(encoding="utf-8")
    )
    counts = {int(weight): int(count) for weight, count in payload["weight_counts"].items()}
    with paths["authenticated_outer_spectrum_csv"].open(
        newline="", encoding="utf-8"
    ) as source:
        csv_counts = {
            int(row["weight"]): int(row["count"]) for row in csv.DictReader(source)
        }
    if counts != csv_counts:
        raise AssertionError("converted RM spectrum does not match the authenticated CSV")
    if counts.get(0) != 1 or sum(counts.values()) != 1 << OUTER_DIMENSION:
        raise AssertionError("RM spectrum has the wrong mass")
    nonzero = {weight: count for weight, count in counts.items() if weight and count}
    if min(nonzero) != 32 or max(nonzero) != OUTER_BITS:
        raise AssertionError("RM spectrum has the wrong support")
    return nonzero


def load_live_spectrum(paths: dict[str, Path]) -> dict[int, int]:
    payload = json.loads(paths["rm2sub_a_spectrum"].read_text(encoding="utf-8"))
    counts = {int(row["weight"]): int(row["count"]) for row in payload["spectrum"]}
    if counts.get(0) != 1 or sum(counts.values()) != 1 << STATE_BITS:
        raise AssertionError("A spectrum has the wrong mass")
    selection = json.loads(paths["rm2sub_map_selection"].read_text(encoding="utf-8"))
    selected = selection["selected"]
    if (
        int(selected["minimum_A_distance"]) != 24
        or not bool(selected["BA_zero"])
        or not bool(selected["distinct_nonzero_columns"])
    ):
        raise AssertionError("selected RM2Sub map does not satisfy the Q1 assumptions")
    return {weight: count for weight, count in counts.items() if weight}


def verify_weight_one_nonactivation(paths: dict[str, Path]) -> Fraction:
    payload = json.loads(
        paths["rm2sub_b_kernel_spectrum"].read_text(encoding="utf-8")
    )
    rows = {int(row["total_weight"]): row for row in payload["by_total_weight"]}
    row = rows[1]
    if int(row["kernel_words"]) != 0 or int(row["shell_size"]) != STEP_BITS:
        raise AssertionError("B kernel does not have exact weight-one nonactivation zero")
    if int(payload["checks"]["minimum_kernel_distance"]) != 6:
        raise AssertionError("B kernel minimum distance changed")
    return Fraction(int(row["kernel_words"]), int(row["shell_size"]))


def load_witnesses(
    path: Path, outer_spectrum: dict[int, int]
) -> tuple[dict[int, int], dict[int, dict[str, object]]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    parameters = payload["parameters"]
    expected = {
        "message_bits": MESSAGE_BITS,
        "output_bits": OUTPUT_BITS,
        "outer_bits": OUTER_BITS,
        "outer_dimension": OUTER_DIMENSION,
        "outer_blocks": OUTER_ROWS,
        "step_bits": STEP_BITS,
        "state_bits": STATE_BITS,
        "epochs_per_region": EPOCHS_PER_REGION,
        "target_distance": BAD_WEIGHT,
    }
    for key, value in expected.items():
        if int(parameters[key]) != value:
            raise AssertionError(f"witness parameter mismatch for {key}")
    rows = {int(row["outer_weight"]): row for row in payload["weight_rows"]}
    if set(rows) != set(outer_spectrum):
        raise AssertionError("witness weights do not equal the nonzero RM spectrum support")
    witness_tenths: dict[int, int] = {}
    for weight, row in rows.items():
        raw = float(row["best_log_surprisal"])
        index = round(10.0 * raw)
        if abs(raw - index / 10.0) >= 1e-12 or not -120 <= index <= 0:
            raise AssertionError(f"weight {weight} has a noncanonical witness")
        witness_tenths[weight] = index
    return witness_tenths, rows


def live_moments_for_zero_and_one(z: arb, spectrum: dict[int, int]) -> tuple[arb, arb]:
    denominator = arb((1 << STATE_BITS) - 1)
    moment_zero = arb(0)
    moment_one = arb(0)
    for weight, count in spectrum.items():
        probability_mass = arb(count) / denominator
        moment_zero += probability_mass * (z**weight)
        kernel_one = (
            arb(weight) * (z ** (weight - 1))
            + arb(STEP_BITS - weight) * (z ** (weight + 1))
        ) / STEP_BITS
        moment_one += probability_mass * kernel_one
    return moment_zero, moment_one


def epoch_matrices(z: arb, live_spectrum: dict[int, int], nonactivation_one: Fraction) -> tuple[Matrix, Matrix]:
    denominator = arb((1 << STATE_BITS) - 1)
    punctured = denominator / arb((1 << STATE_BITS) - 2)
    moment_zero, moment_one = live_moments_for_zero_and_one(z, live_spectrum)
    zero = (
        (arb(1), arb(0)),
        (arb(0), punctured * moment_zero),
    )
    output_factor = z
    averaged = punctured * moment_one
    active = (
        (rational(nonactivation_one) * output_factor, output_factor),
        (averaged / denominator, averaged),
    )
    return zero, active


def region_matrices(zero_epoch: Matrix, active_epoch: Matrix) -> tuple[Matrix, Matrix]:
    powers = [identity_matrix()]
    for _ in range(EPOCHS_PER_REGION):
        powers.append(matrix_multiply(powers[-1], zero_epoch))
    inactive = powers[EPOCHS_PER_REGION]
    active = zero_matrix()
    for epoch in range(EPOCHS_PER_REGION):
        term = matrix_multiply(
            matrix_multiply(powers[epoch], active_epoch),
            powers[EPOCHS_PER_REGION - 1 - epoch],
        )
        active = matrix_add(active, term)
    return inactive, matrix_scale(active, rational(Fraction(1, EPOCHS_PER_REGION)))


def coefficient_moments(
    inactive: Matrix, active: Matrix, requested_weights: set[int]
) -> dict[int, arb]:
    """Return exact-subset averaged transfer moments for requested weights."""
    result: dict[int, arb] = {}
    low = {weight for weight in requested_weights if weight <= OUTER_BITS // 2}
    high = requested_weights - low

    def run(base: Matrix, marked: Matrix, targets: dict[int, int]) -> None:
        if not targets:
            return
        maximum_degree = max(targets)
        coefficients: list[Vector] = [(arb(1), arb(0))] + [
            (arb(0), arb(0)) for _ in range(maximum_degree)
        ]
        for row_index in range(OUTER_BITS):
            upper = min(maximum_degree, row_index + 1)
            updated: list[Vector] = [(arb(0), arb(0)) for _ in range(maximum_degree + 1)]
            for degree in range(upper + 1):
                if degree <= row_index:
                    updated[degree] = vector_add(
                        updated[degree], vector_matrix(coefficients[degree], base)
                    )
                if degree and degree - 1 <= row_index:
                    updated[degree] = vector_add(
                        updated[degree],
                        vector_matrix(coefficients[degree - 1], marked),
                    )
            coefficients = updated
        for degree, output_weight in targets.items():
            vector = coefficients[degree]
            total = vector[0] + vector[1]
            result[output_weight] = total / math.comb(OUTER_BITS, output_weight)

    run(inactive, active, {weight: weight for weight in low})
    run(active, inactive, {OUTER_BITS - weight: weight for weight in high})
    return result


def main() -> None:
    ctx.prec = 256
    manifest, paths = load_manifest()
    outer_spectrum = load_outer_spectrum(paths)
    live_spectrum = load_live_spectrum(paths)
    nonactivation_one = verify_weight_one_nonactivation(paths)
    witness_tenths, discovery_rows = load_witnesses(
        paths["binary64_witness_discovery"], outer_spectrum
    )

    weights_by_witness: dict[int, set[int]] = defaultdict(set)
    for weight, witness in witness_tenths.items():
        weights_by_witness[witness].add(weight)

    contributions: dict[int, arb] = {}
    rows: list[dict[str, object]] = []
    maximum_float_delta = 0.0
    for witness in sorted(weights_by_witness):
        log_surprisal = rational(Fraction(witness, 10))
        surprisal = log_surprisal.exp()
        z = (-surprisal).exp()
        zero_epoch, active_epoch = epoch_matrices(
            z, live_spectrum, nonactivation_one
        )
        inactive_region, active_region = region_matrices(zero_epoch, active_epoch)
        moments = coefficient_moments(
            inactive_region, active_region, weights_by_witness[witness]
        )
        correction = (arb(BAD_WEIGHT) * surprisal).exp()
        for weight in sorted(weights_by_witness[witness]):
            raw_inner = moments[weight] * correction
            inner = raw_inner if raw_inner < arb(1) else arb(1)
            contribution = arb(OUTER_ROWS * outer_spectrum[weight]) * inner
            contributions[weight] = contribution
            pointwise_log2 = log2_interval(contribution)
            binary64 = float(discovery_rows[weight]["pointwise_log2_upper"])
            maximum_float_delta = max(
                maximum_float_delta,
                abs(upper_float(pointwise_log2) - binary64),
            )
            rows.append(
                {
                    "outer_weight": weight,
                    "outer_multiplicity": outer_spectrum[weight],
                    "log_surprisal_exact": f"{witness}/10",
                    "pointwise_log2_interval": str(pointwise_log2),
                    "pointwise_log2_upper": upper_float(pointwise_log2),
                    "pointwise_margin_bits_lower": -upper_float(pointwise_log2),
                    "binary64_pointwise_log2": binary64,
                }
            )

    if set(contributions) != set(outer_spectrum):
        raise AssertionError("outward calculation omitted an RM weight")
    total = sum(contributions.values(), arb(0))
    threshold = rational(Fraction(1, 1 << TARGET_MARGIN_BITS))
    if not total < threshold:
        raise AssertionError("occupation-one failure bound does not clear 40 bits")
    total_log2 = log2_interval(total)
    margin_lower = -upper_float(total_log2)
    if margin_lower <= TARGET_MARGIN_BITS:
        raise AssertionError("reported lower margin does not clear 40 bits")
    if maximum_float_delta >= 1e-6:
        raise AssertionError("outward result does not match binary64 discovery")

    rows.sort(key=lambda row: int(row["outer_weight"]))
    dominant = min(rows, key=lambda row: float(row["pointwise_margin_bits_lower"]))
    payload = {
        "schema": "rm2sub-rm49-q1-outward-certificate-v1",
        "status": "OUTWARD_Q1_CERTIFICATE",
        "claim": {
            "event": "some occupation-one message encodes to output weight at most 13107",
            "failure_probability_upper_interval": str(total),
            "log2_failure_probability_interval": str(total_log2),
            "margin_bits_lower": margin_lower,
            "clears_40_bits": True,
        },
        "construction": {
            "message_bits": MESSAGE_BITS,
            "output_bits": OUTPUT_BITS,
            "bad_weight_upper": BAD_WEIGHT,
            "outer_rows": OUTER_ROWS,
            "outer": "one fixed RM(4,9) [512,256,32] code repeated in every row",
            "routing": "independent uniform row-coordinate permutations and independent uniform permutations in 512 transposed regions",
            "inner": "fixed audited RM2Sub t=64,s=14 A/B pair with independent nonzero field scalar per epoch",
            "epochs_per_region": EPOCHS_PER_REGION,
        },
        "probability_space": "routing permutations and nonzero epoch scalars; the RM outer and RM2Sub A/B maps are fixed",
        "arithmetic": {
            "library": "python-flint Arb",
            "precision_bits": ctx.prec,
            "all_certificate_operations": "positive Arb interval arithmetic",
            "witnesses": "fixed exact tenths authenticated from the binary64 discovery receipt",
            "maximum_pointwise_log2_delta_from_binary64": maximum_float_delta,
        },
        "coverage": {
            "occupation": 1,
            "nonzero_outer_weight_count": len(outer_spectrum),
            "distinct_witness_count": len(weights_by_witness),
            "outer_weights": sorted(outer_spectrum),
        },
        "dominant_row": dominant,
        "weight_rows": rows,
        "inputs": {
            "manifest": str(MANIFEST),
            "manifest_sha256": sha256(MANIFEST),
            "checker_sha256": sha256(Path(__file__).resolve()),
            "files": manifest["files"],
        },
        "scope": "This certificate covers occupation one only. Occupations two through 256 require their own outward certificates before an end-to-end distance theorem is stated.",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status,{payload['status']}")
    print(f"margin_bits_lower,{margin_lower:.12f}")
    print(f"dominant_outer_weight,{dominant['outer_weight']}")
    print(f"distinct_witness_count,{len(weights_by_witness)}")
    print(f"wrote,{OUTPUT}")


if __name__ == "__main__":
    main()

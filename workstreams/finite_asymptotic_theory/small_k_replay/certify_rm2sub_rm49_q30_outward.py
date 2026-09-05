#!/usr/bin/env python3
"""Certify all 496 refined-band compositions at occupation 30."""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

import certify_rm2sub_rm49_q1_outward as q1
import certify_rm2sub_rm49_q30_low_outward as q30


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "RM2SUB_RM49_Q30_INPUT_MANIFEST.json"
OUTPUT = HERE / "rm2sub_rm49_t64_s14_q30_outward.json"
GROUPS = (
    ("low", 1, 95, Fraction(1, 4)),
    ("central", 96, 416, Fraction(1, 2)),
    ("high", 417, 512, Fraction(8677722630069483, 10**16)),
)


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


def compositions(total: int):
    for low in range(total + 1):
        for central in range(total - low + 1):
            yield low, central, total - low - central


def load_witnesses(path: Path) -> tuple[dict[tuple[int, int, int], int], dict[tuple[int, int, int], float], float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    expected = set(compositions(q30.OCCUPATION))
    witnesses: dict[tuple[int, int, int], int] = {}
    binary64: dict[tuple[int, int, int], float] = {}
    for row in payload["composition_rows"]:
        if int(row["occupation"]) != q30.OCCUPATION:
            raise AssertionError("witness receipt contains a different occupation")
        composition = (
            int(row["low_rows"]),
            int(row["central_rows"]),
            int(row["high_rows"]),
        )
        raw = float(row["log_surprisal"])
        witness = round(10.0 * raw)
        if abs(raw - witness / 10.0) >= 1e-12:
            raise AssertionError("composition witness is not on the exact tenth grid")
        witnesses[composition] = witness
        binary64[composition] = float(row["log2_upper_diagnostic"])
    if set(witnesses) != expected or len(witnesses) != 496:
        raise AssertionError("witness table does not cover exactly 496 Q30 compositions")
    return witnesses, binary64, float(payload["interval_log2_upper_diagnostic"])


def band_envelopes(outer_spectrum: dict[int, int]) -> tuple[list[Fraction], list[int]]:
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
    if (
        maximizing_weights[0] != 32
        or maximizing_weights[1] not in (96, 416)
        or maximizing_weights[2] != 420
    ):
        raise AssertionError(
            f"a refined-band density maximum changed: {maximizing_weights}"
        )
    return envelopes, maximizing_weights


def convolve(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    result = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for first, first_value in enumerate(left):
        for second, second_value in enumerate(right):
            result[first + second] += first_value * second_value
    return result


def binomial_distribution(count: int, probability: Fraction) -> list[Fraction]:
    return [
        Fraction(math.comb(count, live), 1)
        * probability**live
        * (1 - probability) ** (count - live)
        for live in range(count + 1)
    ]


def reference_region(
    regions: list[q1.Matrix], composition: tuple[int, int, int]
) -> q1.Matrix:
    distribution = [Fraction(1)]
    for count, (_, _, _, probability) in zip(composition, GROUPS, strict=True):
        distribution = convolve(distribution, binomial_distribution(count, probability))
    if sum(distribution) != 1 or len(distribution) != q30.OCCUPATION + 1:
        raise AssertionError("reference live-count distribution has the wrong mass")
    result = q1.zero_matrix()
    for live, probability in enumerate(distribution):
        result = q1.matrix_add(
            result, q1.matrix_scale(regions[live], q1.rational(probability))
        )
    return result


def location_count(composition: tuple[int, int, int]) -> int:
    remaining = q1.OUTER_ROWS
    result = 1
    for count in composition:
        result *= math.comb(remaining, count)
        remaining -= count
    return result


def main() -> None:
    ctx.prec = 256
    manifest, local_paths = load_manifest()
    _, q1_paths = q1.load_manifest()
    outer_spectrum = q1.load_outer_spectrum(q1_paths)
    live_spectrum = q1.load_live_spectrum(q1_paths)
    nonactivation = q30.load_nonactivation(q1_paths)
    witnesses, binary64_rows, binary64_union = load_witnesses(
        local_paths["binary64_composition_witnesses"]
    )
    envelopes, maximizing_weights = band_envelopes(outer_spectrum)

    by_witness: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    for composition, witness in witnesses.items():
        by_witness[witness].append(composition)

    contributions: dict[tuple[int, int, int], arb] = {}
    rows: list[dict[str, object]] = []
    maximum_pointwise_delta = 0.0
    for witness in sorted(by_witness):
        log_surprisal = q1.rational(Fraction(witness, 10))
        surprisal = log_surprisal.exp()
        z = (-surprisal).exp()
        impulses = q30.impulse_matrices(
            z, live_spectrum, nonactivation, q30.OCCUPATION
        )
        regions = q30.region_matrices(impulses, q30.OCCUPATION)
        correction = (arb(q1.BAD_WEIGHT) * surprisal).exp()
        for composition in by_witness[witness]:
            reference = reference_region(regions, composition)
            powered = q1.matrix_power(reference, q1.OUTER_BITS)
            moment = powered[0][0] + powered[0][1]
            raw_inner = moment * correction
            inner = raw_inner if raw_inner < arb(1) else arb(1)
            density = arb(1)
            for count, envelope in zip(composition, envelopes, strict=True):
                density *= q1.rational(envelope) ** count
            contribution = arb(location_count(composition)) * density * inner
            contributions[composition] = contribution
            log2_value = q1.log2_interval(contribution)
            log2_upper = q1.upper_float(log2_value)
            delta = abs(log2_upper - binary64_rows[composition])
            maximum_pointwise_delta = max(maximum_pointwise_delta, delta)
            rows.append(
                {
                    "low_rows": composition[0],
                    "central_rows": composition[1],
                    "high_rows": composition[2],
                    "log_surprisal_exact": f"{witness}/10",
                    "log2_upper_interval": str(log2_value),
                    "log2_upper": log2_upper,
                    "margin_bits_lower": -log2_upper,
                    "binary64_log2_upper": binary64_rows[composition],
                }
            )

    if set(contributions) != set(compositions(q30.OCCUPATION)):
        raise AssertionError("outward calculation omitted a Q30 composition")
    if maximum_pointwise_delta >= 1e-6:
        raise AssertionError("outward composition does not match binary64 discovery")
    total = sum(contributions.values(), arb(0))
    total_log2 = q1.log2_interval(total)
    total_upper = q1.upper_float(total_log2)
    union_delta = abs(total_upper - binary64_union)
    if union_delta >= 1e-6:
        raise AssertionError("outward Q30 union does not match binary64 discovery")
    if not total < arb(1) / (1 << 1000):
        raise AssertionError("occupation-30 union did not retain ample margin")

    rows.sort(key=lambda row: (row["low_rows"], row["central_rows"], row["high_rows"]))
    dominant = min(rows, key=lambda row: float(row["margin_bits_lower"]))
    payload = {
        "schema": "rm2sub-rm49-q30-outward-certificate-v1",
        "status": "OUTWARD_OCCUPATION_CERTIFICATE",
        "claim": {
            "event": "some occupation-30 message encodes to output weight at most 13107",
            "failure_probability_upper_interval": str(total),
            "log2_failure_probability_interval": str(total_log2),
            "margin_bits_lower": -total_upper,
        },
        "coverage": {
            "occupation": q30.OCCUPATION,
            "composition_count": len(rows),
            "distinct_witness_count": len(by_witness),
            "complete": True,
        },
        "groups": [
            {
                "name": group[0],
                "weights": [group[1], group[2]],
                "reference_probability_exact": f"{group[3].numerator}/{group[3].denominator}",
                "density_envelope_definition": "A_w/(C(512,w)*p^w*(1-p)^(512-w)) at the displayed maximizing weight",
                "density_envelope_log2_interval": str(
                    q1.log2_interval(q1.rational(envelope))
                ),
                "maximizing_weight": maximizing_weight,
            }
            for group, envelope, maximizing_weight in zip(
                GROUPS, envelopes, maximizing_weights, strict=True
            )
        ],
        "arithmetic": {
            "library": "python-flint Arb",
            "precision_bits": ctx.prec,
            "maximum_pointwise_log2_delta_from_binary64": maximum_pointwise_delta,
            "union_log2_delta_from_binary64": union_delta,
        },
        "dominant_composition": dominant,
        "composition_rows": rows,
        "inputs": {
            "manifest": str(MANIFEST),
            "manifest_sha256": sha256(MANIFEST),
            "checker_sha256": sha256(Path(__file__).resolve()),
            "files": manifest["files"],
        },
        "scope": "This certificate covers every composition at occupation 30. Other occupations require separate outward certificates before the full distance theorem is stated.",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status,{payload['status']}")
    print(f"margin_bits_lower,{-total_upper:.12f}")
    print(f"composition_count,{len(rows)}")
    print(f"maximum_pointwise_delta,{maximum_pointwise_delta:.3e}")
    print(f"union_delta,{union_delta:.3e}")
    print(f"wrote,{OUTPUT}")


if __name__ == "__main__":
    main()

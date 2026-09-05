#!/usr/bin/env python3
"""Certify occupation two with exact-shell row tilts and Arb intervals."""

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
MANIFEST = HERE / "RM2SUB_RM49_Q2_SHELL_INPUT_MANIFEST.json"
OUTPUT = HERE / "rm2sub_rm49_t64_s14_q2_shells_outward.json"
OCCUPATION = 2


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


def load_witnesses(
    path: Path, outer_weights: set[int]
) -> tuple[dict[tuple[int, int], tuple[int, int]], dict[tuple[int, int], float], float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    parameters = payload["parameters"]
    expected_parameters = {
        "message_bits": q1.MESSAGE_BITS,
        "output_bits": q1.OUTPUT_BITS,
        "outer_rows": q1.OUTER_ROWS,
        "outer_block_bits": q1.OUTER_BITS,
        "step_bits": q1.STEP_BITS,
        "state_bits": q1.STATE_BITS,
        "bad_weight": q1.BAD_WEIGHT,
        "occupation": OCCUPATION,
    }
    for key, value in expected_parameters.items():
        if int(parameters[key]) != value:
            raise AssertionError(f"witness parameter mismatch for {key}")
    expected_pairs = {
        (first, second) for first in outer_weights for second in outer_weights
    }
    witnesses: dict[tuple[int, int], tuple[int, int]] = {}
    binary64: dict[tuple[int, int], float] = {}
    for row in payload["pair_rows"]:
        pair = (int(row["first_weight"]), int(row["second_weight"]))
        numerator, denominator = str(row["log_surprisal_exact"]).split("/")
        if denominator != "10":
            raise AssertionError("output witness is not an exact tenth")
        output_tenths = int(numerator)
        raw_shift = float(row["logit_shift"])
        shift_quarters = round(4.0 * raw_shift)
        if abs(raw_shift - shift_quarters / 4.0) >= 1e-12:
            raise AssertionError("row-weight witness is not an exact quarter")
        witnesses[pair] = (output_tenths, shift_quarters)
        binary64[pair] = float(row["log2_upper_diagnostic"])
    if set(witnesses) != expected_pairs:
        raise AssertionError("witness table does not cover every ordered weight pair")
    return witnesses, binary64, float(payload["log2_upper_diagnostic"])


def shifted_probability(weight: int, shift_quarters: int) -> arb:
    if weight == q1.OUTER_BITS:
        return arb(1)
    odds = arb(weight) / arb(q1.OUTER_BITS - weight)
    shift = q1.rational(Fraction(shift_quarters, 4)).exp()
    shifted_odds = odds * shift
    return shifted_odds / (arb(1) + shifted_odds)


def shell_density(
    weight: int, count: int, probability: arb
) -> arb:
    reference = (
        arb(math.comb(q1.OUTER_BITS, weight))
        * (probability**weight)
        * ((arb(1) - probability) ** (q1.OUTER_BITS - weight))
    )
    if not reference > 0:
        raise ArithmeticError("shell reference mass is not definitely positive")
    return arb(count) / reference


def reference_region(
    regions: list[q1.Matrix], first_probability: arb, second_probability: arb
) -> q1.Matrix:
    probabilities = (
        (arb(1) - first_probability) * (arb(1) - second_probability),
        first_probability * (arb(1) - second_probability)
        + (arb(1) - first_probability) * second_probability,
        first_probability * second_probability,
    )
    result = q1.zero_matrix()
    for live, probability in enumerate(probabilities):
        result = q1.matrix_add(
            result, q1.matrix_scale(regions[live], probability)
        )
    return result


def main() -> None:
    ctx.prec = 256
    manifest, local_paths = load_manifest()
    _, q1_paths = q1.load_manifest()
    outer_spectrum = q1.load_outer_spectrum(q1_paths)
    live_spectrum = q1.load_live_spectrum(q1_paths)
    nonactivation = q30.load_nonactivation(q1_paths)
    witnesses, binary64_rows, binary64_union = load_witnesses(
        local_paths["binary64_shell_witnesses"], set(outer_spectrum)
    )

    by_output_witness: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for pair, (output_witness, _) in witnesses.items():
        by_output_witness[output_witness].append(pair)

    probability_cache: dict[tuple[int, int], arb] = {}
    density_cache: dict[tuple[int, int], arb] = {}
    contributions: dict[tuple[int, int], arb] = {}
    rows: list[dict[str, object]] = []
    maximum_pointwise_delta = 0.0
    location = math.comb(q1.OUTER_ROWS, OCCUPATION)
    for output_witness in sorted(by_output_witness):
        log_surprisal = q1.rational(Fraction(output_witness, 10))
        surprisal = log_surprisal.exp()
        z = (-surprisal).exp()
        impulses = q30.impulse_matrices(z, live_spectrum, nonactivation, OCCUPATION)
        regions = q30.region_matrices(impulses, OCCUPATION)
        correction = (arb(q1.BAD_WEIGHT) * surprisal).exp()
        for pair in by_output_witness[output_witness]:
            _, shift_quarters = witnesses[pair]
            probabilities: list[arb] = []
            densities: list[arb] = []
            for weight in pair:
                key = (weight, shift_quarters)
                if key not in probability_cache:
                    probability_cache[key] = shifted_probability(
                        weight, shift_quarters
                    )
                    density_cache[key] = shell_density(
                        weight, outer_spectrum[weight], probability_cache[key]
                    )
                probabilities.append(probability_cache[key])
                densities.append(density_cache[key])
            reference = reference_region(regions, probabilities[0], probabilities[1])
            powered = q1.matrix_power(reference, q1.OUTER_BITS)
            moment = powered[0][0] + powered[0][1]
            raw_inner = moment * correction
            inner = raw_inner if raw_inner < arb(1) else arb(1)
            contribution = arb(location) * densities[0] * densities[1] * inner
            contributions[pair] = contribution
            log2_value = q1.log2_interval(contribution)
            log2_upper = q1.upper_float(log2_value)
            delta = abs(log2_upper - binary64_rows[pair])
            maximum_pointwise_delta = max(maximum_pointwise_delta, delta)
            rows.append(
                {
                    "first_weight": pair[0],
                    "second_weight": pair[1],
                    "log_surprisal_exact": f"{output_witness}/10",
                    "logit_shift_exact": f"{shift_quarters}/4",
                    "log2_upper_interval": str(log2_value),
                    "log2_upper": log2_upper,
                    "margin_bits_lower": -log2_upper,
                    "binary64_log2_upper": binary64_rows[pair],
                }
            )

    if set(contributions) != set(witnesses):
        raise AssertionError("outward calculation omitted an ordered weight pair")
    if maximum_pointwise_delta >= 1e-6:
        raise AssertionError("outward pair bound does not match binary64 discovery")
    total = sum(contributions.values(), arb(0))
    total_log2 = q1.log2_interval(total)
    total_upper = q1.upper_float(total_log2)
    union_delta = abs(total_upper - binary64_union)
    if union_delta >= 1e-6:
        raise AssertionError("outward Q2 union does not match binary64 discovery")
    if not total < arb(1) / (1 << q1.TARGET_MARGIN_BITS):
        raise AssertionError("occupation-two failure bound does not clear 40 bits")

    rows.sort(key=lambda row: (row["first_weight"], row["second_weight"]))
    dominant = min(rows, key=lambda row: float(row["margin_bits_lower"]))
    payload = {
        "schema": "rm2sub-rm49-q2-shells-outward-certificate-v1",
        "status": "OUTWARD_OCCUPATION_CERTIFICATE",
        "claim": {
            "event": "some occupation-two message encodes to output weight at most 13107",
            "failure_probability_upper_interval": str(total),
            "log2_failure_probability_interval": str(total_log2),
            "margin_bits_lower": -total_upper,
            "clears_40_bits": True,
        },
        "construction": {
            "message_bits": q1.MESSAGE_BITS,
            "output_bits": q1.OUTPUT_BITS,
            "outer_rows": q1.OUTER_ROWS,
            "outer": "one fixed RM(4,9) [512,256,32] code repeated in every row",
            "routing": "independent uniform row-coordinate permutations and independent uniform permutations in 512 transposed regions",
            "inner": "fixed audited RM2Sub t=64,s=14 A/B pair with independent nonzero field scalar per epoch",
        },
        "reduction": "For each ordered exact-weight pair, apply a positive row-weight Chernoff tilt to both fixed shells and the positive forced-support RM2Sub transfer. Sum all ordered shell pairs and all unordered active-row locations.",
        "coverage": {
            "occupation": OCCUPATION,
            "ordered_weight_pair_count": len(rows),
            "distinct_output_witness_count": len(by_output_witness),
            "complete": True,
        },
        "arithmetic": {
            "library": "python-flint Arb",
            "precision_bits": ctx.prec,
            "maximum_pointwise_log2_delta_from_binary64": maximum_pointwise_delta,
            "union_log2_delta_from_binary64": union_delta,
        },
        "dominant_pair": dominant,
        "pair_rows": rows,
        "inputs": {
            "manifest": str(MANIFEST),
            "manifest_sha256": sha256(MANIFEST),
            "checker_sha256": sha256(Path(__file__).resolve()),
            "files": manifest["files"],
        },
        "scope": "This certificate covers occupation two only. Occupations three through 256 require their outward receipts before the end-to-end theorem is stated.",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status,{payload['status']}")
    print(f"margin_bits_lower,{-total_upper:.12f}")
    print(f"ordered_weight_pair_count,{len(rows)}")
    print(f"maximum_pointwise_delta,{maximum_pointwise_delta:.3e}")
    print(f"union_delta,{union_delta:.3e}")
    print(f"wrote,{OUTPUT}")


if __name__ == "__main__":
    main()

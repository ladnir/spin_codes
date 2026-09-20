#!/usr/bin/env python3
"""Certify the known RM(5,11) occupation-one shells below weight 128."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

from flint import arb, ctx


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

import certify_single_random_constituent_sparse_outward as transfer


COUNTS = HERE / "rm511_low_weight_exact.json"
DIAGNOSTIC = HERE / "rm511_q1_low_randomstepconv_M22_d109_diagnostic.json"
OUTPUT = HERE / "rm511_q1_low_randomstepconv_M22_d109_outward.json"

BLOCK_BITS = 2048
OUTER_ROWS = 1024
MESSAGE_BITS = 1 << 20
MAXIMUM_KNOWN_WEIGHT = 124


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    ctx.prec = 192
    receipt = json.loads(COUNTS.read_text(encoding="utf-8"))
    if receipt.get("status") != "EXACT_INTEGER_FORMULA_CHECK":
        raise AssertionError("low-weight receipt is not accepting")
    counts = {
        int(weight): int(count)
        for weight, count in receipt["coefficients"].items()
    }
    expected_weights = {64, 96, 112, 120, 124}
    if set(counts) != expected_weights:
        raise AssertionError("the exact low-weight shell set changed")

    witness = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    if witness.get("status") != "BINARY64_WITNESS_DIAGNOSTIC":
        raise AssertionError("unexpected witness status")
    parameters = witness["parameters"]
    if int(parameters["inner_memory_bits"]) != transfer.MEMORY:
        raise AssertionError("witness memory differs from transfer memory")
    if int(parameters["distance_cutoff"]) != transfer.D:
        raise AssertionError("witness distance differs from transfer distance")

    rows = witness["weight_rows"]
    if {int(row["weight"]) for row in rows} != expected_weights:
        raise AssertionError("witness does not cover the exact low-weight shells")
    groups: dict[float, list[int]] = {}
    for row in rows:
        groups.setdefault(float(row["log_surprisal"]), []).append(int(row["weight"]))

    pointwise: dict[int, float] = {}
    for group_index, (log_surprisal, weights) in enumerate(sorted(groups.items())):
        transitions, surprisal = transfer.transition_matrices(
            log_surprisal, (transfer.Fraction(1),)
        )
        regions, region_exponents = transfer.coefficients_one_dimension(
            transitions[0], 0, transitions[1], 0, OUTER_ROWS, 1
        )
        coordinates, coordinate_exponents = transfer.coefficients_one_dimension(
            regions[0],
            int(region_exponents[0]),
            regions[1],
            int(region_exponents[1]),
            BLOCK_BITS,
            MAXIMUM_KNOWN_WEIGHT,
        )
        moments = transfer.row_sum_log2_upper(coordinates, coordinate_exponents)
        correction = transfer.correction_log2_upper(surprisal)
        for weight in weights:
            inner = min(0.0, transfer.up(float(moments[weight]) + correction))
            outer = transfer.log2_integer_upper(OUTER_ROWS * counts[weight])
            pointwise[weight] = transfer.up(outer + inner)
        print(f"witness_group,{group_index + 1},{len(groups)}", flush=True)

    total = sum(
        (
            arb(2) ** transfer.a(transfer.fraction_of_float(pointwise[weight]))
            for weight in sorted(pointwise)
        ),
        arb(0),
    )
    aggregate = transfer.upper_float(total.log() / arb(2).log())
    dominant_weight = max(pointwise, key=pointwise.__getitem__)
    payload = {
        "schema": "rm511-exact-low-q1-randomstepconv-outward-v1",
        "status": "OUTWARD_PARTIAL_Q1_CERTIFICATE",
        "parameters": {
            "message_bits": MESSAGE_BITS,
            "output_bits": transfer.N,
            "outer": "one fixed RM(5,11) [2048,1024,64] constituent repeated 1024 times",
            "outer_rows": OUTER_ROWS,
            "distance_cutoff": transfer.D,
            "relative_distance_lower": transfer.D / transfer.N,
            "inner_memory_bits": transfer.MEMORY,
        },
        "probability_space": {
            "outer": "fixed deterministic RM(5,11); no outer setup event",
            "routing": "independent uniform 2048-coordinate permutation per row and uniform 1024-position permutation per transposed region",
            "inner": "independent uniform binary 23-by-23 linear map at each of 2^21 RandomStepConv positions",
        },
        "claim": {
            "occupation": 1,
            "covered_outer_weights": sorted(pointwise),
            "covered_range": "every nonzero RM(5,11) shell below weight 128",
            "partial_q1_log2_upper": aggregate,
            "partial_q1_margin_bits_lower": -aggregate,
            "dominant_weight": dominant_weight,
            "dominant_pointwise_log2_upper": pointwise[dominant_weight],
            "covered_sum_below_2_minus_40": aggregate < -40.0,
        },
        "pointwise_log2_upper": {
            str(weight): pointwise[weight] for weight in sorted(pointwise)
        },
        "sources": [
            {"file": COUNTS.name, "sha256": sha256(COUNTS)},
            {"file": DIAGNOSTIC.name, "sha256": sha256(DIAGNOSTIC)},
        ],
        "scope": [
            "The bound sums every nonzero RM(5,11) weight shell below 128 with outward arithmetic.",
            "The diagnostic supplies only tilt witnesses; the certificate does not require them to be optimal.",
            "Weights 128 and above are not covered. This is not a complete occupation-one or minimum-distance certificate.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Certify the RM(4,9) occupation-one RandomStepConv-M22 contribution."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

from flint import arb, ctx


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
REPOSITORY = WORKSTREAM.parent.parent
sys.path.insert(0, str(WORKSTREAM))

import certify_single_random_constituent_sparse_outward as transfer


SPECTRUM = REPOSITORY / "scripts" / "rm512_256_spectrum.csv"
DIAGNOSTIC = HERE / "rm49_q1_randomstepconv_M22_d109_diagnostic.json"
OUTPUT = HERE / "rm49_q1_randomstepconv_M22_d109_outward.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_spectrum() -> list[int]:
    counts = [0] * (transfer.B + 1)
    with SPECTRUM.open(newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            counts[int(row["weight"])] = int(row["count"])
    if sum(counts) != 1 << transfer.K:
        raise AssertionError("RM spectrum mass is not 2^256")
    if counts[0] != 1 or counts[512] != 1:
        raise AssertionError("RM endpoints changed")
    if next(weight for weight in range(1, 513) if counts[weight]) != 32:
        raise AssertionError("RM minimum distance changed")
    return counts


def main() -> int:
    ctx.prec = 192
    counts = exact_spectrum()
    witness = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    if witness.get("status") != "BINARY64_WITNESS_DIAGNOSTIC":
        raise AssertionError("unexpected witness status")
    parameters = witness["parameters"]
    if int(parameters["memory_bits"]) != transfer.MEMORY:
        raise AssertionError("witness memory differs from transfer memory")
    if int(parameters["distance_cutoff"]) != transfer.D:
        raise AssertionError("witness distance differs from transfer distance")

    rows = witness["weight_rows"]
    observed = {int(row["weight"]) for row in rows}
    expected = {weight for weight in range(1, transfer.B + 1) if counts[weight]}
    if observed != expected:
        raise AssertionError("witness does not cover every nonzero RM shell")
    groups: dict[float, list[int]] = {}
    for row in rows:
        groups.setdefault(float(row["log_surprisal"]), []).append(int(row["weight"]))

    pointwise: dict[int, float] = {}
    for group_index, (log_surprisal, weights) in enumerate(sorted(groups.items())):
        transitions, surprisal = transfer.transition_matrices(log_surprisal, (transfer.Fraction(1),))
        regions, region_exponents = transfer.coefficients_one_dimension(
            transitions[0], 0, transitions[1], 0, transfer.L, 1
        )
        coordinates, coordinate_exponents = transfer.coefficients_one_dimension(
            regions[0], int(region_exponents[0]),
            regions[1], int(region_exponents[1]),
            transfer.B, transfer.B,
        )
        moments = transfer.row_sum_log2_upper(coordinates, coordinate_exponents)
        correction = transfer.correction_log2_upper(surprisal)
        for weight in weights:
            inner = min(0.0, transfer.up(float(moments[weight]) + correction))
            outer = transfer.log2_integer_upper(transfer.L * counts[weight])
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
        "schema": "rm49-exact-q1-randomstepconv-outward-v1",
        "status": "OUTWARD_Q1_CERTIFICATE",
        "parameters": {
            "message_bits": transfer.K * transfer.L,
            "output_bits": transfer.N,
            "outer": "one fixed RM(4,9) [512,256,32] constituent repeated 4096 times",
            "outer_rows": transfer.L,
            "distance_cutoff": transfer.D,
            "relative_distance_lower": transfer.D / transfer.N,
            "inner_memory_bits": transfer.MEMORY,
        },
        "probability_space": {
            "outer": "fixed deterministic RM(4,9); no outer setup event",
            "routing": "independent uniform 512-coordinate permutation per row and uniform 4096-position permutation per transposed region",
            "inner": "independent uniform binary 23-by-23 linear map at each of 2^21 RandomStepConv positions",
        },
        "claim": {
            "occupation": 1,
            "q1_log2_upper": aggregate,
            "q1_margin_bits_lower": -aggregate,
            "dominant_weight": dominant_weight,
            "dominant_pointwise_log2_upper": pointwise[dominant_weight],
            "closes_40_bits": aggregate < -40.0,
        },
        "sources": [
            {"file": str(SPECTRUM.relative_to(REPOSITORY)).replace("\\", "/"), "sha256": sha256(SPECTRUM)},
            {"file": DIAGNOSTIC.name, "sha256": sha256(DIAGNOSTIC)},
        ],
        "scope": [
            "The bound sums every nonzero RM(4,9) weight shell with outward arithmetic.",
            "The diagnostic supplies only tilt witnesses; the certificate does not require them to be optimal.",
            "Only occupation one is certified. No full minimum-distance claim follows from this receipt alone.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

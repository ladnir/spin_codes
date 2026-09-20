#!/usr/bin/env python3
"""Combine the outward sparse and dense one-shot constituent receipts."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

from flint import arb, ctx

from certify_single_random_constituent_dense_outward import (
    B,
    D,
    K,
    L,
    MEMORY,
    N,
    a,
    exact_caps,
    fraction_of_float,
    upper_float,
)


WORKSTREAM = Path(__file__).resolve().parent
SPARSE = WORKSTREAM / "single_random_constituent_B512_sparse_outward_s22.json"
DENSE = WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_cover_outward_s22.json"
OUTPUT = WORKSTREAM / "single_random_constituent_B512_combined_outward_s22.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exp2_of_upper_float(value: float) -> arb:
    return arb(2) ** a(fraction_of_float(value))


def main() -> None:
    ctx.prec = 192
    sparse = json.loads(SPARSE.read_text(encoding="utf-8"))
    dense = json.loads(DENSE.read_text(encoding="utf-8"))
    if sparse.get("status") != "OUTWARD_SPARSE_CERTIFICATE":
        raise AssertionError("sparse receipt is not accepting")
    if dense.get("status") != "OUTWARD_DENSE_CERTIFICATE":
        raise AssertionError("dense receipt is not accepting")
    if sparse["claim"]["maximum_occupation"] + 1 != dense["claim"]["minimum_occupation"]:
        raise AssertionError("occupation cover has a gap or overlap")
    if dense["claim"]["maximum_occupation"] != L:
        raise AssertionError("occupation cover does not reach all outer rows")

    sparse_log2 = float(sparse["claim"]["aggregate_log2_upper"])
    dense_log2 = float(dense["claim"]["aggregate_log2_upper"])
    conditional_mass = exp2_of_upper_float(sparse_log2) + exp2_of_upper_float(dense_log2)
    conditional_log2 = upper_float(conditional_mass.log() / arb(2).log())

    _, setup_failure = exact_caps()
    total_failure = a(setup_failure) + conditional_mass
    total_log2 = upper_float(total_failure.log() / arb(2).log())
    if total_log2 >= -40:
        raise AssertionError("combined one-shot failure does not close at 40 bits")

    payload = {
        "schema": "single-random-constituent-combined-outward-v1",
        "status": "OUTWARD_FINITE_DISTANCE_CERTIFICATE",
        "claim": {
            "minimum_distance_at_least": D,
            "relative_distance_lower": D / N,
            "conditional_bad_probability_log2_upper": conditional_log2,
            "conditional_margin_bits_lower": -conditional_log2,
            "constituent_event_failure_upper": upper_float(a(setup_failure)),
            "constituent_event_margin_bits_lower": -upper_float(a(setup_failure).log() / arb(2).log()),
            "unconditional_setup_or_bad_log2_upper": total_log2,
            "unconditional_margin_bits_lower": -total_log2,
            "comparison_to_2^-40": True,
        },
        "parameters": {
            "message_bits": K * L,
            "output_bits": N,
            "outer_block_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "distance_cutoff": D,
            "memory_bits": MEMORY,
        },
        "probability_space": {
            "outer": (
                "sample one uniform binary 256-by-512 matrix and reuse it in all 4096 outer rows; no rejection test is performed"
            ),
            "routing": (
                "independently sample one uniform 512-coordinate permutation per outer row and one uniform 4096-position permutation per transposed region"
            ),
            "inner": (
                "independently sample one uniform binary 23-by-23 linear map at each of the 2^21 RandomStepConv positions; all maps are fixed and shared by every message"
            ),
        },
        "proof_accounting": {
            "sparse_occupations": [1, 159],
            "dense_occupations": [160, L],
            "outer_event": (
                "full row rank and the integer shell caps generated with per-shell target 2^-51"
            ),
            "failure_union": (
                "Pr[outer event fails] plus the conditional expected number of nonzero messages below the distance cutoff"
            ),
        },
        "component_receipts": [
            {"file": SPARSE.name, "sha256": sha256(SPARSE)},
            {"file": DENSE.name, "sha256": sha256(DENSE)},
        ],
        "limitations": [
            "The theorem concerns the RandomStepConv-M22 comparison ensemble, not RM2Sub.",
            "The information-theoretic inner samples a fresh 23-by-23 binary matrix at every position; no compact implementation claim is made.",
            "The theorem is a finite k=2^20 result. An asymptotic family requires a separate parameter schedule.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()

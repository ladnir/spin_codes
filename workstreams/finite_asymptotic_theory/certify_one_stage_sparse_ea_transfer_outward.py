#!/usr/bin/env python3
"""Attach one-stage sparse-EA caps to the 10.9% SPIN transfer.

The verifier recomputes occupations one and two with the realized-spectrum
caps.  For occupations at least three, it checks that each new band majorant
is no larger than the majorant used by the frozen outward certificates.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

import certify_single_random_constituent_sparse_outward as sparse
from certify_single_random_constituent_dense_outward import (
    B,
    D,
    K,
    L,
    MEMORY,
    N,
    P_CENTRAL,
    P_DEFECT,
    a,
    exact_band_majorant,
    exact_caps,
    fraction_of_float,
    upper_float,
)


WORKSTREAM = Path(__file__).resolve().parent
CAPS = WORKSTREAM / "one_stage_sparse_ea_K256_B512_r33_spectrum_caps_outward.json"
FROZEN_SPARSE = WORKSTREAM / "single_random_constituent_B512_sparse_outward_s22.json"
FROZEN_DENSE = WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_cover_outward_s22.json"
OUTPUT = WORKSTREAM / "one_stage_sparse_ea_K256_B512_r33_randomstepconv_M22_d109_outward.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exp2_upper(log2_value: float) -> arb:
    return arb(2) ** a(fraction_of_float(log2_value))


def main() -> None:
    ctx.prec = 192
    cap_payload = json.loads(CAPS.read_text(encoding="utf-8"))
    if cap_payload.get("status") != "OUTWARD_CAP_CERTIFICATE":
        raise AssertionError("cap receipt is not accepting")
    caps = [int(value) for value in cap_payload["caps"]]
    if len(caps) != B + 1:
        raise AssertionError("cap vector has the wrong length")
    if any(caps[weight] for weight in [0, *range(1, 42), *range(471, B + 1)]):
        raise AssertionError("unexpected cap outside weights 42 through 470")

    frozen_caps, _ = exact_caps()
    new_majorants = {
        "low": exact_band_majorant(caps, 42, 79, P_DEFECT),
        "central": exact_band_majorant(caps, 80, 432, P_CENTRAL),
        "high": exact_band_majorant(caps, 433, 470, 1 - P_DEFECT),
    }
    frozen_majorants = {
        "low": exact_band_majorant(frozen_caps, 42, 79, P_DEFECT),
        "central": exact_band_majorant(frozen_caps, 80, 432, P_CENTRAL),
        "high": exact_band_majorant(frozen_caps, 433, 470, 1 - P_DEFECT),
    }
    domination = {
        name: new_majorants[name] <= frozen_majorants[name]
        for name in new_majorants
    }
    if not all(domination.values()):
        raise AssertionError("new caps exceed a frozen band majorant")

    q1, q1_details = sparse.q1_bound(caps)
    q2, q2_details = sparse.q2_bound(caps)
    frozen_sparse = json.loads(FROZEN_SPARSE.read_text(encoding="utf-8"))
    frozen_dense = json.loads(FROZEN_DENSE.read_text(encoding="utf-8"))
    if frozen_sparse.get("status") != "OUTWARD_SPARSE_CERTIFICATE":
        raise AssertionError("frozen sparse receipt is not accepting")
    if frozen_dense.get("status") != "OUTWARD_DENSE_CERTIFICATE":
        raise AssertionError("frozen dense receipt is not accepting")
    q3_159 = float(frozen_sparse["claim"]["q3_159_log2_upper"])
    q160_L = float(frozen_dense["claim"]["aggregate_log2_upper"])
    conditional_terms = [q1, q2, q3_159, q160_L]
    conditional_mass = sum(
        (exp2_upper(value) for value in conditional_terms), arb(0)
    )
    conditional_log2 = upper_float(conditional_mass.log() / arb(2).log())

    setup_log2 = float(cap_payload["claim"]["total_setup_failure_upper_log2"])
    setup_failure = exp2_upper(setup_log2)
    total_failure = setup_failure + conditional_mass
    total_log2 = upper_float(total_failure.log() / arb(2).log())
    print(
        json.dumps(
            {
                "q1_log2_upper": q1,
                "q2_log2_upper": q2,
                "q3_159_log2_upper": q3_159,
                "q160_4096_log2_upper": q160_L,
                "conditional_log2_upper": conditional_log2,
                "setup_log2_upper": setup_log2,
                "total_log2_upper": total_log2,
            },
            indent=2,
        ),
        flush=True,
    )
    if total_log2 >= -40:
        raise AssertionError("combined finite distance bound misses 40 bits")

    payload = {
        "schema": "one-stage-sparse-ea-randomstepconv-finite-outward-v1",
        "status": "OUTWARD_FINITE_DISTANCE_CERTIFICATE",
        "parameters": {
            "message_bits": K * L,
            "output_bits": N,
            "outer_block_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "sparse_right_degree": 33,
            "distance_cutoff": D,
            "relative_distance_lower": D / N,
            "inner_memory_bits": MEMORY,
            "rank_test_attempts": 16,
        },
        "probability_space": {
            "outer_attempt": "sample 512 independent uniform weight-33 row vectors in F_2^256, then apply the zero-initialized accumulator",
            "outer_setup": "sample independent attempts, accept the first full-rank constituent, and abort after 16 failures",
            "outer_reuse": "reuse the accepted constituent in all 4096 outer rows",
            "routing": "independently sample one uniform 512-coordinate permutation per outer row and one uniform 4096-position permutation per transposed region",
            "inner": "independently sample one uniform binary 23-by-23 linear map at each of the 2^21 RandomStepConv positions; fix all maps for the code",
        },
        "claim": {
            "minimum_distance_at_least": D,
            "relative_distance_lower": D / N,
            "q1_log2_upper": q1,
            "q2_log2_upper": q2,
            "q3_159_log2_upper": q3_159,
            "q160_4096_log2_upper": q160_L,
            "conditional_bad_probability_log2_upper": conditional_log2,
            "conditional_margin_bits_lower": -conditional_log2,
            "setup_failure_log2_upper": setup_log2,
            "setup_margin_bits_lower": -setup_log2,
            "setup_or_bad_log2_upper": total_log2,
            "overall_margin_bits_lower": -total_log2,
            "comparison_to_2^-40": True,
        },
        "transfer_checks": {
            "band_majorant_domination": domination,
            "q1_recomputed_with_exact_caps": q1_details,
            "q2_recomputed_with_exact_caps": q2_details,
            "q3_159_reused_by_band_monotonicity": True,
            "q160_4096_reused_by_band_monotonicity": True,
        },
        "sources": [
            {"file": path.name, "sha256": sha256(path)}
            for path in (CAPS, FROZEN_SPARSE, FROZEN_DENSE)
        ],
        "scope": [
            "The distance statement is unconditional over bounded rank-tested outer setup, routing, and RandomStepConv setup.",
            "The certificate applies to RandomStepConv-M22, not RM2Sub.",
            "The routing includes independent uniform row-coordinate and region permutations; no structured-permutation replacement is proved here.",
            "The encoder reuses one accepted outer constituent; it does not resample a constituent per outer row.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()

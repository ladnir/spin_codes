#!/usr/bin/env python3
"""Combine the conditional BA-128 cap event with the certified SPIN transfer.

This verifier does not prove the accumulator singular-value lemma.  Subject
to that lemma, it checks the BA cap-event receipt, reuses the already outward-
certified cap-conditional sparse and dense RandomStepConv-M22 receipts, and
combines their failure probabilities with outward Arb arithmetic.
"""

from __future__ import annotations

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
EVENT = WORKSTREAM / "ba128_pair_singular_random_cap_event_outward_B512.json"
SPARSE = WORKSTREAM / "single_random_constituent_B512_sparse_outward_s22.json"
DENSE = WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_cover_outward_s22.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba128_randomstepconv_combined_outward_s22.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exp2_of_upper_float(value: float) -> arb:
    return arb(2) ** a(fraction_of_float(value))


def main() -> None:
    ctx.prec = 256
    event = json.loads(EVENT.read_text(encoding="utf-8"))
    sparse = json.loads(SPARSE.read_text(encoding="utf-8"))
    dense = json.loads(DENSE.read_text(encoding="utf-8"))
    if event.get("status") != "CONDITIONAL_OUTWARD_CERTIFICATE":
        raise AssertionError("BA cap-event receipt is not conditionally accepting")
    if event["parameters"]["cap_source"] != "certify_single_random_constituent_dense_outward.exact_caps":
        raise AssertionError("BA event uses an unexpected cap vector")
    if event["parameters"]["accumulator_stages"] != 128:
        raise AssertionError("BA event uses an unexpected stage count")
    if sparse.get("status") != "OUTWARD_SPARSE_CERTIFICATE":
        raise AssertionError("sparse conditional receipt is not accepting")
    if dense.get("status") != "OUTWARD_DENSE_CERTIFICATE":
        raise AssertionError("dense conditional receipt is not accepting")
    if sparse["claim"]["maximum_occupation"] + 1 != dense["claim"]["minimum_occupation"]:
        raise AssertionError("occupation cover has a gap or overlap")
    if dense["claim"]["maximum_occupation"] != L:
        raise AssertionError("occupation cover does not reach all outer rows")

    caps, _ = exact_caps()
    cap_bytes = ",".join(str(int(cap)) for cap in caps).encode("ascii")
    if hashlib.sha256(cap_bytes).hexdigest() != event["parameters"]["cap_vector_sha256"]:
        raise AssertionError("BA event cap-vector hash mismatch")

    # Every imported float is already an outward upper endpoint.  Converting
    # its exact binary64 value to a rational therefore preserves the bound.
    event_failure = a(
        fraction_of_float(float(event["claim"]["cap_event_failure_upper"]))
    )
    sparse_mass = exp2_of_upper_float(float(sparse["claim"]["aggregate_log2_upper"]))
    dense_mass = exp2_of_upper_float(float(dense["claim"]["aggregate_log2_upper"]))
    conditional_mass = sparse_mass + dense_mass
    total_failure = event_failure + conditional_mass
    conditional_log2 = upper_float(conditional_mass.log() / arb(2).log())
    total_log2 = upper_float(total_failure.log() / arb(2).log())
    if total_log2 >= -40:
        raise AssertionError("conditional BA-128 theorem does not close at 40 bits")

    payload = {
        "schema": "ebch128x4-ba128-randomstepconv-combined-v1",
        "status": "CONDITIONAL_OUTWARD_FINITE_DISTANCE_CERTIFICATE",
        "claim": {
            "minimum_distance_at_least": D,
            "relative_distance_lower": D / N,
            "cap_conditional_bad_probability_log2_upper": conditional_log2,
            "cap_conditional_margin_bits_lower": -conditional_log2,
            "constituent_event_failure_upper": upper_float(event_failure),
            "constituent_event_margin_bits_lower": -upper_float(event_failure.log() / arb(2).log()),
            "total_failure_log2_upper": total_log2,
            "total_margin_bits_lower": -total_log2,
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
            "accumulator_stages": 128,
            "singular_value_upper_hypothesis": "63/1000",
        },
        "outer_constituent": {
            "base": "direct sum of four copies of the fixed extended BCH [128,64,22] code",
            "stage": "independent uniform permutation of 512 coordinates followed by the zero-state binary prefix accumulator",
            "reuse": "the 128 sampled stage permutations define one constituent, repeated in all 4096 outer rows",
        },
        "probability_space": {
            "outer": "128 independent uniform elements of S_512 sampled once",
            "routing": "independently sample one uniform 512-coordinate permutation per outer row and one uniform 4096-position permutation per transposed region",
            "inner": "independently sample one uniform binary 23-by-23 linear map at each of the 2^21 RandomStepConv positions; all maps are fixed and shared by every message",
        },
        "hypothesis": (
            "in stationary L2, the second singular values of the one-word and ordered-distinct-pair type operators for one permuted accumulator are each at most 63/1000"
        ),
        "proof_accounting": {
            "outer_event": "the exact integer cap vector used by the certified random-[512,256] transfer",
            "conditional_transfer": "the existing outward sparse and dense receipts depend only on that cap vector",
            "failure_union": "BA-128 cap-event failure plus the cap-conditional expected number of low-weight nonzero messages",
        },
        "component_receipts": [
            {"file": EVENT.name, "sha256": sha256(EVENT), "role": "conditional BA-128 cap event"},
            {"file": SPARSE.name, "sha256": sha256(SPARSE), "role": "cap-conditional occupations 1 through 159"},
            {"file": DENSE.name, "sha256": sha256(DENSE), "role": "cap-conditional occupations 160 through 4096"},
        ],
        "limitations": [
            "The accumulator singular-value hypothesis is an open proof obligation; this is not yet an unconditional construction certificate.",
            "The theorem concerns RandomStepConv-M22, not RM2Sub.",
            "The 128 accumulator stages are a proof fallback, not a performance-competitive design.",
            "The inner samples a fresh 23-by-23 binary matrix at every position; no compact inner implementation claim is made.",
            "The theorem is finite at k=2^20; it does not by itself provide an asymptotic schedule.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()

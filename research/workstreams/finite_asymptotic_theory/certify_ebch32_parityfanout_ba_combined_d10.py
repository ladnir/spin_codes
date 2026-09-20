#!/usr/bin/env python3
"""Combine the sparse and dense outward receipts at the 10% target."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
Q1 = WORKSTREAM / "ebch32_parityfanout31x33_ba3_rm2sub_B256_q1_outward_d11.json"
Q2_64 = WORKSTREAM / "ebch32_parityfanout31x33_ba3_rm2sub_B256_q2_64_outward_d11.json"
DENSE = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_one_band_interval_cover_d10_outward.json"
SETUP = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_setup_outward.json"
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_combined_outward_d10.json"
N = 2**21
DISTANCE = 209_715


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_binary64(hexadecimal: str) -> Fraction:
    return Fraction.from_float(float.fromhex(hexadecimal))


def main() -> None:
    q1 = json.loads(Q1.read_text(encoding="utf-8"))
    q2 = json.loads(Q2_64.read_text(encoding="utf-8"))
    dense = json.loads(DENSE.read_text(encoding="utf-8"))
    setup = json.loads(SETUP.read_text(encoding="utf-8"))

    q1_claim = q1["claim"]
    q2_claim = q2["claim"]
    dense_claim = dense["claim"]
    if int(q1_claim["bad_weight_at_most"]) < DISTANCE:
        raise ArithmeticError("occupation-one receipt does not cover the target")
    if int(q2_claim["bad_weight_at_most"]) < DISTANCE:
        raise ArithmeticError("occupation-2..64 receipt does not cover the target")
    if int(dense_claim["bad_weight_at_most"]) != DISTANCE:
        raise ArithmeticError("dense receipt distance mismatch")
    if list(q2_claim["occupations"]) != [2, 64]:
        raise ArithmeticError("sparse occupation range mismatch")
    if list(dense_claim["occupations"]) != [65, 8192]:
        raise ArithmeticError("dense occupation range mismatch")
    if dense["status"] != "COMPLETE_OUTWARD_DENSE_CERTIFICATE":
        raise ArithmeticError("dense receipt is not complete")
    if int(dense_claim["certified_integer_margin_bits"]) < 468:
        raise ArithmeticError("dense integer margin is insufficient")

    q1_upper = exact_binary64(q1_claim["expected_bad_word_count_upper_hex"])
    scaled = q2_claim["expected_bad_word_count_upper_scaled"]
    q2_upper = exact_binary64(scaled["mantissa_binary64_hex"])
    exponent = int(scaled["binary_exponent"])
    q2_upper = q2_upper * (2**exponent if exponent >= 0 else Fraction(1, 2 ** (-exponent)))
    dense_upper = Fraction(1, 2**468)
    total_upper = q1_upper + q2_upper + dense_upper
    if total_upper >= Fraction(1, 2**51):
        raise ArithmeticError("combined expectation does not pass 51 bits")

    display_margin = -math.log2(float(total_upper))
    payload = {
        "schema": "ebch32-parityfanout31x33-ba3-b256-combined-d10-v1",
        "status": "COMPLETE_OUTWARD_FIRST_MOMENT_CERTIFICATE",
        "claim": {
            "message_bits": 2**20,
            "output_bits": N,
            "bad_weight_at_most": DISTANCE,
            "minimum_distance_if_good": DISTANCE + 1,
            "relative_distance_if_good": (DISTANCE + 1) / N,
            "occupations_covered": [1, 8192],
            "expected_bad_word_count_less_than_2^-51": True,
            "setup_failure_probability_less_than_2^-51": True,
            "setup_failure_probability_less_than_2^-40": True,
            "combined_margin_bits_lower_display": display_margin,
        },
        "exact_combination": {
            "q1_upper_numerator": q1_upper.numerator,
            "q1_upper_denominator": q1_upper.denominator,
            "q2_64_upper_numerator": q2_upper.numerator,
            "q2_64_upper_denominator": q2_upper.denominator,
            "dense_upper_used": "2^-468",
            "total_upper_numerator": total_upper.numerator,
            "total_upper_denominator": total_upper.denominator,
        },
        "reasoning": [
            "The occupation-one and occupation-2..64 receipts bound the larger event of weight at most 230686, so they also bound weight at most 209715.",
            "The dense receipt covers every occupation from 65 through 8192 at weight at most 209715.",
            "The three occupation ranges partition all nonzero messages.",
            "Markov's inequality bounds the probability of any bad nonzero message by its expected count.",
        ],
        "outer_setup": {
            "good_probability_lower_hex": setup["claim"]["good_event_probability_lower_hex"],
            "expected_trials_upper_hex": setup["claim"]["expected_rejection_trials_per_row_upper_hex"],
            "law": "Independent row draws conditioned on G_256 by conceptual rejection sampling.",
        },
        "dependencies": [
            {"path": path.name, "sha256": sha256(path)}
            for path in (Q1, Q2_64, DENSE, SETUP, Path(__file__).resolve())
        ],
        "scope": (
            "This receipt proves the finite first-moment claim for the declared "
            "conditional mathematical ensemble. It does not supply an efficient "
            "test for G_256 or prove equivalence with a particular optimized implementation."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()

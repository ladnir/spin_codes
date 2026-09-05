#!/usr/bin/env python3
"""Bind bounded exact rejection setup to the finite 10% distance receipt."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
SETUP_RECEIPT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_setup_outward.json"
DISTANCE_RECEIPT = (
    WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_combined_outward_d10.json"
)
OUTPUT = (
    WORKSTREAM
    / "ebch32_parityfanout31x33_ba3_B256_bounded_setup_combined_outward_d10.json"
)

ROW_COUNT = 8192
ATTEMPTS_PER_ROW = 6


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ratio_payload(value: Fraction) -> dict[str, object]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "display_log2_upper": math.log2(float(value)),
        "display_margin_bits": -math.log2(float(value)),
    }


def certify() -> dict[str, object]:
    setup = json.loads(SETUP_RECEIPT.read_text(encoding="utf-8"))
    distance = json.loads(DISTANCE_RECEIPT.read_text(encoding="utf-8"))
    if setup["status"] != "OUTWARD_CERTIFICATE":
        raise ArithmeticError("setup receipt is not an outward certificate")
    if setup["claim"]["permitted_final_nonzero_weights"] != [24, 232]:
        raise ArithmeticError("setup acceptance window mismatch")
    if setup["parameters"]["row_length"] != 256:
        raise ArithmeticError("setup row length mismatch")
    if distance["status"] != "COMPLETE_OUTWARD_FIRST_MOMENT_CERTIFICATE":
        raise ArithmeticError("distance receipt is not complete")
    claim = distance["claim"]
    if (
        claim["output_bits"] != 1 << 21
        or claim["bad_weight_at_most"] != 209_715
        or claim["occupations_covered"] != [1, ROW_COUNT]
    ):
        raise ArithmeticError("distance claim mismatch")

    # The stored hexadecimal binary64 value is itself an outward upper bound.
    # Convert that exact dyadic value to a rational before exponentiating it.
    rejected_draw_upper = Fraction.from_float(
        float.fromhex(setup["claim"]["expected_rejected_word_count_upper_hex"])
    )
    if not 0 < rejected_draw_upper < 1:
        raise ArithmeticError("invalid rejected-draw probability upper bound")

    exact = distance["exact_combination"]
    conditional_bad_upper = Fraction(
        int(exact["total_upper_numerator"]),
        int(exact["total_upper_denominator"]),
    )
    if not conditional_bad_upper < Fraction(1, 1 << 51):
        raise ArithmeticError("conditional distance margin is below 51 bits")

    per_row_abort_upper = rejected_draw_upper**ATTEMPTS_PER_ROW
    any_abort_upper = ROW_COUNT * per_row_abort_upper
    total_failure_upper = any_abort_upper + conditional_bad_upper
    if not total_failure_upper < Fraction(1, 1 << 45):
        raise ArithmeticError("bounded setup and distance do not retain 45 bits")

    return {
        "schema": "ebch32-parityfanout31x33-ba3-b256-bounded-setup-d10-v1",
        "status": "COMPLETE_BOUNDED_SETUP_AND_DISTANCE_CERTIFICATE",
        "claim": {
            "message_bits": 1 << 20,
            "output_bits": 1 << 21,
            "bad_weight_at_most": 209_715,
            "minimum_distance_if_successful_and_good": 209_716,
            "attempts_per_row": ATTEMPTS_PER_ROW,
            "probability_of_setup_abort_or_bad_distance_less_than_2^-45": True,
            "probability_of_setup_abort_or_bad_distance_less_than_2^-40": True,
            "combined_margin_bits_lower_display": -math.log2(
                float(total_failure_upper)
            ),
        },
        "exact_bounds": {
            "rejected_single_draw_upper": ratio_payload(rejected_draw_upper),
            "per_row_abort_upper": ratio_payload(per_row_abort_upper),
            "any_row_abort_union_upper": ratio_payload(any_abort_upper),
            "conditional_bad_distance_upper": ratio_payload(conditional_bad_upper),
            "abort_or_bad_distance_upper": ratio_payload(total_failure_upper),
        },
        "setup_algorithm": {
            "row_count": ROW_COUNT,
            "attempts_per_row": ATTEMPTS_PER_ROW,
            "draw_law": (
                "For each row and attempt, independently sample ParityFanout-31x33 "
                "and two uniform accumulator interleavers."
            ),
            "exact_acceptance_test": (
                "Enumerate all 2^128 messages of the candidate linear row code; "
                "accept exactly when every nonzero output has weight in [24,232]."
            ),
            "selection": "Use the first accepted draw for each row.",
            "failure": (
                "Return setup failure as soon as a row has no accepted draw in six attempts."
            ),
        },
        "reasoning": [
            "Markov's inequality bounds rejection of one unconditioned row draw by the certified expected number of forbidden nonzero words.",
            "Six independent rejected draws are required for one row to abort.",
            "A union bound over 8192 rows gives the displayed setup-abort bound.",
            "Conditioned on no abort, first-success selection gives independent rows with exactly the conditional law used by the distance certificate.",
            "The unconditional probability of setup abort or a bad accepted code is at most the sum of the abort bound and the conditional distance bound.",
        ],
        "complexity_scope": {
            "online_encoding": "O(N) bit operations for the fixed finite construction.",
            "setup": (
                "Finite and exact but impractical: at most 6*8192*2^128 row-message "
                "evaluations. This receipt does not claim linear-time setup."
            ),
        },
        "dependencies": [
            {"path": SETUP_RECEIPT.name, "sha256": sha256(SETUP_RECEIPT)},
            {"path": DISTANCE_RECEIPT.name, "sha256": sha256(DISTANCE_RECEIPT)},
            {"path": Path(__file__).name, "sha256": sha256(Path(__file__))},
        ],
    }


def main() -> None:
    payload = certify()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()

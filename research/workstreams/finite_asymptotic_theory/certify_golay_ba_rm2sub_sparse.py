#!/usr/bin/env python3
"""Transfer the random-outer sparse certificates to the selected BA outer."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from mpmath import iv

from certify_golay_ba_rm2sub_joint_interval import as_interval, load_segments


iv.dps = 100
LIKELIHOOD_EXCESS = Fraction(1281, 100000)
BLOCK_CONSTANT = Fraction(39, 4)


def upper_float(value) -> float:
    return math.nextafter(float(value.b), math.inf)


def fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": format(float(value), ".18g"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--majorant",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "golay_ba3_concave_majorant.json"
        ),
    )
    parser.add_argument(
        "--small-receipt",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "rm2sub_dense_small_exact_d11.json"
        ),
    )
    parser.add_argument(
        "--fixed-receipt",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "rm2sub_uniform_fixed_occupation_d11.json"
        ),
    )
    parser.add_argument(
        "--weight-coupled-fixed-receipt",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "golay_ba3_rm2sub_weight_coupled_fixed_d11.json"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "golay_ba3_rm2sub_sparse_d11.json"
        ),
    )
    args = parser.parse_args()

    segments = load_segments(args.majorant)
    ln2 = iv.log(iv.mpf(2))
    endpoint_checks = []
    maximum_upper = -math.inf
    maximum_claim_margin_upper = -math.inf
    claimed_excess = as_interval(LIKELIHOOD_EXCESS)
    for segment in segments:
        for endpoint in (segment.lower, segment.upper):
            x = as_interval(endpoint)
            entropy = -x * iv.log(x) - (1 - x) * iv.log(1 - x)
            majorant = (
                as_interval(segment.slope) * x
                + as_interval(segment.intercept)
            )
            excess = majorant - entropy + ln2 / 2
            upper = upper_float(excess)
            claim_margin_upper = upper_float(excess - claimed_excess)
            maximum_upper = max(maximum_upper, upper)
            maximum_claim_margin_upper = max(
                maximum_claim_margin_upper, claim_margin_upper
            )
            endpoint_checks.append(
                {
                    "segment": segment.index,
                    "weight": str(endpoint),
                    "likelihood_excess_upper_natural": upper,
                }
            )
    likelihood_proved = maximum_claim_margin_upper < 0.0

    small = json.loads(args.small_receipt.read_text(encoding="utf-8"))
    if not small.get("all_checks_pass"):
        raise ValueError("the four-state small-density receipt did not pass")
    coefficient = small["log_bound"]["coefficient_upper_bound_natural"]
    random_small_decay = -Fraction(
        int(coefficient["numerator"]), int(coefficient["denominator"])
    )
    structured_small_decay = random_small_decay - LIKELIHOOD_EXCESS

    fixed = json.loads(args.fixed_receipt.read_text(encoding="utf-8"))
    if not all(fixed["assertions"].values()):
        raise ValueError("the fixed-occupation receipt did not pass")
    gap_record = fixed["exact_bounds"]["gap_bits_per_active_row_lower"]
    random_fixed_gap = Fraction(
        int(gap_record["numerator"]), int(gap_record["denominator"])
    )
    log_two_record = fixed["exact_bounds"]["log_two_lower"]
    log_two_lower = Fraction(
        int(log_two_record["numerator"]),
        int(log_two_record["denominator"]),
    )
    structured_fixed_gap = (
        random_fixed_gap - LIKELIHOOD_EXCESS / log_two_lower
    )
    weight_coupled = json.loads(
        args.weight_coupled_fixed_receipt.read_text(encoding="utf-8")
    )
    if weight_coupled.get("status") != "proved" or not all(
        weight_coupled["assertions"].values()
    ):
        raise ValueError("the weight-coupled fixed receipt did not pass")
    selected_constant = weight_coupled["witness"]["block_constant"]
    receipt_constant = Fraction(
        int(selected_constant["numerator"]),
        int(selected_constant["denominator"]),
    )
    if receipt_constant != BLOCK_CONSTANT:
        raise ValueError("the fixed receipt uses a different block constant")
    fixed_schedule_decay = float(
        weight_coupled["bounds"][
            "certified_schedule_decay_lower_natural_per_active_row"
        ]
    )

    # The growing sparse support costs at most ln(2)/c per active row and
    # output bit.  The exact elementary bound ln(2)<7/10 suffices here.
    growing_support_margin = (
        structured_small_decay - Fraction(7, 10) / BLOCK_CONSTANT
    )

    # The Q=1,2 certificates accept c=18/5 before transfer.  Hence their
    # certified random per-row gaps exceed 5/18.  Subtract the same uniform
    # likelihood excess and test c=39/4.
    one_two_gap = Fraction(5, 18) - LIKELIHOOD_EXCESS / log_two_lower
    one_two_schedule_decay = BLOCK_CONSTANT * one_two_gap - 1

    assertions = {
        "majorant_likelihood_excess_below_1281_over_100000": likelihood_proved,
        "structured_small_decay_positive": structured_small_decay > 0,
        "growing_sparse_support_margin_positive": growing_support_margin > 0,
        "fixed_Q_ge_3_schedule_decay_positive": fixed_schedule_decay > 0.0,
        "fixed_Q_1_2_schedule_decay_positive": one_two_schedule_decay > 0,
    }
    if not all(assertions.values()):
        raise AssertionError(f"sparse transfer assertion failed: {assertions}")

    payload = {
        "schema": "golay-ba3-rm2sub-sparse-transfer-v2",
        "status": "proved",
        "claim": (
            "With B=(39/4) log2(N)+O(1), the weight-coupled receipt covers "
            "every fixed Q>=3.  The BA likelihood excess below 0.01281 "
            "transfers the Q=1,2 and growing-sparse RM2Sub bounds."
        ),
        "selected_block_constant": fraction_record(BLOCK_CONSTANT),
        "likelihood_excess": {
            "claimed_upper": fraction_record(LIKELIHOOD_EXCESS),
            "largest_outward_endpoint_upper": maximum_upper,
            "largest_outward_claim_margin_upper": (
                maximum_claim_margin_upper
            ),
            "convexity_rule": (
                "On each affine majorant segment, majorant(x)-h(x)+ln(2)/2 "
                "is convex, so its maximum occurs at an endpoint."
            ),
            "endpoint_checks": endpoint_checks,
        },
        "exact_bounds": {
            "random_small_decay_natural": fraction_record(random_small_decay),
            "structured_small_decay_natural": fraction_record(
                structured_small_decay
            ),
            "growing_sparse_support_margin_natural": fraction_record(
                growing_support_margin
            ),
            "random_fixed_gap_bits_per_active_row": fraction_record(
                random_fixed_gap
            ),
            "structured_fixed_gap_bits_per_active_row": fraction_record(
                structured_fixed_gap
            ),
            "fixed_Q_ge_3_weight_coupled_schedule_decay_natural": (
                fixed_schedule_decay
            ),
            "fixed_Q_1_2_gap_bits_per_active_row": fraction_record(one_two_gap),
            "fixed_Q_1_2_schedule_decay": fraction_record(
                one_two_schedule_decay
            ),
        },
        "assertions": assertions,
        "scope": {
            "fixed_occupations": "every fixed Q>=1",
            "growing_sparse_occupations": "Q->infinity and Q/L<=10^-4",
            "outer_selection_factor": (
                "The B^2 factor per active row is exp(o(QB)) for growing Q "
                "and is absorbed into o_Q(B) for each fixed Q."
            ),
        },
    }
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                key: value
                for key, value in payload.items()
                if key != "likelihood_excess"
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

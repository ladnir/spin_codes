#!/usr/bin/env python3
"""Outward certificate for the weight-coupled fixed-occupation bound."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from mpmath import iv

from certify_golay_ba_rm2sub_joint_interval import as_interval, load_segments


iv.dps = 100
ZERO = Fraction(0)
ONE = Fraction(1)
DELTA = Fraction(11, 100)
RESET_PROBABILITY = Fraction(1, (1 << 19) - 1)
LIVE_ONE_PROBABILITY = Fraction(1 << 18, (1 << 19) - 1)

# One common witness works for every fixed Q >= 3 and every BA weight.
SPACING_TILT = Fraction(127, 250)
CHERNOFF_TILT_PER_ROW = Fraction(133, 125)
COLLATZ_LIVE_COORDINATE = Fraction(3, 1600)
BLOCK_CONSTANT = Fraction(39, 4)


def upper_float(value) -> float:
    return math.nextafter(float(value.b), math.inf)


def fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": format(float(value), ".18g"),
    }


def entropy(value: Fraction):
    x = as_interval(value)
    return -x * iv.log(x) - (1 - x) * iv.log(1 - x)


def collatz_lines() -> tuple[tuple[Fraction, Fraction], ...]:
    """Return A_i,B_i for m(u)=max_i(A_i+B_i u)."""
    s = SPACING_TILT
    v = COLLATZ_LIVE_COORDINATE
    r = RESET_PROBABILITY
    return (
        (ONE, s * v),
        (s, r / v + s * (1 - r)),
    )


def propose_fugacity(weight: Fraction) -> Fraction:
    """Use binary64 only to propose a rational positive witness."""
    x = float(weight)
    lines = tuple((float(a), float(b)) for a, b in collatz_lines())
    candidates = [x * a / (b * (1 - x)) for a, b in lines]
    (a0, b0), (a1, b1) = lines
    if b0 != b1:
        crossing = (a1 - a0) / (b0 - b1)
        if crossing > 0:
            candidates.append(crossing)

    def objective(u: float) -> float:
        norm = max(a + b * u for a, b in lines)
        return -x * math.log(u) + math.log(norm)

    selected = min(candidates, key=objective)
    return Fraction(format(selected, ".17g"))


def spacing_exponent_candidates() -> list[tuple[str, object]]:
    """Return all candidates for G_3(s,tau)."""
    s = SPACING_TILT
    tau = CHERNOFF_TILT_PER_ROW
    d = (1 - s) / s
    critical = Fraction(4, 3) / tau - 1 / d
    if not ZERO < critical < ONE:
        raise AssertionError("the selected G_3 critical point is not interior")
    return [
        ("left_endpoint", as_interval(ZERO)),
        (
            "right_endpoint",
            -as_interval(tau)
            + as_interval(Fraction(4, 3)) * iv.log(1 / as_interval(s)),
        ),
        (
            "interior_critical_point",
            -as_interval(tau) * as_interval(critical)
            + as_interval(Fraction(4, 3))
            * iv.log(1 + as_interval(d) * as_interval(critical)),
        ),
    ]


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
        "--output",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "golay_ba3_rm2sub_weight_coupled_fixed_d11.json"
        ),
    )
    args = parser.parse_args()

    segments = load_segments(args.majorant)
    lines = collatz_lines()
    g_candidates = spacing_exponent_candidates()
    delta_term = (
        as_interval(DELTA)
        * as_interval(CHERNOFF_TILT_PER_ROW)
        / as_interval(LIVE_ONE_PROBABILITY)
    )
    support_term = iv.log(iv.mpf(2)) / as_interval(BLOCK_CONSTANT)

    checks = []
    largest_objective_upper = -math.inf
    largest_schedule_margin_upper = -math.inf
    for segment in segments:
        midpoint = (segment.lower + segment.upper) / 2
        fugacity = propose_fugacity(midpoint)
        norm = max(a + b * fugacity for a, b in lines)
        for endpoint in (segment.lower, segment.upper):
            x = as_interval(endpoint)
            outer = (
                as_interval(segment.slope) * x
                + as_interval(segment.intercept)
                - entropy(endpoint)
            )
            row_term = (
                outer
                - x * iv.log(as_interval(fugacity))
                + iv.log(as_interval(norm))
            )
            objective_upper = max(
                upper_float(row_term + g_value + delta_term)
                for _, g_value in g_candidates
            )
            schedule_margin_upper = max(
                upper_float(row_term + g_value + delta_term + support_term)
                for _, g_value in g_candidates
            )
            largest_objective_upper = max(
                largest_objective_upper, objective_upper
            )
            largest_schedule_margin_upper = max(
                largest_schedule_margin_upper, schedule_margin_upper
            )
            checks.append(
                {
                    "segment": segment.index,
                    "weight": str(endpoint),
                    "fugacity": str(fugacity),
                    "objective_upper_natural_per_active_row": objective_upper,
                    "schedule_margin_upper_natural_per_active_row": (
                        schedule_margin_upper
                    ),
                }
            )

    assertions = {
        "spacing_tilt_between_zero_and_one": (
            ZERO < SPACING_TILT < ONE
        ),
        "all_fugacities_positive": all(
            Fraction(record["fugacity"]) > 0 for record in checks
        ),
        "fixed_occupation_objective_negative": largest_objective_upper < 0,
        "block_schedule_pays_support_choice": (
            largest_schedule_margin_upper < 0
        ),
    }
    if not all(assertions.values()):
        raise AssertionError(f"weight-coupled assertion failed: {assertions}")

    payload = {
        "schema": "golay-ba3-rm2sub-weight-coupled-fixed-v1",
        "status": "proved",
        "claim": (
            "For every fixed Q>=3, the weight-coupled BA/RM2Sub transfer "
            "and B=(39/4) log2(N)+O(1) give a negative first-moment "
            "exponent after choosing the Q active rows."
        ),
        "scope": {
            "occupation_quantifier": "every fixed integer Q>=3",
            "distance": "11/100",
            "majorant_segments": len(segments),
            "convexity_rule": (
                "For fixed fugacity on one affine BA segment, the row "
                "objective is convex in the row weight. Its maximum is at "
                "a segment endpoint."
            ),
            "quantifier_order": (
                "For each BA weight segment choose one fugacity, then bound "
                "both endpoints outward; the spacing, Chernoff, and Collatz "
                "witnesses are common to all segments and every fixed Q>=3."
            ),
        },
        "witness": {
            "spacing_tilt": fraction_record(SPACING_TILT),
            "chernoff_tilt_per_active_row": fraction_record(
                CHERNOFF_TILT_PER_ROW
            ),
            "collatz_vector": ["1", str(COLLATZ_LIVE_COORDINATE)],
            "block_constant": fraction_record(BLOCK_CONSTANT),
        },
        "bounds": {
            "largest_objective_upper_natural_per_active_row": (
                largest_objective_upper
            ),
            "certified_decay_lower_natural_per_active_row": (
                -largest_objective_upper
            ),
            "largest_schedule_margin_upper_natural_per_active_row": (
                largest_schedule_margin_upper
            ),
            "certified_schedule_decay_lower_natural_per_active_row": (
                -largest_schedule_margin_upper
            ),
            "g3_candidate_uppers": {
                name: upper_float(value) for name, value in g_candidates
            },
            "support_cost_upper_natural_per_active_row": upper_float(
                support_term
            ),
        },
        "checks": checks,
        "assertions": assertions,
    }
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": payload["status"],
                "checks": len(checks),
                "bounds": payload["bounds"],
                "assertions": assertions,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Compute exact combinatorial bounds for the locator incidence reduction."""

from __future__ import annotations

import json
import math
from fractions import Fraction
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "generated"


def main() -> None:
    cap_summary = json.loads((GENERATED / "exact_bch_lp_caps.json").read_text())
    application_a38_log2 = cap_summary["shells"][0][
        "sufficient_application_cap_log2"
    ]

    h38_certificate = json.loads(
        (GENERATED / "max_h_38_scaled_rational.certificate.json").read_text()
    )
    h38_cap = int(h38_certificate["physical_lattice_cap"])
    assert h38_cap % 128 == 0
    lp_locator_cap = 19 * h38_cap // 128

    # If g'=0, the perfectness of GF(256) gives g=h^2 with deg(h)<=9 and
    # h(0)=1.  Any nine nonzero agreement points determine h uniquely.
    degenerate_pair_numerator = comb(255, 9)
    degenerate_pair_denominator = comb(37, 9)
    degenerate_locator_cap = (
        degenerate_pair_numerator // degenerate_pair_denominator
    )

    # For the 24-point incidence, six conditions remain after interpolating
    # the 18 free coefficients.  This is the random-rank reference value, not
    # a proved count for the smooth incidence variety.
    random_rank_locator_bound = Fraction(
        comb(255, 24), (256**6) * comb(37, 24)
    )
    application_locator_log2 = application_a38_log2 - math.log2(Fraction(3968, 19))
    smooth_factor_budget = 2**application_locator_log2 / float(
        random_rank_locator_bound
    )
    smooth_factor_budget_after_degenerate = (
        2**application_locator_log2 - degenerate_locator_cap
    ) / float(random_rank_locator_bound)

    result = {
        "classification": {
            "identities_and_degenerate_bound": "exact",
            "random_rank_reference": "heuristic comparator only",
            "application_threshold": "derived from binary64 RandomStepConv coefficients",
        },
        "counting_relations": {
            "L_37": "A_37(P_punctured)/255 = 19*h_38/128",
            "A_38_C": "3968*L_37/19",
        },
        "exact_lp_locator_cap": lp_locator_cap,
        "exact_lp_locator_cap_log2": math.log2(lp_locator_cap),
        "application_locator_threshold_log2": application_locator_log2,
        "remaining_gap_bits": math.log2(lp_locator_cap) - application_locator_log2,
        "jacobian_theorem": {
            "statement": (
                "Every admissible 24-set associated with g' != 0 has Jacobian rank 6."
            ),
            "reason": (
                "After invertible row and column scalings, any six columns "
                "with g'(t) != 0 form the Vandermonde matrix on "
                "1,t,t^2,t^3,t^4,t^5. A nonzero g' has degree at most 16, "
                "so at least eight of 24 columns survive."
            ),
        },
        "degenerate_family_g_prime_zero": {
            "reduction": "g=h^2, deg(h)<=9, h(0)=1, h(u)=u^73 on 37 points",
            "factorial_moment_numerator": degenerate_pair_numerator,
            "factorial_moment_denominator": degenerate_pair_denominator,
            "locator_count_upper": degenerate_locator_cap,
            "locator_count_upper_log2": math.log2(degenerate_locator_cap),
        },
        "smooth_24_point_target": {
            "random_rank_bound_exact": (
                f"{random_rank_locator_bound.numerator}/"
                f"{random_rank_locator_bound.denominator}"
            ),
            "random_rank_bound_log2": math.log2(float(random_rank_locator_bound)),
            "factor_budget_before_application_threshold": smooth_factor_budget,
            "smooth_factor_budget_after_degenerate_cap": (
                smooth_factor_budget_after_degenerate
            ),
        },
    }
    output = GENERATED / "locator_incidence_bounds.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

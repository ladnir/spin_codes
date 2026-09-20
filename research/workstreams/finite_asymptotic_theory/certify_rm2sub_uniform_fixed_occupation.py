#!/usr/bin/env python3
"""Exact rational certificate for the uniform fixed-occupation RM2Sub bound.

The proof is in RM2SUB_UNIFORM_FIXED_OCCUPATION.md.  This checker uses only
integer and rational arithmetic.  Decimal fields are displays of exact
rational bounds; they are not inputs to the certificate.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import isqrt
from pathlib import Path


def atanh_log_bounds(x: Fraction, terms: int) -> tuple[Fraction, Fraction]:
    """Bound log((1+x)/(1-x)) by a truncated positive series."""
    if not Fraction(0) < x < Fraction(1):
        raise ValueError("x must lie strictly between zero and one")
    if terms < 1:
        raise ValueError("terms must be positive")
    lower = 2 * sum(
        (x ** (2 * index + 1)) / (2 * index + 1)
        for index in range(terms)
    )
    first_tail_power = 2 * terms + 1
    tail_upper = (
        2
        * x**first_tail_power
        / first_tail_power
        / (1 - x * x)
    )
    return lower, lower + tail_upper


def rational_sqrt_upper(value: Fraction, decimal_digits: int) -> Fraction:
    """Return a rational decimal-grid upper bound on sqrt(value)."""
    if value < 0:
        raise ValueError("cannot take the square root of a negative number")
    scale = 10**decimal_digits
    scaled_numerator = value.numerator * scale * scale
    root = isqrt(scaled_numerator // value.denominator)
    if root * root * value.denominator < scaled_numerator:
        root += 1
    return Fraction(root, scale)


def fraction_record(value: Fraction, decimal_digits: int = 18) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal_display": f"{float(value):.{decimal_digits}g}",
    }


def build_certificate() -> dict[str, object]:
    state_bits = 19
    state_count = (1 << state_bits) - 1
    reset_probability = Fraction(1, state_count)
    live_one_probability = Fraction(1 << (state_bits - 1), state_count)
    delta = Fraction(11, 100)
    tau = Fraction(11, 10)
    y = Fraction(1, 2)
    minimum_occupation = 3
    block_constant = 9

    # log(2)=log((1+1/3)/(1-1/3)).
    log_two_lower, log_two_upper = atanh_log_bounds(Fraction(1, 3), 14)

    # For Q >= 3, G_Q decreases with Q.  At Q=3 the maximizer is
    # x=(4/3)/(11/10)-1=7/33, and 1+x=40/33.
    log_ratio_lower, log_ratio_upper = atanh_log_bounds(Fraction(7, 73), 10)
    del log_ratio_lower
    g_three_upper = -Fraction(7, 30) + Fraction(4, 3) * log_ratio_upper

    radicand = 2 * reset_probability + reset_probability**2 / 4
    sqrt_upper = rational_sqrt_upper(radicand, 40)
    perron_upper = (
        1 - reset_probability / 4 + sqrt_upper / 2
    )
    if perron_upper < 1:
        raise AssertionError("Perron upper bound must be at least one")

    # log(rho_*) <= rho_*-1.  All remaining operations are exact rational
    # upper or lower bounds in the indicated directions.
    natural_objective_upper = (
        (perron_upper - 1)
        + g_three_upper
        + delta * tau / live_one_probability
    )
    objective_bits_per_active_row_upper = natural_objective_upper / log_two_lower
    gap_per_active_row_lower = (
        Fraction(1, 2) - objective_bits_per_active_row_upper
    )
    strict_block_threshold_upper = 1 / gap_per_active_row_lower
    schedule_decay_exponent_per_active_row_lower = (
        block_constant * gap_per_active_row_lower - 1
    )

    assertions = {
        "log_two_bracket_is_strict": log_two_lower < log_two_upper,
        "perron_square_root_is_upper": (
            sqrt_upper * sqrt_upper >= radicand
        ),
        "gap_is_positive": gap_per_active_row_lower > 0,
        "block_constant_exceeds_threshold": (
            Fraction(block_constant) > strict_block_threshold_upper
        ),
        "schedule_decay_is_positive": (
            schedule_decay_exponent_per_active_row_lower > 0
        ),
    }
    if not all(assertions.values()):
        raise AssertionError(f"certificate assertion failed: {assertions}")

    return {
        "schema": "rm2sub-uniform-fixed-occupation-v1",
        "status": "EXACT_RATIONAL_CERTIFICATE",
        "scope": {
            "state_bits": state_bits,
            "delta": "11/100",
            "minimum_occupation": minimum_occupation,
            "occupation_quantifier": "every fixed integer Q >= 3",
            "selected_y": "1/2",
            "selected_tau": "11/10",
            "selected_block_constant": block_constant,
        },
        "exact_bounds": {
            "log_two_lower": fraction_record(log_two_lower),
            "log_two_upper": fraction_record(log_two_upper),
            "log_40_over_33_upper": fraction_record(log_ratio_upper),
            "g_three_upper": fraction_record(g_three_upper),
            "sqrt_radicand_upper": fraction_record(sqrt_upper),
            "perron_upper": fraction_record(perron_upper),
            "natural_objective_per_active_row_upper": fraction_record(
                natural_objective_upper
            ),
            "objective_bits_per_active_row_upper": fraction_record(
                objective_bits_per_active_row_upper
            ),
            "gap_bits_per_active_row_lower": fraction_record(
                gap_per_active_row_lower
            ),
            "strict_log2_block_threshold_upper": fraction_record(
                strict_block_threshold_upper
            ),
            "schedule_decay_exponent_per_active_row_lower": fraction_record(
                schedule_decay_exponent_per_active_row_lower
            ),
        },
        "assertions": assertions,
        "interpretation": (
            "For each fixed Q>=3, B=9 log2(N)+O(1) gives "
            "Pr[Z_{floor(0.11N),Q}>0] <= "
            "N^{-epsilon Q+o_Q(1)}, where epsilon exceeds the certified "
            "schedule-decay lower bound."
        ),
        "limitations": [
            "The continuum approximation is proved separately for each fixed Q.",
            "This certificate does not allow Q to grow with L or N.",
            "Occupations Q=1 and Q=2 use their sharper existing certificates.",
            "The structured-outer transfer still requires a proved uniform spectrum comparison.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

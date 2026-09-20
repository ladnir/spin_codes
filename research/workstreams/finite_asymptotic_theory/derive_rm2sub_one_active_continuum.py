#!/usr/bin/env python3
"""Evaluate the proved one-active RM2Sub continuum formula.

The formula is derived in RM2SUB_ONE_ACTIVE_CONTINUUM.md.  This program uses
binary64 arithmetic to report convenient constants; it is not an
outward-rounded certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from math import isqrt
from pathlib import Path

from scipy.optimize import minimize_scalar


def continuum_values(theta: float, *, state_bits: int, delta: float) -> dict[str, float]:
    state_count = (1 << state_bits) - 1
    live_one_probability = (1 << (state_bits - 1)) / state_count
    gamma = live_one_probability * theta
    live_region = math.exp(-gamma)
    averaged_suffix = -math.expm1(-gamma) / gamma if gamma else 1.0
    lower_diagonal = (2.0 - 1.0 / state_count) * live_region
    discriminant = (
        (1.0 - lower_diagonal) ** 2
        + 4.0 * averaged_suffix**2 / state_count
    )
    spectral_radius = (
        1.0 + lower_diagonal + math.sqrt(discriminant)
    ) / 2.0
    objective = (
        math.log2(spectral_radius)
        + delta * theta / math.log(2.0)
    )
    return {
        "theta": theta,
        "gamma": gamma,
        "live_region_factor": live_region,
        "averaged_suffix_factor": averaged_suffix,
        "spectral_radius": spectral_radius,
        "objective_bits_per_outer_coordinate": objective,
    }


def optimize_delta(delta: float, *, state_bits: int) -> dict[str, float]:
    result = minimize_scalar(
        lambda theta: continuum_values(
            theta, state_bits=state_bits, delta=delta
        )["objective_bits_per_outer_coordinate"],
        bounds=(0.0, 32.0),
        method="bounded",
        options={"xatol": 1e-14},
    )
    values = continuum_values(
        float(result.x), state_bits=state_bits, delta=delta
    )
    eta = 0.5 - values["objective_bits_per_outer_coordinate"]
    values.update(
        {
            "delta": delta,
            "one_active_gap_bits_per_outer_coordinate": eta,
            "strict_log2_block_threshold": 1.0 / eta,
            "optimizer_success": bool(result.success),
        }
    )
    return values


def rational_theta_log2_certificate(*, state_bits: int) -> dict[str, object]:
    """Certify c=18/5 at delta=5501/50000 using rational inequalities."""
    state_count = (1 << state_bits) - 1
    delta = Fraction(5501, 50000)
    live_one_probability = Fraction(1 << (state_bits - 1), state_count)

    # ln(2) = 2 atanh(1/3).  Every omitted term is positive, so this partial
    # sum is an exact rational lower bound.
    ln2_lower = 2 * sum(
        Fraction(1, 3) ** (2 * index + 1) / (2 * index + 1)
        for index in range(13)
    )
    radicand_upper = (
        Fraction(1, 4 * state_count * state_count)
        + Fraction(1, state_count) / (ln2_lower * ln2_lower)
    )
    scale = 10**30
    scaled_square = (
        radicand_upper.numerator * scale * scale
        + radicand_upper.denominator
        - 1
    ) // radicand_upper.denominator
    sqrt_numerator = isqrt(scaled_square)
    if sqrt_numerator * sqrt_numerator < scaled_square:
        sqrt_numerator += 1
    sqrt_upper = Fraction(sqrt_numerator, scale)
    lambda_upper = (
        1
        - Fraction(1, 4 * state_count)
        + sqrt_upper / 2
    )

    # At theta=ln(2)/p, gamma=ln(2), a=1/2, and
    # log_2(lambda) <= (lambda-1)/ln(2).
    log2_lambda_upper = (lambda_upper - 1) / ln2_lower
    eta_lower = (
        Fraction(1, 2)
        - delta / live_one_probability
        - log2_lambda_upper
    )
    c = Fraction(18, 5)
    polynomial_exponent_lower = c * eta_lower - 1
    if polynomial_exponent_lower <= 0:
        raise AssertionError("rational c=18/5 certificate did not close")

    def fraction_record(value: Fraction) -> dict[str, object]:
        return {
            "numerator": str(value.numerator),
            "denominator": str(value.denominator),
            "binary64_display": float(value),
        }

    return {
        "status": "EXACT_RATIONAL_INEQUALITY_CERTIFICATE",
        "delta": "5501/50000",
        "theta": "ln(2)/p",
        "ln2_lower_terms": 13,
        "ln2_lower": fraction_record(ln2_lower),
        "sqrt_upper_decimal_scale": scale,
        "lambda_upper": fraction_record(lambda_upper),
        "one_active_gap_lower": fraction_record(eta_lower),
        "block_constant": "18/5",
        "polynomial_decay_exponent_lower": fraction_record(
            polynomial_exponent_lower
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state-bits", type=int, default=19)
    parser.add_argument(
        "--deltas",
        type=float,
        nargs="+",
        default=[0.09, 0.10, 0.105, 0.1085, 0.11, 0.11002],
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    state_count = (1 << args.state_bits) - 1
    live_one_probability = (1 << (args.state_bits - 1)) / state_count
    rows = [
        optimize_delta(delta, state_bits=args.state_bits)
        for delta in args.deltas
    ]
    payload = {
        "schema": "rm2sub-one-active-continuum-v1",
        "status": "PROVED_FORMULA_BINARY64_OPTIMIZATION_NOT_OUTWARD_ROUNDED",
        "state_bits": args.state_bits,
        "nonzero_state_count": state_count,
        "live_coordinate_one_probability": live_one_probability,
        "limiting_matrices": {
            "inactive": [["1", "0"], ["0", "a"]],
            "marked": [
                ["0", "f"],
                ["f/M", "(1-1/M)a"],
            ],
            "definitions": {
                "M": "2^s-1",
                "gamma": "theta*2^(s-1)/M",
                "a": "exp(-gamma)",
                "f": "(1-exp(-gamma))/gamma",
            },
        },
        "optimization": rows,
        "rational_certificate_delta_0_11002": rational_theta_log2_certificate(
            state_bits=args.state_bits
        ),
        "selected_delta_0_11002_polynomial_exponents": {
            str(c): c * rows[-1]["one_active_gap_bits_per_outer_coordinate"] - 1.0
            for c in (3.6, 3.75, 4.0, 5.0, 8.0, 17.0)
        },
        "arithmetic": (
            "IEEE-754 binary64 with scipy.optimize.minimize_scalar; "
            "reported constants are reproducible diagnostics, not directed-rounding certificates"
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

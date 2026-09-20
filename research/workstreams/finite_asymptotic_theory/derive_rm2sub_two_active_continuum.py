#!/usr/bin/env python3
"""Evaluate the q=2 RM2Sub continuum formula and exact safe constant."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from math import isqrt
from pathlib import Path

from scipy.optimize import minimize_scalar


def continuum_values(theta: float, *, state_bits: int, delta: float) -> dict[str, object]:
    state_count = (1 << state_bits) - 1
    reset_probability = 1.0 / state_count
    survival_probability = 1.0 - reset_probability
    live_one_probability = (1 << (state_bits - 1)) / state_count
    gamma = live_one_probability * theta
    live_region = math.exp(-gamma)
    if gamma:
        one_impulse_average = -math.expm1(-gamma) / gamma
        gap_average = 2.0 * (gamma - 1.0 + live_region) / gamma**2
        edge_average = (
            2.0 * (1.0 - (1.0 + gamma) * live_region) / gamma**2
        )
    else:
        one_impulse_average = 1.0
        gap_average = 1.0
        edge_average = 1.0

    k0 = [[1.0, 0.0], [0.0, live_region]]
    k1 = [
        [0.0, one_impulse_average],
        [reset_probability * one_impulse_average, survival_probability * live_region],
    ]
    k2 = [
        [reset_probability * gap_average, survival_probability * edge_average],
        [
            survival_probability * reset_probability * edge_average,
            reset_probability * edge_average
            + survival_probability**2 * live_region,
        ],
    ]
    total = [
        [k0[row][column] + 2.0 * k1[row][column] + k2[row][column] for column in range(2)]
        for row in range(2)
    ]
    trace = total[0][0] + total[1][1]
    determinant = total[0][0] * total[1][1] - total[0][1] * total[1][0]
    spectral_radius = (
        trace + math.sqrt(trace**2 - 4.0 * determinant)
    ) / 2.0
    objective = math.log2(spectral_radius) + delta * theta / math.log(2.0)
    return {
        "theta": theta,
        "gamma": gamma,
        "live_region_factor": live_region,
        "one_impulse_average": one_impulse_average,
        "two_impulse_gap_average": gap_average,
        "two_impulse_edge_average": edge_average,
        "K0": k0,
        "K1": k1,
        "K2": k2,
        "K0_plus_2K1_plus_K2": total,
        "spectral_radius": spectral_radius,
        "objective_bits_per_outer_coordinate": objective,
    }


def optimize_delta(delta: float, *, state_bits: int) -> dict[str, object]:
    result = minimize_scalar(
        lambda theta: continuum_values(theta, state_bits=state_bits, delta=delta)[
            "objective_bits_per_outer_coordinate"
        ],
        bounds=(0.0, 64.0),
        method="bounded",
        options={"xatol": 1e-14},
    )
    values = continuum_values(float(result.x), state_bits=state_bits, delta=delta)
    eta = 1.0 - float(values["objective_bits_per_outer_coordinate"])
    values.update(
        {
            "delta": delta,
            "two_active_gap_bits_per_outer_coordinate": eta,
            "strict_log2_block_threshold": 2.0 / eta,
            "optimizer_success": bool(result.success),
        }
    )
    return values


def fraction_record(value: Fraction) -> dict[str, object]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "binary64_display": float(value),
    }


def rational_theta_log4_certificate(*, state_bits: int) -> dict[str, object]:
    """Certify c=18/5 at delta=5501/50000 with exact integers."""
    state_count = (1 << state_bits) - 1
    reset = Fraction(1, state_count)
    survive = 1 - reset
    live_one_probability = Fraction(1 << (state_bits - 1), state_count)
    delta = Fraction(5501, 50000)

    terms = 13
    ln2_lower = 2 * sum(
        Fraction(1, 3) ** (2 * index + 1) / (2 * index + 1)
        for index in range(terms)
    )
    last_index = terms - 1
    remainder_upper = (
        2
        * Fraction(1, 3) ** (2 * last_index + 3)
        / (2 * last_index + 3)
        / (1 - Fraction(1, 9))
    )
    ln2_upper = ln2_lower + remainder_upper
    gamma_lower = 2 * ln2_lower
    gamma_upper = 2 * ln2_upper

    def f(gamma: Fraction) -> Fraction:
        return Fraction(3, 4) / gamma

    def gap(gamma: Fraction) -> Fraction:
        return 2 * (gamma - Fraction(3, 4)) / (gamma * gamma)

    def edge(gamma: Fraction) -> Fraction:
        return 2 * (Fraction(3, 4) - gamma / 4) / (gamma * gamma)

    # On this narrow interval, f and edge decrease, while gap increases.
    f_lower, f_upper = f(gamma_upper), f(gamma_lower)
    gap_lower, gap_upper = gap(gamma_lower), gap(gamma_upper)
    edge_lower, edge_upper = edge(gamma_upper), edge(gamma_lower)

    def entries(
        fv: Fraction, gapv: Fraction, edgev: Fraction
    ) -> tuple[Fraction, Fraction, Fraction, Fraction]:
        live_region = Fraction(1, 4)
        return (
            1 + reset * gapv,
            2 * fv + survive * edgev,
            2 * reset * fv + survive * reset * edgev,
            live_region
            + 2 * survive * live_region
            + reset * edgev
            + survive * survive * live_region,
        )

    lower = entries(f_lower, gap_lower, edge_lower)
    upper = entries(f_upper, gap_upper, edge_upper)
    x_lower, _, _, w_lower = lower
    x_upper, y_upper, z_upper, w_upper = upper
    difference_upper = max(abs(x_upper - w_lower), abs(w_upper - x_lower))
    radicand_upper = difference_upper**2 + 4 * y_upper * z_upper

    scale = 10**40
    scaled_square = (
        radicand_upper.numerator * scale * scale
        + radicand_upper.denominator
        - 1
    ) // radicand_upper.denominator
    sqrt_numerator = isqrt(scaled_square)
    if sqrt_numerator * sqrt_numerator < scaled_square:
        sqrt_numerator += 1
    sqrt_upper = Fraction(sqrt_numerator, scale)
    lambda_upper = (x_upper + w_upper + sqrt_upper) / 2

    # At theta=2 ln(2)/p, the Chernoff contribution is exactly 2 delta/p.
    # Also, log_2(lambda) <= (lambda-1)/ln(2).
    log2_lambda_upper = (lambda_upper - 1) / ln2_lower
    eta_lower = 1 - 2 * delta / live_one_probability - log2_lambda_upper
    block_constant = Fraction(18, 5)
    decay_exponent_lower = block_constant * eta_lower - 2
    if decay_exponent_lower <= 0:
        raise AssertionError("rational c=18/5 q=2 certificate did not close")

    return {
        "status": "EXACT_RATIONAL_INEQUALITY_CERTIFICATE",
        "delta": "5501/50000",
        "theta": "2*ln(2)/p",
        "ln2_lower_terms": terms,
        "ln2_lower": fraction_record(ln2_lower),
        "ln2_upper": fraction_record(ln2_upper),
        "sqrt_upper_decimal_scale": scale,
        "lambda_upper": fraction_record(lambda_upper),
        "two_active_gap_lower": fraction_record(eta_lower),
        "block_constant": "18/5",
        "polynomial_decay_exponent_lower": fraction_record(decay_exponent_lower),
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

    rows = [optimize_delta(delta, state_bits=args.state_bits) for delta in args.deltas]
    payload = {
        "schema": "rm2sub-two-active-continuum-v1",
        "status": "PROVED_FORMULA_BINARY64_OPTIMIZATION_NOT_OUTWARD_ROUNDED",
        "state_bits": args.state_bits,
        "nonzero_state_count": (1 << args.state_bits) - 1,
        "same_epoch_collision_probability": "(t-1)/(L-1)",
        "optimization": rows,
        "rational_certificate_delta_0_11002": rational_theta_log4_certificate(
            state_bits=args.state_bits
        ),
        "selected_delta_0_11002_polynomial_exponents": {
            str(c): c * float(rows[-1]["two_active_gap_bits_per_outer_coordinate"]) - 2.0
            for c in (3.6, 3.75, 4.0, 5.0, 8.0, 17.0)
        },
        "arithmetic": (
            "IEEE-754 binary64 for optimized decimals; the c=18/5 receipt uses "
            "only exact rational and integer inequalities"
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Certify the RM2Sub four-state envelope for 0 < alpha <= 10^-4.

All polynomial coefficients and Bernstein sign tests use Fraction.  The only
floating-point values in the JSON receipt are decimal renderings of exact
rational quantities that are also included as numerator/denominator pairs.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from analyze_rm2sub_dense_occupation import (
    DEFAULT_SELECTION,
    association_matrix,
    state_code_weights,
    verify_transpose_columns,
)


Polynomial = list[Fraction]


def trim(poly: Polynomial) -> Polynomial:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def add(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [Fraction(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += value
    return trim(result)


def subtract(left: Polynomial, right: Polynomial) -> Polynomial:
    return add(left, [-value for value in right])


def scale(poly: Polynomial, factor: Fraction) -> Polynomial:
    return trim([factor * value for value in poly])


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        if left_value == 0:
            continue
        for right_index, right_value in enumerate(right):
            if right_value != 0:
                result[left_index + right_index] += left_value * right_value
    return trim(result)


def powers(base: Polynomial, maximum: int) -> list[Polynomial]:
    result = [[Fraction(1)]]
    for _ in range(maximum):
        result.append(multiply(result[-1], base))
    return result


def weighted_sum(
    polynomials: list[Polynomial], weights: list[int]
) -> Polynomial:
    result = [Fraction(0)]
    for poly, weight in zip(polynomials, weights, strict=True):
        result = add(result, scale(poly, Fraction(weight)))
    return result


def valuation(poly: Polynomial) -> int:
    for index, value in enumerate(poly):
        if value != 0:
            return index
    raise ValueError("the zero polynomial has no valuation")


def bernstein_coefficients_at_zero(
    poly: Polynomial, endpoint: Fraction
) -> list[Fraction]:
    """Convert a power polynomial on [0, endpoint] to Bernstein form."""
    degree = len(poly) - 1
    scaled_power = [
        coefficient * endpoint**index
        for index, coefficient in enumerate(poly)
    ]
    coefficients: list[Fraction] = []
    for bernstein_index in range(degree + 1):
        value = Fraction(0)
        for power_index in range(bernstein_index + 1):
            value += (
                scaled_power[power_index]
                * math.comb(bernstein_index, power_index)
                / math.comb(degree, power_index)
            )
        coefficients.append(value)
    return coefficients


def rational_payload(value: Fraction) -> dict[str, int | str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": format(float(value), ".17g"),
    }


def exponential_bounds(
    argument: Fraction, terms: int
) -> tuple[Fraction, Fraction]:
    """Return exact Taylor lower and geometric-tail upper bounds for exp."""
    partial = Fraction(1)
    term = Fraction(1)
    for index in range(1, terms + 1):
        term *= argument / index
        partial += term
    next_term = term * argument / (terms + 1)
    tail_ratio = argument / (terms + 2)
    upper = partial + next_term / (1 - tail_ratio)
    return partial, upper


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "rm2sub_dense_small_exact_d11.json"
        ),
    )
    args = parser.parse_args()

    payload = json.loads(args.selection.read_text(encoding="utf-8"))
    selected = payload["selected"]
    generator_words = [
        int(value, 16) for value in selected["A_generator_words_hex"]
    ]
    columns = [int(value, 16) for value in selected["B_columns_hex"]]
    verify_transpose_columns(
        generator_words=generator_words,
        columns=columns,
        output_bits=len(columns),
    )
    weights = state_code_weights(generator_words, len(columns))
    classes_array, counts_array, association_array = association_matrix(weights)
    classes = [int(value) for value in classes_array]
    counts = [int(value) for value in counts_array]
    association = [[int(value) for value in row] for row in association_array]

    state_space = 1 << len(generator_words)
    live_states = state_space - 1
    alpha = [Fraction(0), Fraction(1)]
    one = [Fraction(1)]
    schedule = Fraction(8, 5)
    beta = scale(alpha, schedule / 2)
    z = subtract(one, scale(alpha, schedule))
    input_complement = subtract(one, beta)
    u = add(input_complement, multiply(beta, z))
    v = add(beta, multiply(input_complement, z))
    signed = subtract(input_complement, multiply(beta, z))
    one_minus_p = subtract(one, scale(alpha, schedule))

    u_powers = powers(u, 128)
    v_powers = powers(v, 128)
    signed_powers = powers(signed, 128)
    z_powers = powers(z, 128)
    one_minus_p_powers = powers(one_minus_p, 128)

    live_by_class = [
        multiply(u_powers[128 - weight], v_powers[weight])
        for weight in classes
    ]
    syndrome_by_class = [
        multiply(u_powers[128 - weight], signed_powers[weight])
        for weight in classes
    ]
    zero_to_zero = scale(
        weighted_sum(syndrome_by_class, counts), Fraction(1, state_space)
    )
    zero_total = u_powers[128]
    zero_to_deterministic = subtract(zero_total, zero_to_zero)

    transformed_live: list[Polynomial] = []
    for row in association:
        transformed_live.append(weighted_sum(live_by_class, row))
    paired_total = [Fraction(0)]
    for syndrome_poly, transformed_poly in zip(
        syndrome_by_class, transformed_live, strict=True
    ):
        paired_total = add(
            paired_total, multiply(syndrome_poly, transformed_poly)
        )
    paired_total = scale(paired_total, Fraction(1, state_space))
    paired_nonzero = subtract(
        paired_total, multiply(zero_to_zero, live_by_class[0])
    )

    uniform_live = scale(
        subtract(weighted_sum(live_by_class, counts), live_by_class[0]),
        Fraction(1, live_states),
    )
    kernel_probability = scale(
        weighted_sum(
            [one_minus_p_powers[weight] for weight in classes], counts
        ),
        Fraction(1, state_space),
    )
    activation_probability = subtract(one, kernel_probability)
    correction = multiply(
        activation_probability,
        add(
            subtract(one, z_powers[128]),
            [Fraction(1024, live_states)],
        ),
    )
    radius_target = subtract(one, scale(alpha, Fraction(99)))
    test_live = Fraction(1, 1024)

    zero_inequality = subtract(
        add(
            zero_to_zero,
            scale(zero_to_deterministic, test_live),
        ),
        radius_target,
    )
    uniform_inequality = subtract(
        add(uniform_live, correction), radius_target
    )
    deterministic_inequality = subtract(
        add(
            paired_nonzero,
            multiply(zero_to_deterministic, correction),
        ),
        multiply(zero_to_deterministic, radius_target),
    )

    endpoint = Fraction(1, 10_000)
    named_polynomials = {
        "zero_row": zero_inequality,
        "uniform_and_punctured_rows": uniform_inequality,
        "deterministic_row_after_clearing_activation_mass": (
            deterministic_inequality
        ),
    }
    checks: dict[str, object] = {}
    all_negative = True
    for name, poly in named_polynomials.items():
        order = valuation(poly)
        quotient = poly[order:]
        bernstein = bernstein_coefficients_at_zero(quotient, endpoint)
        maximum = max(bernstein)
        minimum = min(bernstein)
        strictly_negative = maximum < 0
        all_negative = all_negative and strictly_negative
        checks[name] = {
            "original_degree": len(poly) - 1,
            "zero_order": order,
            "quotient_degree": len(quotient) - 1,
            "bernstein_coefficient_count": len(bernstein),
            "minimum_bernstein_coefficient": rational_payload(minimum),
            "maximum_bernstein_coefficient": rational_payload(maximum),
            "all_bernstein_coefficients_strictly_negative": (
                strictly_negative
            ),
        }

    # Natural-log upper bound for the complete exponent.  For
    # 0 < alpha <= endpoint, use log(1-x) <= -x and
    # -log(1-schedule*alpha) <= schedule*alpha/(1-schedule*endpoint).
    # The displayed decimal is diagnostic; the inequality itself is reduced
    # to an exact rational comparison after using log(2) < 7/10 and
    # log(8/5) > 47/100.
    log_rational_upper = (
        Fraction(7, 20)
        - Fraction(47, 100)
        - 1
        + Fraction(111, 100)
        * schedule
        / (1 - schedule * endpoint)
        - Fraction(99, 128)
    )
    exponent_negative = log_rational_upper < 0
    exp_seven_tenths_lower, _ = exponential_bounds(Fraction(7, 10), 12)
    _, exp_forty_seven_hundredths_upper = exponential_bounds(
        Fraction(47, 100), 12
    )
    logarithm_bounds_verified = (
        exp_seven_tenths_lower > 2
        and exp_forty_seven_hundredths_upper < Fraction(8, 5)
    )

    result = {
        "schema": "rm2sub-dense-small-exact-v1",
        "status": "EXACT_RATIONAL_BERNSTEIN_CERTIFICATE",
        "selection": str(args.selection),
        "delta": {"numerator": 11, "denominator": 100},
        "alpha_interval": {
            "lower": "0 (open)",
            "upper": rational_payload(endpoint),
        },
        "schedule": {
            "candidate_probability_over_alpha": rational_payload(schedule),
            "z": "1-(8/5)alpha",
        },
        "collatz_test_vector": [
            "1",
            "1/1024",
            "1/1024",
            f"{live_states}/{live_states - 1}/1024",
        ],
        "perron_radius_upper_bound": "1-99alpha",
        "polynomial_checks": checks,
        "log_bound": {
            "method": [
                "log(2) < 7/10",
                "log(8/5) > 47/100",
                "log(1-x) <= -x",
                "-log(1-x) <= x/(1-x_max)",
            ],
            "coefficient_upper_bound_natural": rational_payload(
                log_rational_upper
            ),
            "strictly_negative": exponent_negative,
            "logarithm_bounds_verified_by_exact_exponential_series": (
                logarithm_bounds_verified
            ),
        },
        "all_checks_pass": (
            all_negative and exponent_negative and logarithm_bounds_verified
        ),
        "limitations": [
            "This receipt certifies only 0 < alpha <= 10^-4.",
            "A separate certificate is required for alpha >= 10^-4.",
        ],
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if not result["all_checks_pass"]:
        raise SystemExit("the exact small-density certificate failed")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Certify the dominant-character likelihood variance bound.

The exact algebraic reduction is described in the companion probe.  Arb
encloses the one-dimensional reference laws.  Directed binary64 operations
then upper-bound every positive pairwise excess and their complete sum.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import json
import math
import time
from pathlib import Path

from flint import arb, arb_poly, ctx, fmpq
import numpy as np

from analyze_dual_walk_and_accumulator_energy import (
    accumulator_joint_counts,
    krawtchouk_row,
)
from certify_primal_schur_diagonal_outward import (
    add_lower_scalar,
    add_upper_scalar,
    coefficient_lower,
    coefficient_upper,
    down,
    positive_divide_upper,
    positive_multiply_lower,
    positive_multiply_upper,
    positive_power_upper,
    positive_sum_upper,
    rational_interval,
    up,
)
from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "dominant_character_likelihood_outward.json"


@lru_cache(maxsize=None)
def cached_joint_counts(output_bits: int) -> list[dict[int, int]]:
    return accumulator_joint_counts(output_bits)


def arb_float_interval(value: arb) -> tuple[np.float64, np.float64]:
    if not value.is_finite():
        raise ArithmeticError(f"non-finite Arb value: {value}")
    low = math.nextafter(float(value.lower()), -math.inf)
    high = math.nextafter(float(value.upper()), math.inf)
    if high < 0:
        raise ArithmeticError(f"negative probability ratio: {value}")
    return np.float64(max(0.0, low)), np.float64(max(0.0, high))


def reference_ratio_intervals(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    joint = cached_joint_counts(output_bits)
    krawtchouk_values = krawtchouk_row(output_bits, shell_weight)
    shell_size = math.comb(output_bits, shell_weight)
    one_polynomial = arb_poly(
        [
            arb(
                fmpq(
                    sum(
                        count
                        * krawtchouk_values[derivative_weight]
                        for derivative_weight, count in joint[level].items()
                    ),
                    shell_size,
                )
            )
            for level in range(output_bits + 1)
        ]
    )
    pair_polynomial = arb_poly(
        [
            arb(
                fmpq(
                    sum(
                        count
                        * krawtchouk_values[derivative_weight] ** 2
                        for derivative_weight, count in joint[level].items()
                    ),
                    shell_size * shell_size,
                )
            )
            for level in range(output_bits + 1)
        ]
    )
    denominator = math.comb(message_bits, right_degree)
    one_low = np.empty(message_bits + 1, dtype=np.float64)
    one_high = np.empty(message_bits + 1, dtype=np.float64)
    pair_low = np.empty(message_bits + 1, dtype=np.float64)
    pair_high = np.empty(message_bits + 1, dtype=np.float64)
    random_probability = arb(fmpq(shell_size, 1 << output_bits))
    for weight in range(message_bits + 1):
        numerator = krawtchouk(message_bits, right_degree, weight)
        if numerator == denominator:
            one_value = arb(1) / random_probability if shell_weight == 0 else arb(0)
            pair_value = arb(1) / random_probability
        elif numerator == -denominator:
            alternating_weight = (output_bits + 1) // 2
            one_value = (
                arb(1) / random_probability
                if shell_weight == alternating_weight
                else arb(0)
            )
            if 2 * shell_weight < alternating_weight:
                pair_value = arb(0)
            else:
                pair_value = pair_polynomial(arb(-1))
        else:
            bias = arb(fmpq(numerator, denominator))
            one_value = one_polynomial(bias)
            pair_value = pair_polynomial(bias)
        one_low[weight], one_high[weight] = arb_float_interval(one_value)
        pair_low[weight], pair_high[weight] = arb_float_interval(pair_value)
    return one_low, one_high, pair_low, pair_high


def likelihood_excess_upper(delta_upper: np.ndarray, exponent: int) -> np.ndarray:
    """Bound (1+delta)^exponent-1 without rounding a tiny delta into one."""
    scaled = positive_multiply_upper(
        np.full_like(delta_upper, np.float64(exponent)), delta_upper
    )
    result = np.empty_like(delta_upper)
    small = scaled < 0.5
    if np.any(small):
        denominator = down(np.float64(1) - scaled[small])
        result[small] = positive_divide_upper(scaled[small], denominator)
    if np.any(~small):
        base = up(np.float64(1) + delta_upper[~small])
        with np.errstate(over="ignore", invalid="ignore"):
            powered = positive_power_upper(base, exponent)
        result[~small] = up(powered - np.float64(1))
    return result


def certify(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weights: list[int],
    progress_every: int,
) -> tuple[list[dict[str, object]], int, float]:
    denominator = math.comb(message_bits, right_degree)
    absolute_bias_intervals = [
        rational_interval(
            abs(krawtchouk(message_bits, right_degree, weight)), denominator
        )
        for weight in range(message_bits + 1)
    ]
    absolute_low = np.asarray(
        [row[0] for row in absolute_bias_intervals], dtype=np.float64
    )
    absolute_high = np.asarray(
        [row[1] for row in absolute_bias_intervals], dtype=np.float64
    )
    reference = {
        shell: reference_ratio_intervals(
            message_bits, output_bits, right_degree, shell
        )
        for shell in shell_weights
    }
    random_probability_intervals = {
        shell: rational_interval(
            math.comb(output_bits, shell), 1 << output_bits
        )
        for shell in shell_weights
    }
    positive_excess = {shell: np.float64(0) for shell in shell_weights}
    pair_types = 0
    started = time.perf_counter()

    for n11 in range(message_bits + 1):
        first_factor = math.comb(message_bits, n11)
        remaining = message_bits - n11
        for n10 in range(remaining + 1):
            second_factor = math.comb(remaining, n10)
            last = remaining - n10
            n01 = np.arange(last + 1, dtype=np.int64)
            first_weight = n10 + n11
            second_weight = n01 + n11
            difference_weight = n10 + n01
            valid = (
                (first_weight != 0)
                & (second_weight != 0)
                & (difference_weight != 0)
            )
            if not np.any(valid):
                continue
            n01 = n01[valid]
            second_weight = second_weight[valid]
            difference_weight = difference_weight[valid]
            first_weights = np.full_like(second_weight, first_weight)
            multiplicity_upper = np.asarray(
                [
                    coefficient_upper(
                        first_factor
                        * second_factor
                        * math.comb(last, int(value))
                    )
                    for value in n01
                ],
                dtype=np.float64,
            )
            weights_upper = np.ldexp(multiplicity_upper, -2 * message_bits)

            nearest = np.stack(
                (
                    absolute_high[first_weights],
                    absolute_high[second_weight],
                    absolute_high[difference_weight],
                )
            )
            dominant = np.argmax(nearest, axis=0)
            first_residual = np.where(
                dominant == 0,
                absolute_high[second_weight],
                absolute_high[first_weights],
            )
            second_residual = np.where(
                dominant == 2,
                absolute_high[second_weight],
                absolute_high[difference_weight],
            )
            # Correct the middle case: when beta(y) is retained, the residuals
            # are beta(x) and beta(x+y).
            second_residual = np.where(
                dominant == 1,
                absolute_high[difference_weight],
                second_residual,
            )
            residual_upper = up(first_residual + second_residual)
            dominant_upper = np.where(
                dominant == 0,
                absolute_high[first_weights],
                np.where(
                    dominant == 1,
                    absolute_high[second_weight],
                    absolute_high[difference_weight],
                ),
            )
            denominator_lower = down(np.float64(1) - dominant_upper)
            delta_upper = np.full_like(denominator_lower, math.inf)
            regular = denominator_lower > 0
            delta_upper[regular] = positive_divide_upper(
                residual_upper[regular], denominator_lower[regular]
            )
            likelihood_excess = likelihood_excess_upper(
                delta_upper, output_bits
            )

            for shell in shell_weights:
                one_low, one_high, _pair_low, pair_high = reference[shell]
                first_low = one_low[first_weights]
                second_low = one_low[second_weight]
                first_high = one_high[first_weights]
                second_high = one_high[second_weight]
                base_high = np.where(
                    dominant == 0,
                    first_high,
                    np.where(
                        dominant == 1,
                        second_high,
                        pair_high[difference_weight],
                    ),
                )
                with np.errstate(over="ignore", invalid="ignore"):
                    increment = positive_multiply_upper(
                        base_high, likelihood_excess
                    )
                    candidate_upper = up(base_high + increment)
                candidate_upper[~np.isfinite(likelihood_excess)] = math.inf
                probability_lower = np.float64(
                    random_probability_intervals[shell][0]
                )
                marginal_first = positive_divide_upper(
                    first_high, np.full_like(first_high, probability_lower)
                )
                marginal_second = positive_divide_upper(
                    second_high, np.full_like(second_high, probability_lower)
                )
                pair_upper = np.minimum(
                    candidate_upper, np.minimum(marginal_first, marginal_second)
                )
                impossible = (first_high == 0) | (second_high == 0)
                pair_upper[impossible] = 0
                independent_lower = np.maximum(
                    down(first_low * second_low), 0
                )
                excess_upper = np.maximum(
                    up(pair_upper - independent_lower), 0
                )
                terms = positive_multiply_upper(weights_upper, excess_upper)
                positive_excess[shell] = add_upper_scalar(
                    positive_excess[shell], positive_sum_upper(terms)
                )
            pair_types += len(n01)
        if progress_every and (n11 % progress_every == 0 or n11 == message_bits):
            print(
                f"progress,n11,{n11},{message_bits},pair_types,{pair_types},elapsed_seconds,{time.perf_counter()-started:.3f}",
                flush=True,
            )

    rows = []
    for shell in shell_weights:
        one_low, one_high, _pair_low, _pair_high = reference[shell]
        mean_lower = np.float64(0)
        mean_upper = np.float64(0)
        for weight in range(1, message_bits + 1):
            low_weight = np.ldexp(
                coefficient_lower(math.comb(message_bits, weight)), -message_bits
            )
            high_weight = np.ldexp(
                coefficient_upper(math.comb(message_bits, weight)), -message_bits
            )
            mean_lower = add_lower_scalar(
                mean_lower,
                positive_multiply_lower(low_weight, one_low[weight]),
            )
            mean_upper = add_upper_scalar(
                mean_upper,
                positive_multiply_upper(
                    np.asarray([high_weight]), np.asarray([one_high[weight]])
                )[0],
            )
        random_mean_upper = np.ldexp(
            coefficient_upper(math.comb(output_bits, shell)),
            message_bits - output_bits,
        )
        extra_upper = positive_divide_upper(
            positive_multiply_upper(
                np.asarray([random_mean_upper]),
                np.asarray([positive_excess[shell]]),
            ),
            np.asarray([mean_lower]),
        )[0]
        variance_ratio_upper = add_upper_scalar(np.float64(1), extra_upper)
        rows.append(
            {
                "shell_weight": shell,
                "mean_ratio_lower": float(mean_lower),
                "mean_ratio_upper": float(mean_upper),
                "positive_pair_excess_normalized_upper": float(
                    positive_excess[shell]
                ),
                "variance_to_mean_upper": float(variance_ratio_upper),
                "passes_factor_512": bool(variance_ratio_upper <= 512),
            }
        )
    return rows, pair_types, time.perf_counter() - started


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument("--right-degree", type=int, default=33)
    parser.add_argument("--shell-weight", type=int, action="append", default=[])
    parser.add_argument("--shell-min", type=int)
    parser.add_argument("--shell-max", type=int)
    parser.add_argument("--precision", type=int, default=1536)
    parser.add_argument("--progress-every", type=int, default=16)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    shells = list(args.shell_weight)
    if args.shell_min is not None or args.shell_max is not None:
        if args.shell_min is None or args.shell_max is None:
            parser.error("shell-min and shell-max must be supplied together")
        shells.extend(range(args.shell_min, args.shell_max + 1))
    shells = sorted(set(shells or [42]))
    if any(not 0 <= shell <= args.output_bits for shell in shells):
        parser.error("invalid shell weight")
    ctx.prec = args.precision
    ctx.threads = 1
    rows, pair_types, elapsed = certify(
        args.message_bits,
        args.output_bits,
        args.right_degree,
        shells,
        args.progress_every,
    )
    payload = {
        "schema": "pure-ea-dominant-character-likelihood-outward-v1",
        "status": "OUTWARD_ARB_AND_BINARY64_CERTIFICATE",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "arb_precision_bits": args.precision,
        },
        "ordered_pair_types": pair_types,
        "elapsed_seconds": elapsed,
        "shells": rows,
        "claim": {
            "all_requested_shells_pass_factor_512": all(
                row["passes_factor_512"] for row in rows
            ),
            "maximum_variance_to_mean_upper": max(
                row["variance_to_mean_upper"] for row in rows
            ),
        },
        "proved_reduction": [
            "For each ordered nonzero distinct pair, the reference row law retains one bias of maximum absolute value.",
            "The omitted biases give a pointwise likelihood ratio at most (1+delta)^n.",
            "The pair probability is also at most either marginal probability.",
            "Replacing each covariance contribution by the positive part of its certified upper endpoint gives a valid variance upper bound.",
        ],
        "scope": [
            "Arb encloses every one-word and correlated-reference polynomial evaluation.",
            "Binary64 endpoints are rounded outward after each basic operation; nonnegative reductions include a summation-error allowance.",
            "The receipt covers exactly the requested shells and does not infer any uncomputed shell by symmetry.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

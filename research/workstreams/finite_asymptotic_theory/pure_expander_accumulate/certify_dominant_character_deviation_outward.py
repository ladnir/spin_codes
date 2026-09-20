#!/usr/bin/env python3
"""Certify a cancellation-free dominant-character variance bound.

The likelihood comparison in ``certify_dominant_character_likelihood_outward``
is sharp enough on the defect shells, but subtracting two ratios near one in
binary64 loses the central-shell covariance scale.  This verifier performs
that subtraction algebraically.  Arb first encloses the deviations of each
reference ratio from one.  Every subsequent operation combines nonnegative
upper endpoints with directed binary64 arithmetic.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
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
from certify_dominant_character_likelihood_outward import (
    arb_float_interval,
    likelihood_excess_upper,
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
    positive_sum_upper,
    rational_interval,
    up,
)
from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "dominant_character_deviation_outward.json"


@lru_cache(maxsize=None)
def cached_joint_counts(output_bits: int) -> list[dict[int, int]]:
    return accumulator_joint_counts(output_bits)


def nonnegative_arb_upper(value: arb) -> np.float64:
    if not value.is_finite():
        raise ArithmeticError(f"non-finite Arb value: {value}")
    high = math.nextafter(float(value.upper()), math.inf)
    if high < 0:
        raise ArithmeticError(f"negative upper endpoint: {value}")
    return np.float64(max(0.0, high))


def reference_intervals(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return one-word ratio endpoints and deviations from one.

    The fourth array contains deviations for the unbiased correlated-pair
    reference.  Its correlation parameter is the row-character bias.
    """

    joint = cached_joint_counts(output_bits)
    krawtchouk_values = krawtchouk_row(output_bits, shell_weight)
    shell_size = math.comb(output_bits, shell_weight)
    one_polynomial = arb_poly(
        [
            arb(
                fmpq(
                    sum(
                        count * krawtchouk_values[derivative_weight]
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
                        count * krawtchouk_values[derivative_weight] ** 2
                        for derivative_weight, count in joint[level].items()
                    ),
                    shell_size * shell_size,
                )
            )
            for level in range(output_bits + 1)
        ]
    )
    denominator = math.comb(message_bits, right_degree)
    random_probability = arb(fmpq(shell_size, 1 << output_bits))
    one_low = np.empty(message_bits + 1, dtype=np.float64)
    one_high = np.empty(message_bits + 1, dtype=np.float64)
    one_deviation = np.empty(message_bits + 1, dtype=np.float64)
    pair_deviation = np.empty(message_bits + 1, dtype=np.float64)
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
        one_deviation[weight] = nonnegative_arb_upper(abs(one_value - 1))
        pair_deviation[weight] = nonnegative_arb_upper(abs(pair_value - 1))
    return one_low, one_high, one_deviation, pair_deviation


def add_arrays_upper(*terms: np.ndarray) -> np.ndarray:
    result = np.zeros_like(terms[0])
    for term in terms:
        result = up(result + term)
    return result


def product_arrays_upper(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    with np.errstate(over="ignore", invalid="ignore"):
        return positive_multiply_upper(left, right)


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
    absolute_high = np.asarray(
        [row[1] for row in absolute_bias_intervals], dtype=np.float64
    )
    reference = {
        shell: reference_intervals(
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

            # If x+y is the all-ones message, the reference support remains
            # exact.  The actual-to-reference atom ratio is at most
            # 1+|beta(x)|.  Pairs with x or y all-ones have zero covariance,
            # because their shell indicator is deterministic.
            complement_pair = (
                (dominant == 2) & (difference_weight == message_bits)
            )
            delta_upper[complement_pair] = absolute_high[
                first_weights[complement_pair]
            ]
            likelihood_excess = likelihood_excess_upper(
                delta_upper, output_bits
            )
            deterministic_marginal = (
                (first_weights == message_bits) | (second_weight == message_bits)
            )

            for shell in shell_weights:
                one_low, one_high, one_deviation, pair_deviation = reference[shell]
                dx = one_deviation[first_weights]
                dy = one_deviation[second_weight]
                dz = pair_deviation[difference_weight]
                ell = likelihood_excess

                # Retaining beta(x) gives
                #   u_x(1+ell)-u_x u_y
                # <= d_y+d_x d_y+ell+d_x ell.
                excess_x = add_arrays_upper(
                    dy,
                    product_arrays_upper(dx, dy),
                    ell,
                    product_arrays_upper(dx, ell),
                )
                excess_y = add_arrays_upper(
                    dx,
                    product_arrays_upper(dy, dx),
                    ell,
                    product_arrays_upper(dy, ell),
                )
                # Retaining beta(x+y) gives
                #   b_z(1+ell)-u_x u_y
                # <= d_z+d_x+d_y+d_x d_y+ell+d_z ell.
                excess_z = add_arrays_upper(
                    dz,
                    dx,
                    dy,
                    product_arrays_upper(dx, dy),
                    ell,
                    product_arrays_upper(dz, ell),
                )
                candidate = np.where(
                    dominant == 0,
                    excess_x,
                    np.where(dominant == 1, excess_y, excess_z),
                )
                candidate[deterministic_marginal] = 0

                probability_lower = np.float64(
                    random_probability_intervals[shell][0]
                )
                marginal_first = positive_divide_upper(
                    one_high[first_weights],
                    np.full_like(first_weights, probability_lower, dtype=np.float64),
                )
                marginal_second = positive_divide_upper(
                    one_high[second_weight],
                    np.full_like(second_weight, probability_lower, dtype=np.float64),
                )
                excess_upper = np.minimum(
                    candidate, np.minimum(marginal_first, marginal_second)
                )
                impossible = (
                    (one_high[first_weights] == 0)
                    | (one_high[second_weight] == 0)
                )
                excess_upper[impossible] = 0
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
        one_low, one_high, _one_deviation, _pair_deviation = reference[shell]
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


def certify_worker(
    arguments: tuple[int, int, int, list[int], int, int]
) -> tuple[list[dict[str, object]], int, float]:
    message_bits, output_bits, right_degree, shells, progress_every, precision = arguments
    ctx.prec = precision
    ctx.threads = 1
    return certify(
        message_bits,
        output_bits,
        right_degree,
        shells,
        progress_every,
    )


def partition_shells(shells: list[int], workers: int) -> list[list[int]]:
    count = min(workers, len(shells))
    return [shells[index::count] for index in range(count)]


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
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="independent shell workers; each worker uses one FLINT thread",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    shells = list(args.shell_weight)
    if args.shell_min is not None or args.shell_max is not None:
        if args.shell_min is None or args.shell_max is None:
            parser.error("shell-min and shell-max must be supplied together")
        shells.extend(range(args.shell_min, args.shell_max + 1))
    shells = sorted(set(shells or [256]))
    if any(not 0 <= shell <= args.output_bits for shell in shells):
        parser.error("invalid shell weight")
    if args.workers < 1:
        parser.error("workers must be positive")
    started = time.perf_counter()
    shell_groups = partition_shells(shells, args.workers)
    if len(shell_groups) == 1:
        rows, pair_types, _worker_elapsed = certify_worker(
            (
                args.message_bits,
                args.output_bits,
                args.right_degree,
                shell_groups[0],
                args.progress_every,
                args.precision,
            )
        )
        worker_elapsed = [_worker_elapsed]
    else:
        worker_arguments = [
            (
                args.message_bits,
                args.output_bits,
                args.right_degree,
                group,
                args.progress_every,
                args.precision,
            )
            for group in shell_groups
        ]
        with ProcessPoolExecutor(max_workers=len(shell_groups)) as executor:
            worker_results = list(executor.map(certify_worker, worker_arguments))
        rows = sorted(
            [row for result, _types, _elapsed in worker_results for row in result],
            key=lambda row: int(row["shell_weight"]),
        )
        pair_type_counts = {types for _result, types, _elapsed in worker_results}
        if len(pair_type_counts) != 1:
            raise AssertionError("workers covered different pair-type domains")
        pair_types = pair_type_counts.pop()
        worker_elapsed = [elapsed for _result, _types, elapsed in worker_results]
    elapsed = time.perf_counter() - started
    payload = {
        "schema": "pure-ea-dominant-character-deviation-outward-v1",
        "status": "OUTWARD_ARB_AND_BINARY64_CERTIFICATE",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "arb_precision_bits": args.precision,
            "shell_workers": len(shell_groups),
            "flint_threads_per_worker": 1,
        },
        "ordered_pair_types": pair_types,
        "elapsed_seconds": elapsed,
        "worker_elapsed_seconds": worker_elapsed,
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
            "For each ordered nonzero distinct pair, retain one row-character bias of maximum absolute value.",
            "Arb encloses each one-word and correlated-reference ratio and its absolute deviation from one.",
            "The pair excess over independent marginals is bounded by a sum of nonnegative deviation and likelihood terms; no near-one subtraction occurs in binary64.",
            "Pairs with a deterministic all-ones marginal have zero covariance. Complement pairs use their exact two-atom support and a likelihood factor 1+abs(beta(x)).",
        ],
        "scope": [
            "Arb encloses every reference polynomial evaluation and deviation.",
            "Binary64 combines only nonnegative upper endpoints and advances each basic operation outward.",
            "The receipt covers exactly the requested shells and assumes no realized-spectrum complement symmetry.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

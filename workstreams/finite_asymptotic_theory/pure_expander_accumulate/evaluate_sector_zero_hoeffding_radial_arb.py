#!/usr/bin/env python3
"""Evaluate a sector-zero Hoeffding order range with Arb intervals.

The fixed-order radial formula is exact.  A range starting at order one gives
a certified lower bound on the sector-zero diagonal because all omitted
fixed-order quadratic forms are nonnegative.  A strict subrange is only an
enclosure of that range sum.  Neither is an upper bound until every omitted
order has been included or separately bounded.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

from flint import arb, arb_mat, arb_poly, ctx, fmpq

from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "sector_zero_hoeffding_radial_arb.json"


def q(value: fmpq) -> arb:
    return arb(value)


def require_finite(value: arb, context: str) -> arb:
    if not value.is_finite():
        raise ArithmeticError(f"non-finite Arb value at {context}: {value}")
    return value


def coefficient(
    output_bits: int, order: int, level: int, probability: arb
) -> arb:
    result = arb(0)
    low = max(0, level - (output_bits - order))
    high = min(order, level)
    for selected in range(low, high + 1):
        remainder = level - selected
        result += (
            (-1) ** selected
            * math.comb(order, selected)
            * math.comb(output_bits - order, remainder)
            * probability**remainder
            * (1 - probability) ** (output_bits - order - remainder)
        )
    return result


def radial_walk(
    message_bits: int, right_degree: int, maximum_order: int
) -> list[list[arb]]:
    denominator = math.comb(message_bits, right_degree)
    shell = [arb(0) for _ in range(message_bits + 1)]
    shell[0] = arb(1)
    result = [
        [
            shell[weight] / math.comb(message_bits, weight)
            for weight in range(message_bits + 1)
        ]
    ]
    for _step in range(maximum_order):
        following = [arb(0) for _ in range(message_bits + 1)]
        for weight, mass in enumerate(shell):
            low = max(0, right_degree - (message_bits - weight))
            high = min(weight, right_degree)
            for overlap in range(low, high + 1):
                next_weight = weight + right_degree - 2 * overlap
                transition = q(
                    fmpq(
                        math.comb(weight, overlap)
                        * math.comb(
                            message_bits - weight, right_degree - overlap
                        ),
                        denominator,
                    )
                )
                following[next_weight] += mass * transition
        shell = following
        result.append(
            [
                shell[weight] / math.comb(message_bits, weight)
                for weight in range(message_bits + 1)
            ]
        )
    return result


def krawtchouk_matrix(message_bits: int) -> arb_mat:
    return arb_mat(
        [
            [
                krawtchouk(message_bits, source_weight, frequency_weight)
                for source_weight in range(message_bits + 1)
            ]
            for frequency_weight in range(message_bits + 1)
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument("--right-degree", type=int, default=33)
    parser.add_argument("--level", type=int, required=True)
    parser.add_argument("--minimum-order", type=int, default=1)
    parser.add_argument("--maximum-order", type=int, default=8)
    parser.add_argument("--precision", type=int, default=384)
    parser.add_argument(
        "--transform-mode",
        choices=("scalar", "batched", "polynomial"),
        default="scalar",
        help=(
            "evaluate radial transforms separately, in one dense matrix "
            "product, or by the exact row-XOR polynomial identity"
        ),
    )
    parser.add_argument("--include-convolution-terms", action="store_true")
    parser.add_argument("--suppress-json", action="store_true")
    parser.add_argument("--progress-every", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not 0 <= args.level <= args.output_bits:
        parser.error("invalid level")
    if not 1 <= args.minimum_order <= args.maximum_order <= args.output_bits:
        parser.error("invalid order range")
    if args.precision < 128:
        parser.error("precision must be at least 128 bits")
    ctx.prec = args.precision
    ctx.threads = 1
    started = time.perf_counter()

    denominator = math.comb(args.message_bits, args.right_degree)
    bias_rational = [
        fmpq(
            krawtchouk(args.message_bits, args.right_degree, weight),
            denominator,
        )
        for weight in range(args.message_bits + 1)
    ]
    biases = [q(value) for value in bias_rational]
    probabilities = [(1 - value) / 2 for value in biases]
    walk = radial_walk(
        args.message_bits,
        args.right_degree,
        (
            args.output_bits
            if args.transform_mode == "polynomial"
            else args.maximum_order
        ),
    )
    transform = (
        None
        if args.transform_mode == "polynomial"
        else krawtchouk_matrix(args.message_bits)
    )
    walk_polynomials = (
        [
            arb_poly(
                [
                    walk[power][frequency_weight]
                    for power in range(args.output_bits + 1)
                ]
            )
            for frequency_weight in range(args.message_bits + 1)
        ]
        if args.transform_mode == "polynomial"
        else None
    )
    multiplicities = [
        math.comb(args.message_bits, weight)
        for weight in range(args.message_bits + 1)
    ]
    order_rows = []
    partial = arb(0)
    for order in range(args.minimum_order, args.maximum_order + 1):
        shell_coefficients = (
            [
                coefficient(args.output_bits, order, args.level, probability)
                for probability in probabilities
            ]
            if args.transform_mode != "polynomial"
            else None
        )
        if shell_coefficients is not None:
            for weight, value in enumerate(shell_coefficients):
                require_finite(
                    value, f"order {order}, coefficient weight {weight}"
                )
        alternating = arb(0)
        convolution_rows = []
        if args.transform_mode == "batched":
            radial_functions = arb_mat(args.message_bits + 1, order + 1)
            for weight in range(args.message_bits + 1):
                power = arb(1)
                for bias_power in range(order + 1):
                    radial_functions[weight, bias_power] = (
                        shell_coefficients[weight] * power
                    )
                    power *= biases[weight]
            fouriers = transform * radial_functions
            convolution_powers = range(order, -1, -1)
        elif args.transform_mode == "polynomial":
            remainder_bits = args.output_bits - order
            dyadic_scale = arb(2) ** (args.message_bits - remainder_bits)
            reversed_coefficients = arb_poly(
                [
                    dyadic_scale
                    * math.comb(remainder_bits, remainder_bits - index)
                    * krawtchouk(
                        args.output_bits,
                        args.level,
                        order + remainder_bits - index,
                    )
                    for index in range(remainder_bits + 1)
                ]
            )
            polynomial_fouriers = [
                reversed_coefficients * walk_polynomial
                for walk_polynomial in walk_polynomials
            ]
            for frequency_weight, polynomial_fourier in enumerate(
                polynomial_fouriers
            ):
                for bias_power in range(order + 1):
                    require_finite(
                        polynomial_fourier[remainder_bits + bias_power],
                        (
                            f"order {order}, bias power {bias_power}, "
                            f"polynomial Fourier weight {frequency_weight}"
                        ),
                    )
            fouriers = None
            convolution_powers = range(order + 1)
        else:
            fouriers = None
            convolution_powers = range(order + 1)
        for convolution_power in convolution_powers:
            bias_power = order - convolution_power
            if args.transform_mode == "polynomial":
                fourier = None
                fourier_column = None
            elif fouriers is None:
                radial_function = arb_mat(
                    [
                        [
                            shell_coefficients[weight]
                            * biases[weight] ** bias_power
                        ]
                        for weight in range(args.message_bits + 1)
                    ]
                )
                fourier = transform * radial_function
                fourier_column = 0
            else:
                fourier = fouriers
                fourier_column = bias_power
            if args.transform_mode != "polynomial":
                for frequency_weight in range(args.message_bits + 1):
                    require_finite(
                        fourier[frequency_weight, fourier_column],
                        f"order {order}, convolution power {convolution_power}, Fourier weight {frequency_weight}",
                    )
            quadratic = arb(0)
            for frequency_weight in range(args.message_bits + 1):
                walk_value = walk[convolution_power][frequency_weight]
                require_finite(
                    walk_value,
                    f"order {order}, convolution power {convolution_power}, walk weight {frequency_weight}",
                )
                fourier_value = (
                    polynomial_fouriers[frequency_weight][
                        remainder_bits + bias_power
                    ]
                    if args.transform_mode == "polynomial"
                    else fourier[frequency_weight, fourier_column]
                )
                # python-flint's generic Arb ``x ** 2`` path can return NaN
                # when x is a zero-centered interval. Direct multiplication
                # computes the same square and preserves a finite enclosure.
                square = fourier_value * fourier_value
                summand = (
                    multiplicities[frequency_weight]
                    * walk_value
                    * square
                )
                require_finite(
                    summand,
                    (
                        f"order {order}, convolution power {convolution_power}, "
                        f"quadratic weight {frequency_weight}; walk={walk_value}; "
                        f"fourier={fourier_value}"
                    ),
                )
                quadratic += summand
                require_finite(
                    quadratic,
                    f"order {order}, convolution power {convolution_power}, quadratic prefix {frequency_weight}",
                )
            require_finite(
                quadratic,
                f"order {order}, convolution power {convolution_power}, quadratic",
            )
            signed = (
                (-1) ** (order - convolution_power)
                * math.comb(order, convolution_power)
                * quadratic
            )
            alternating += signed
            require_finite(
                alternating,
                f"order {order}, convolution power {convolution_power}, alternating sum",
            )
            convolution_rows.append(
                {
                    "convolution_power": convolution_power,
                    "quadratic": str(quadratic),
                    "signed_term": str(signed),
                }
            )
        aggregate = (
            math.comb(args.output_bits, order)
            * alternating
            / 4**order
        )
        require_finite(aggregate, f"order {order}, aggregate")
        if aggregate < 0:
            raise ArithmeticError(
                f"order {order} contradicts proved nonnegativity: {aggregate}"
            )
        nonnegative_aggregate = (aggregate + abs(aggregate)) / 2
        partial += nonnegative_aggregate
        scaled = (
            (arb(2) ** args.message_bits)
            * nonnegative_aggregate
            / math.comb(args.output_bits, args.level)
        )
        order_rows.append(
            {
                "order": order,
                "aggregate": str(aggregate),
                "nonnegative_aggregate": str(nonnegative_aggregate),
                "scaled_diagonal_contribution": str(scaled),
                "convolution_terms": (
                    convolution_rows if args.include_convolution_terms else None
                ),
            }
        )
        if args.progress_every and (
            order % args.progress_every == 0 or order == args.maximum_order
        ):
            print(
                f"order,{order},{args.maximum_order},scaled,{scaled},elapsed_seconds,{time.perf_counter()-started:.3f}",
                flush=True,
            )
    scaled_partial = (
        (arb(2) ** args.message_bits)
        * partial
        / math.comb(args.output_bits, args.level)
    )
    payload = {
        "schema": "pure-ea-sector-zero-hoeffding-radial-arb-v1",
        "status": (
            "OUTWARD_FIXED_ORDER_PARTIAL_SUM"
            if args.minimum_order == 1
            else "OUTWARD_FIXED_ORDER_RANGE_SUM"
        ),
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "level": args.level,
            "minimum_order": args.minimum_order,
            "maximum_order": args.maximum_order,
            "precision_bits": args.precision,
            "transform_mode": args.transform_mode,
        },
        "orders": order_rows,
        "scaled_diagonal_range_sum": str(scaled_partial),
        "scaled_diagonal_partial_sum": (
            str(scaled_partial) if args.minimum_order == 1 else None
        ),
        "elapsed_seconds": time.perf_counter() - started,
        "scope": [
            "Each fixed-order aggregate is evaluated by the exact radial identity with Arb intervals.",
            "The Arb interval encloses the exact sum over the reported order range.",
            "A range starting at order one is a lower bound on the complete sector-zero diagonal because omitted orders are nonnegative.",
            "No upper bound or variance certificate follows until every omitted order is included or bounded.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if args.suppress_json:
        print(
            f"receipt,{args.output},range,{args.minimum_order},{args.maximum_order},scaled,{scaled_partial}",
            flush=True,
        )
    else:
        print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

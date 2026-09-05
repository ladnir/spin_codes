#!/usr/bin/env python3
"""Exact-rational counterexample to iid-lane entrywise domination.

The candidate is MultiBlockStripeFresh-32 with one 256-bit checkpoint epoch
per lane slice.  Compare two fair candidates:

* shared group: the active block lanes are {0,1}, followed by a uniform cyclic
  shift in each outer-coordinate region;
* iid model: the two candidates choose their lanes independently and uniformly.

At the rational Chernoff variable z=63/64, this script computes both (0,1)
entries exactly and verifies that the shared-group entry is larger.  Therefore
the iid region matrix does not dominate every shared-group matrix entrywise.
"""

from __future__ import annotations

import json
import math
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path


OUTPUT = Path(
    "constructions/riffle_multiblockstripe_fieldcheckpoint/"
    "receipts/shared_group_adjacent_pair_exact_counterexample.json"
)

Matrix = tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[row][column] + right[row][column] for column in range(2))
        for row in range(2)
    )  # type: ignore[return-value]


def matrix_scale(matrix: Matrix, scalar: Fraction) -> Matrix:
    return tuple(
        tuple(matrix[row][column] * scalar for column in range(2))
        for row in range(2)
    )  # type: ignore[return-value]


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum(left[row][inner] * right[inner][column] for inner in range(2))
            for column in range(2)
        )
        for row in range(2)
    )  # type: ignore[return-value]


def truncated_power(base: list[Fraction], exponent: int, degree: int) -> list[Fraction]:
    result = [Fraction(0) for _ in range(degree + 1)]
    result[0] = Fraction(1)
    for _ in range(exponent):
        updated = [Fraction(0) for _ in range(degree + 1)]
        for left_degree in range(degree + 1):
            for right_degree in range(degree - left_degree + 1):
                updated[left_degree + right_degree] += (
                    result[left_degree] * base[right_degree]
                )
        result = updated
    return result


def lane_prefix_polynomials(visits: int, z: Fraction) -> tuple[list[Fraction], list[Fraction]]:
    start_zero = [Fraction(0) for _ in range(visits + 1)]
    start_one = [Fraction(0) for _ in range(visits + 1)]
    for support in range(1 << visits):
        parity = 0
        prefix_weight = 0
        for visit in range(visits):
            parity ^= (support >> visit) & 1
            prefix_weight += parity
        impulses = support.bit_count()
        start_zero[impulses] += z**prefix_weight
        start_one[impulses] += z ** (visits - prefix_weight)
    return start_zero, start_one


def impulse_epoch_matrices(z: Fraction, maximum_impulses: int = 2) -> list[Matrix]:
    state_bits = 64
    epoch_bits = 256
    visits = epoch_bits // state_bits
    start_zero, start_one = lane_prefix_polynomials(visits, z)
    even_zero = [value if degree % 2 == 0 else Fraction(0) for degree, value in enumerate(start_zero)]
    odd_one = [value if degree % 2 == 1 else Fraction(0) for degree, value in enumerate(start_one)]
    lane_scale = Fraction(1, 1 << visits)

    def scaled(values: list[Fraction], scale: Fraction) -> list[Fraction]:
        return [value * scale for value in values[: maximum_impulses + 1]]

    zero_power = truncated_power(scaled(start_zero, lane_scale), state_bits, maximum_impulses)
    even_power = truncated_power(scaled(even_zero, lane_scale), state_bits, maximum_impulses)
    endpoint_power = truncated_power(
        scaled([even_zero[i] + odd_one[i] for i in range(visits + 1)], lane_scale),
        state_bits,
        maximum_impulses,
    )
    all_power = truncated_power(
        scaled(
            [start_zero[i] + start_one[i] for i in range(visits + 1)],
            lane_scale / 2,
        ),
        state_bits,
        maximum_impulses,
    )

    live_states = (1 << state_bits) - 1
    matrices = []
    for impulses in range(maximum_impulses + 1):
        support_probability = Fraction(math.comb(epoch_bits, impulses), 1 << epoch_bits)
        zero_to_zero = even_power[impulses] / support_probability
        zero_total = zero_power[impulses] / support_probability
        live_to_zero = (
            endpoint_power[impulses] - even_power[impulses]
        ) / (live_states * support_probability)
        live_total = (
            (1 << state_bits) * all_power[impulses] - zero_power[impulses]
        ) / (live_states * support_probability)
        matrices.append(
            (
                (zero_to_zero, zero_total - zero_to_zero),
                (live_to_zero, live_total - live_to_zero),
            )
        )
    return matrices


def fair_candidate_matrices(impulses: list[Matrix]) -> list[Matrix]:
    return [
        impulses[0],
        matrix_scale(matrix_add(impulses[0], impulses[1]), Fraction(1, 2)),
        matrix_scale(
            matrix_add(matrix_add(impulses[0], matrix_scale(impulses[1], 2)), impulses[2]),
            Fraction(1, 4),
        ),
    ]


def pattern_product(lanes: list[Matrix], active: set[int], q: int) -> Matrix:
    product: Matrix = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)))
    for lane in range(q):
        product = matrix_multiply(product, lanes[1 if lane in active else 0])
    return product


def adjacent_shared_transfer(lanes: list[Matrix], q: int) -> Matrix:
    total: Matrix = ((Fraction(0), Fraction(0)), (Fraction(0), Fraction(0)))
    for shift in range(q):
        active = {shift, (shift + 1) % q}
        total = matrix_add(total, pattern_product(lanes, active, q))
    return matrix_scale(total, Fraction(1, q))


def iid_pair_transfer(lanes: list[Matrix], q: int) -> Matrix:
    # Exponential-generating-function coefficient through degree two.
    coefficients: list[Matrix] = [
        ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1))),
        ((Fraction(0), Fraction(0)), (Fraction(0), Fraction(0))),
        ((Fraction(0), Fraction(0)), (Fraction(0), Fraction(0))),
    ]
    weighted = [lanes[0], lanes[1], matrix_scale(lanes[2], Fraction(1, 2))]
    for _ in range(q):
        updated: list[Matrix] = [
            ((Fraction(0), Fraction(0)), (Fraction(0), Fraction(0)))
            for _ in range(3)
        ]
        for total in range(3):
            for here in range(total + 1):
                updated[total] = matrix_add(
                    updated[total],
                    matrix_multiply(coefficients[total - here], weighted[here]),
                )
        coefficients = updated
    return matrix_scale(coefficients[2], Fraction(2, q * q))


def main() -> None:
    q = 32
    z = Fraction(63, 64)
    lanes = fair_candidate_matrices(impulse_epoch_matrices(z))
    shared = adjacent_shared_transfer(lanes, q)
    iid = iid_pair_transfer(lanes, q)
    difference = shared[0][1] - iid[0][1]
    if difference <= 0:
        raise AssertionError("the proposed adjacent-pair counterexample is not positive")

    getcontext().prec = 50
    decimal_difference = Decimal(difference.numerator) / Decimal(difference.denominator)
    payload = {
        "schema": "riffle-multiblockstripe-shared-group-exact-counterexample-v1",
        "candidate": "Riffle MultiBlockStripeFresh-32 FieldCheckpoint",
        "comparison": "shared adjacent pair minus iid-lane pair, matrix entry (zero,live)",
        "q": q,
        "active_lane_pattern": [0, 1],
        "z": {"numerator": z.numerator, "denominator": z.denominator},
        "difference": {
            "sign": "positive",
            "decimal_50_digits": str(decimal_difference),
            "numerator_bits": difference.numerator.bit_length(),
            "denominator_bits": difference.denominator.bit_length(),
        },
        "conclusion": (
            "The iid-lane region matrix does not entrywise dominate every "
            "shared-group region matrix."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"difference,{decimal_difference}")
    print(f"numerator_bits,{difference.numerator.bit_length()}")
    print(f"denominator_bits,{difference.denominator.bit_length()}")
    print(f"wrote,{OUTPUT}")


if __name__ == "__main__":
    main()

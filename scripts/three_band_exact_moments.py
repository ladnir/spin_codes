#!/usr/bin/env python3
"""Exact/outward moment helpers for the fixed three-band proof.

For ``d`` independent uniform weight-``w`` subsets of ``[n]`` and rational
``q``, inclusion-exclusion gives the exact diagonal OR moment

    sum_t C(n,t) q^(n-t) (1-q)^t
          (C(n-t,w) / C(n,w))^d.

The cell-cap proof uses its ``d``-th root.  We round that root upward to a
dyadic rational by an integer comparison, so every subsequent greedy
transport calculation is an exact rational upper bound.
"""

from __future__ import annotations

import math
from fractions import Fraction


def cluster_diagonal_exact(
    columns: int, weight: int, cluster_size: int, pole: Fraction
) -> Fraction:
    """Return the exact fixed-weight diagonal OR moment."""

    if not (0 <= weight <= columns and cluster_size >= 1 and 0 < pole < 1):
        raise ValueError("invalid diagonal-moment arguments")
    denominator = math.comb(columns, weight)
    total = Fraction(0)
    for avoided in range(columns - weight + 1):
        total += (
            math.comb(columns, avoided)
            * pole ** (columns - avoided)
            * (1 - pole) ** avoided
            * Fraction(math.comb(columns - avoided, weight), denominator)
            ** cluster_size
        )
    return total


def _ceil_nth_root_ratio(numerator: int, denominator: int, degree: int) -> int:
    """Least integer m satisfying m**degree * denominator >= numerator."""

    if numerator == 0:
        return 0
    if numerator < 0 or denominator <= 0 or degree <= 0:
        raise ValueError("invalid rational root arguments")
    lower = 0
    upper = 1 << ((numerator.bit_length() + degree - 1) // degree)
    while upper**degree * denominator < numerator:
        upper <<= 1
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if middle**degree * denominator >= numerator:
            upper = middle
        else:
            lower = middle
    return upper


def dyadic_nth_root_upper(value: Fraction, degree: int, bits: int) -> Fraction:
    """Return an exact dyadic upper bound on ``value**(1/degree)``."""

    if value < 0 or degree <= 0 or bits <= 0:
        raise ValueError("invalid dyadic root arguments")
    scale = 1 << bits
    numerator = value.numerator * scale**degree
    root = _ceil_nth_root_ratio(numerator, value.denominator, degree)
    result = Fraction(root, scale)
    if result**degree < value:
        raise AssertionError("dyadic root was not rounded upward")
    if root and Fraction(root - 1, scale) ** degree >= value:
        raise AssertionError("dyadic root is not the least grid upper bound")
    return result


def factor_table_upper(
    columns: int,
    cluster_size: int,
    pole: Fraction,
    *,
    root_bits: int = 160,
) -> tuple[Fraction, ...]:
    """Outward dyadic table for the rooted diagonal cluster moments."""

    return tuple(
        dyadic_nth_root_upper(
            cluster_diagonal_exact(columns, weight, cluster_size, pole),
            cluster_size,
            root_bits,
        )
        for weight in range(columns + 1)
    )


def greedy_factor_bound_upper(
    states: list[tuple[int, int, int]],
    caps: list[int],
    diagonal_mass: dict[int, int],
    by_total: dict[int, list[int]],
    factors: tuple[
        tuple[Fraction, ...], tuple[Fraction, ...], tuple[Fraction, ...]
    ],
) -> Fraction:
    """Exact upper version of the cell-cap greedy transport relaxation."""

    total = Fraction(0)
    for weight, mass in diagonal_mass.items():
        if not mass:
            continue
        scored = sorted(
            (
                factors[0][states[index][0]]
                * factors[1][states[index][1]]
                * factors[2][states[index][2]],
                caps[index],
            )
            for index in by_total[weight]
            if caps[index]
        )
        remaining = mass
        diagonal_value = Fraction(0)
        for value, cap in reversed(scored):
            take = min(remaining, cap)
            diagonal_value += take * value
            remaining -= take
            if not remaining:
                break
        if remaining:
            raise ValueError(f"cell caps miss total-weight mass {weight}")
        total += diagonal_value
    return total

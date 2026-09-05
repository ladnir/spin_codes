#!/usr/bin/env python3
"""Generate and verify binary degree-6/3 EC certificates near memory 80.

Exact uniform-slice blocks cover small supports.  Positive-coefficient blocks
cover both outer ranges.  The degree-three conditional point-mass bound covers
the central range.  Generation uses floating point only to select fixed
decimal markers.  Verification reads those markers, uses entry-scaled
outward-rounded binary64 for the exact matrix recurrence, and evaluates the
remaining scalar bounds with Arb arithmetic.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Any

import numpy as np
from flint import arb, ctx

from binary_biregular_diagnostic import (
    degree_three_parity_shell_counts,
    geometric_blocks,
    saddle_biregular_block_logterm,
)
from binary_biregular_ec_certificate import (
    coefficient_block_arb,
    endpoint_mass,
    full_support_marker,
    full_support_term_arb,
)
from check_binary_biregular_ec_d6_m80_certificate import (
    SCHEMA,
    validate_certificate_structure,
)
from ea_certificate import decimal_marker, hamming_ball_bound_arb
from regular_ec_certificate import (
    uniform_slice_transfer_matrices_arb,
    zero_matrix,
)
from positive_matrix_upper import (
    ScaledUpperMatrix,
    add_scaled_upper,
    add_upper,
    align_upper,
    normalize_upper,
    power_scaled_upper,
    rational_upper,
    row_sum_upper,
    scale_scaled_upper,
    multiply_scaled_upper,
)
from positive_matrix_factor_upper import (
    FactorUpperMatrix,
    multiply_factor,
    normalize_factor,
    power_factor,
    row_sum_factor,
    scale_factor,
    sum_factors,
)
from positive_matrix_row_upper import (
    RowUpperMatrix,
    multiply_rows,
    normalize_rows,
    power_rows,
    row_sum_rows,
    scale_rows,
    sum_rows,
)
from positive_matrix_layered_upper import (
    LayeredUpperMatrix,
    layered_from_matrix,
    layered_identity,
    multiply_layered,
    multiply_scale_layered_rational,
    power_layered,
    row_sum_layers,
    scale_layered_rational,
    sum_layered,
)


DEFAULT_EXACT_BLOCKS = (
    (1, 2, "0.9999"),
    (3, 7, "0.99985"),
    (8, 15, "0.99965"),
    (16, 31, "0.99925"),
    (32, 47, "0.99895"),
    (48, 63, "0.99855"),
    (64, 79, "0.99815"),
    (80, 95, "0.99775366"),
    (96, 111, "0.99732763"),
    (112, 127, "0.99691451"),
    (128, 143, "0.99650139"),
    (144, 159, "0.99608827"),
    (160, 175, "0.99567515"),
    (176, 191, "0.99526203"),
    (192, 223, "0.99464135"),
    (224, 255, "0.99381531"),
    (256, 383, "0.992"),
)


@dataclass(frozen=True)
class DegreeSixECResult:
    success: bool
    total_bound: arb
    security_bits: arb
    exact_bound: arb
    outer_bound: arb
    central_bound: arb
    full_support_bound: arb
    largest_exact_block: tuple[int, int]
    largest_exact_block_bound: arb
    exact_blocks: int
    outer_blocks: int
    central_blocks: int


def exact_block_arb(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, lo: int, hi: int, output_marker: str,
) -> arb:
    """Sum one exact support block while reusing its slice family."""
    if right_degree != 3:
        raise ValueError("the direct shell formula requires right degree three")
    z = arb(output_marker)
    if not (z > 0 and z <= 1):
        raise ValueError("invalid exact output marker")
    slices = uniform_slice_transfer_matrices_arb(
        length=k // right_degree,
        max_weight=hi,
        output_marker=z,
        memory=memory,
    )
    total = arb(0)
    for support in range(lo, hi + 1):
        denominator, counts = degree_three_parity_shell_counts(
            left_vertices=k,
            support_size=support,
        )
        region = zero_matrix(memory + 1)
        for weight, count in counts.items():
            region += slices[weight] * (arb(count) / denominator)
        tail = endpoint_mass(region, left_degree, memory) * z ** (-cutoff)
        total += arb(math.comb(k, support)) * tail
    return total


def wrapping_input_matrices_upper(
    output_marker: str, memory: int
) -> tuple[ScaledUpperMatrix, ScaledUpperMatrix]:
    """Return directed upper matrices for one zero or one convolution input."""
    nearest = float(Decimal(output_marker))
    z_upper = math.nextafter(nearest, math.inf)
    size = memory + 1
    zero = np.zeros((size, size), dtype=np.float64)
    one = np.zeros((size, size), dtype=np.float64)
    for state in range(memory - 1):
        zero[state, 0] = one[state, 0] = z_upper / 2.0
        zero[state, state + 1] = one[state, state + 1] = 0.5
    zero[memory - 1, 0] = z_upper
    one[memory - 1, memory] = 1.0
    zero[memory, memory] = 1.0
    one[memory, 0] = z_upper
    return (
        normalize_upper(ScaledUpperMatrix(zero, 0)),
        normalize_upper(ScaledUpperMatrix(one, 0)),
    )


def combine_uniform_slice_transfers_upper(
    left: tuple[int, list[ScaledUpperMatrix]],
    right: tuple[int, list[ScaledUpperMatrix]],
    max_weight: int,
) -> tuple[int, list[ScaledUpperMatrix]]:
    """Concatenate two families of directed upper slice matrices."""
    left_length, left_slices = left
    right_length, right_slices = right
    total_length = left_length + right_length
    maximum = min(max_weight, total_length)
    left_choose = [
        math.comb(left_length, weight)
        for weight in range(min(max_weight, left_length) + 1)
    ]
    right_choose = [
        math.comb(right_length, weight)
        for weight in range(min(max_weight, right_length) + 1)
    ]
    total_choose = [math.comb(total_length, weight) for weight in range(maximum + 1)]
    result: list[ScaledUpperMatrix] = []
    for weight in range(maximum + 1):
        accumulator: ScaledUpperMatrix | None = None
        denominator = total_choose[weight]
        lo = max(0, weight - right_length)
        hi = min(weight, left_length, len(left_slices) - 1)
        for left_weight in range(lo, hi + 1):
            right_weight = weight - left_weight
            if right_weight >= len(right_slices):
                continue
            probability = rational_upper(
                left_choose[left_weight] * right_choose[right_weight],
                denominator,
            )
            product = multiply_scaled_upper(
                left_slices[left_weight], right_slices[right_weight]
            )
            term = scale_scaled_upper(product, probability)
            if accumulator is None:
                accumulator = term
            else:
                exponent = max(accumulator.exponent, term.exponent)
                accumulator = ScaledUpperMatrix(
                    add_upper(
                        align_upper(accumulator, exponent),
                        align_upper(term, exponent),
                    ),
                    exponent,
                )
        if accumulator is None:
            raise ArithmeticError("empty uniform-slice convolution")
        result.append(normalize_upper(accumulator))
    return total_length, result


def uniform_slice_transfer_matrices_upper(
    *, length: int, max_weight: int, output_marker: str, memory: int,
) -> list[ScaledUpperMatrix]:
    """Return directed upper conditional slices through ``max_weight``."""
    zero, one = wrapping_input_matrices_upper(output_marker, memory)
    identity = ScaledUpperMatrix(np.eye(memory + 1, dtype=np.float64), 0)
    result: tuple[int, list[ScaledUpperMatrix]] = (0, [identity])
    power: tuple[int, list[ScaledUpperMatrix]] = (1, [zero, one])
    remaining = length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers_upper(
                result, power, max_weight
            )
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers_upper(
                power, power, max_weight
            )
    return result[1]


def arb_from_binary64(value: float) -> arb:
    """Convert a finite binary64 number to the same exact Arb value."""
    numerator, denominator = value.as_integer_ratio()
    return arb(numerator) / denominator


def exact_block_upper(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, lo: int, hi: int, output_marker: str,
) -> arb:
    """Bound an exact block using directed upper matrices and Arb scalars."""
    if right_degree != 3:
        raise ValueError("the direct shell formula requires right degree three")
    z = arb(output_marker)
    if not (z > 0 and z <= 1):
        raise ValueError("invalid exact output marker")
    slices = uniform_slice_transfer_matrices_upper(
        length=k // right_degree,
        max_weight=hi,
        output_marker=output_marker,
        memory=memory,
    )
    total = arb(0)
    for support in range(lo, hi + 1):
        denominator, counts = degree_three_parity_shell_counts(
            left_vertices=k,
            support_size=support,
        )
        region: ScaledUpperMatrix | None = None
        for weight, count in counts.items():
            term = scale_scaled_upper(
                slices[weight], rational_upper(count, denominator)
            )
            region = add_scaled_upper(region, term)
        if region is None:
            raise ArithmeticError("empty degree-three parity shell")
        transform = power_scaled_upper(region, left_degree)
        mantissa, exponent = row_sum_upper(transform, memory)
        tail_upper = arb_from_binary64(mantissa) * arb(2) ** exponent
        total += (
            arb(math.comb(k, support))
            * tail_upper
            * z ** (-cutoff)
        )
    return total


def wrapping_input_matrices_layered(
    output_marker: str, memory: int
) -> tuple[LayeredUpperMatrix, LayeredUpperMatrix]:
    """Return layered upper matrices for one zero or one convolution input."""
    z_upper = math.nextafter(float(Decimal(output_marker)), math.inf)
    size = memory + 1
    zero = np.zeros((size, size), dtype=np.float64)
    one = np.zeros((size, size), dtype=np.float64)
    for state in range(memory - 1):
        zero[state, 0] = one[state, 0] = z_upper / 2.0
        zero[state, state + 1] = one[state, state + 1] = 0.5
    zero[memory - 1, 0] = z_upper
    one[memory - 1, memory] = 1.0
    zero[memory, memory] = 1.0
    one[memory, 0] = z_upper
    return layered_from_matrix(zero), layered_from_matrix(one)


def combine_uniform_slice_transfers_layered(
    left: tuple[int, list[LayeredUpperMatrix]],
    right: tuple[int, list[LayeredUpperMatrix]],
    max_weight: int,
) -> tuple[int, list[LayeredUpperMatrix]]:
    """Concatenate two conditional layered slice families."""
    left_length, left_slices = left
    right_length, right_slices = right
    total_length = left_length + right_length
    maximum = min(max_weight, total_length)
    left_choose = [
        math.comb(left_length, weight)
        for weight in range(min(max_weight, left_length) + 1)
    ]
    right_choose = [
        math.comb(right_length, weight)
        for weight in range(min(max_weight, right_length) + 1)
    ]
    total_choose = [
        math.comb(total_length, weight) for weight in range(maximum + 1)
    ]
    result: list[LayeredUpperMatrix] = []
    for weight in range(maximum + 1):
        terms: list[LayeredUpperMatrix] = []
        lo = max(0, weight - right_length)
        hi = min(weight, left_length, len(left_slices) - 1)
        for left_weight in range(lo, hi + 1):
            right_weight = weight - left_weight
            if right_weight >= len(right_slices):
                continue
            terms.append(multiply_scale_layered_rational(
                left_slices[left_weight],
                right_slices[right_weight],
                left_choose[left_weight] * right_choose[right_weight],
                total_choose[weight],
            ))
        result.append(sum_layered(terms))
    return total_length, result


def uniform_slice_transfer_matrices_layered(
    *, length: int, max_weight: int, output_marker: str, memory: int,
) -> list[LayeredUpperMatrix]:
    """Return entry-scaled conditional slices through ``max_weight``."""
    zero, one = wrapping_input_matrices_layered(output_marker, memory)
    result: tuple[int, list[LayeredUpperMatrix]] = (
        0, [layered_identity(memory + 1)]
    )
    power: tuple[int, list[LayeredUpperMatrix]] = (1, [zero, one])
    remaining = length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers_layered(
                result, power, max_weight
            )
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers_layered(
                power, power, max_weight
            )
    return result[1]


def exact_block_layered_upper(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, lo: int, hi: int, output_marker: str,
) -> arb:
    """Bound an exact block with independently scaled entry layers."""
    if right_degree != 3:
        raise ValueError("the direct shell formula requires right degree three")
    z = arb(output_marker)
    if not (z > 0 and z <= 1):
        raise ValueError("invalid exact output marker")
    slices = uniform_slice_transfer_matrices_layered(
        length=k // right_degree,
        max_weight=hi,
        output_marker=output_marker,
        memory=memory,
    )
    total = arb(0)
    for support in range(lo, hi + 1):
        denominator, counts = degree_three_parity_shell_counts(
            left_vertices=k, support_size=support
        )
        region = sum_layered([
            scale_layered_rational(slices[weight], count, denominator)
            for weight, count in counts.items()
        ])
        transform = power_layered(region, left_degree)
        tail_upper = sum((
            arb_from_binary64(mantissa) * arb(2) ** exponent
            for mantissa, exponent in row_sum_layers(transform, memory)
        ), arb(0))
        total += (
            arb(math.comb(k, support))
            * tail_upper
            * z ** (-cutoff)
        )
    return total


def wrapping_input_matrices_factor(
    output_marker: str, memory: int
) -> tuple[FactorUpperMatrix, FactorUpperMatrix]:
    """Return factor-enclosed matrices for one convolution input bit."""
    z_upper = math.nextafter(float(Decimal(output_marker)), math.inf)
    size = memory + 1
    zero = np.zeros((size, size), dtype=np.float64)
    one = np.zeros((size, size), dtype=np.float64)
    for state in range(memory - 1):
        zero[state, 0] = one[state, 0] = z_upper / 2.0
        zero[state, state + 1] = one[state, state + 1] = 0.5
    zero[memory - 1, 0] = z_upper
    one[memory - 1, memory] = 1.0
    zero[memory, memory] = 1.0
    one[memory, 0] = z_upper
    return (
        normalize_factor(FactorUpperMatrix(zero, 0)),
        normalize_factor(FactorUpperMatrix(one, 0)),
    )


def combine_uniform_slice_transfers_factor(
    left: tuple[int, list[FactorUpperMatrix]],
    right: tuple[int, list[FactorUpperMatrix]],
    max_weight: int,
) -> tuple[int, list[FactorUpperMatrix]]:
    """Concatenate two conditional slice families with global factors."""
    left_length, left_slices = left
    right_length, right_slices = right
    total_length = left_length + right_length
    maximum = min(max_weight, total_length)
    left_choose = [
        math.comb(left_length, weight)
        for weight in range(min(max_weight, left_length) + 1)
    ]
    right_choose = [
        math.comb(right_length, weight)
        for weight in range(min(max_weight, right_length) + 1)
    ]
    total_choose = [math.comb(total_length, weight) for weight in range(maximum + 1)]
    result: list[FactorUpperMatrix] = []
    for weight in range(maximum + 1):
        terms: list[FactorUpperMatrix] = []
        lo = max(0, weight - right_length)
        hi = min(weight, left_length, len(left_slices) - 1)
        for left_weight in range(lo, hi + 1):
            right_weight = weight - left_weight
            if right_weight >= len(right_slices):
                continue
            probability = rational_upper(
                left_choose[left_weight] * right_choose[right_weight],
                total_choose[weight],
            )
            terms.append(scale_factor(
                multiply_factor(
                    left_slices[left_weight], right_slices[right_weight]
                ),
                probability,
            ))
        result.append(sum_factors(terms))
    return total_length, result


def uniform_slice_transfer_matrices_factor(
    *, length: int, max_weight: int, output_marker: str, memory: int,
) -> list[FactorUpperMatrix]:
    """Return scaled conditional slices with global error factors."""
    zero, one = wrapping_input_matrices_factor(output_marker, memory)
    result: tuple[int, list[FactorUpperMatrix]] = (
        0, [FactorUpperMatrix(np.eye(memory + 1), 0)]
    )
    power: tuple[int, list[FactorUpperMatrix]] = (1, [zero, one])
    remaining = length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers_factor(
                result, power, max_weight
            )
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers_factor(
                power, power, max_weight
            )
    return result[1]


def exact_block_factor_upper(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, lo: int, hi: int, output_marker: str,
) -> arb:
    """Bound an exact block with global matrix factors and Arb scalars."""
    if right_degree != 3:
        raise ValueError("the direct shell formula requires right degree three")
    z = arb(output_marker)
    slices = uniform_slice_transfer_matrices_factor(
        length=k // right_degree,
        max_weight=hi,
        output_marker=output_marker,
        memory=memory,
    )
    total = arb(0)
    for support in range(lo, hi + 1):
        denominator, counts = degree_three_parity_shell_counts(
            left_vertices=k, support_size=support
        )
        region = sum_factors([
            scale_factor(slices[weight], rational_upper(count, denominator))
            for weight, count in counts.items()
        ])
        transform = power_factor(region, left_degree)
        mantissa, exponent, factor, error = row_sum_factor(
            transform, memory
        )
        tail_upper = (
            (
                arb_from_binary64(mantissa)
                * arb_from_binary64(factor)
                + arb_from_binary64(error)
            )
            * arb(2) ** exponent
        )
        total += (
            arb(math.comb(k, support))
            * tail_upper
            * z ** (-cutoff)
        )
    return total


def wrapping_input_matrices_rows(
    output_marker: str, memory: int
) -> tuple[RowUpperMatrix, RowUpperMatrix]:
    """Return row-scaled matrices for one convolution input bit."""
    z_upper = math.nextafter(float(Decimal(output_marker)), math.inf)
    size = memory + 1
    zero = np.zeros((size, size), dtype=np.float64)
    one = np.zeros((size, size), dtype=np.float64)
    for state in range(memory - 1):
        zero[state, 0] = one[state, 0] = z_upper / 2.0
        zero[state, state + 1] = one[state, state + 1] = 0.5
    zero[memory - 1, 0] = z_upper
    one[memory - 1, memory] = 1.0
    zero[memory, memory] = 1.0
    one[memory, 0] = z_upper
    exponents = np.zeros(size, dtype=np.int64)
    return (
        normalize_rows(RowUpperMatrix(zero, exponents.copy())),
        normalize_rows(RowUpperMatrix(one, exponents.copy())),
    )


def combine_uniform_slice_transfers_rows(
    left: tuple[int, list[RowUpperMatrix]],
    right: tuple[int, list[RowUpperMatrix]],
    max_weight: int,
) -> tuple[int, list[RowUpperMatrix]]:
    """Concatenate two conditional row-scaled slice families."""
    left_length, left_slices = left
    right_length, right_slices = right
    total_length = left_length + right_length
    maximum = min(max_weight, total_length)
    left_choose = [
        math.comb(left_length, weight)
        for weight in range(min(max_weight, left_length) + 1)
    ]
    right_choose = [
        math.comb(right_length, weight)
        for weight in range(min(max_weight, right_length) + 1)
    ]
    total_choose = [math.comb(total_length, weight) for weight in range(maximum + 1)]
    result: list[RowUpperMatrix] = []
    for weight in range(maximum + 1):
        terms: list[RowUpperMatrix] = []
        lo = max(0, weight - right_length)
        hi = min(weight, left_length, len(left_slices) - 1)
        for left_weight in range(lo, hi + 1):
            right_weight = weight - left_weight
            if right_weight >= len(right_slices):
                continue
            probability = rational_upper(
                left_choose[left_weight] * right_choose[right_weight],
                total_choose[weight],
            )
            terms.append(scale_rows(
                multiply_rows(
                    left_slices[left_weight], right_slices[right_weight]
                ),
                probability,
            ))
        result.append(sum_rows(terms))
    return total_length, result


def uniform_slice_transfer_matrices_rows(
    *, length: int, max_weight: int, output_marker: str, memory: int,
) -> list[RowUpperMatrix]:
    """Return row-scaled conditional slices through ``max_weight``."""
    zero, one = wrapping_input_matrices_rows(output_marker, memory)
    size = memory + 1
    result: tuple[int, list[RowUpperMatrix]] = (
        0,
        [RowUpperMatrix(
            np.eye(size, dtype=np.float64), np.zeros(size, dtype=np.int64)
        )],
    )
    power: tuple[int, list[RowUpperMatrix]] = (1, [zero, one])
    remaining = length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers_rows(
                result, power, max_weight
            )
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers_rows(
                power, power, max_weight
            )
    return result[1]


def exact_block_rows_upper(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, lo: int, hi: int, output_marker: str,
) -> arb:
    """Bound an exact block with row-scaled binary64 and Arb scalars."""
    if right_degree != 3:
        raise ValueError("the direct shell formula requires right degree three")
    z = arb(output_marker)
    if not (z > 0 and z <= 1):
        raise ValueError("invalid exact output marker")
    slices = uniform_slice_transfer_matrices_rows(
        length=k // right_degree,
        max_weight=hi,
        output_marker=output_marker,
        memory=memory,
    )
    total = arb(0)
    for support in range(lo, hi + 1):
        denominator, counts = degree_three_parity_shell_counts(
            left_vertices=k, support_size=support
        )
        region = sum_rows([
            scale_rows(slices[weight], rational_upper(count, denominator))
            for weight, count in counts.items()
        ])
        transform = power_rows(region, left_degree)
        mantissa, exponent, factor, error = row_sum_rows(transform, memory)
        tail_upper = (
            (
                arb_from_binary64(mantissa)
                * arb_from_binary64(factor)
                + arb_from_binary64(error)
            )
            * arb(2) ** exponent
        )
        total += (
            arb(math.comb(k, support))
            * tail_upper
            * z ** (-cutoff)
        )
    return total


def degree_three_odd_variance_arb(p: arb) -> arb:
    """Return Var(Z_1) for a three-slot group conditioned on odd parity."""
    return 3 * p**2 * (1 - p) ** 2 / (4 * p**2 - 6 * p + 3) ** 2


def degree_three_dense_low_half_block_arb(
    *, k: int, left_degree: int, cutoff: int, lo: int, hi: int,
) -> arb:
    """Bound one degree-three central block in the low half."""
    region_length = k // 3
    n = left_degree * region_length
    p_lo = arb(lo) / k
    p_hi = arb(hi) / k
    q_zero_lo = (1 + (1 - 2 * p_lo) ** 3) / 2
    variance_lo = degree_three_odd_variance_arb(p_lo)
    local_mass = (arb.pi() / (8 * region_length * variance_lo)).sqrt()
    if not local_mass < 1:
        raise ValueError("degree-three local factor is not below one")
    choose_hi = arb(math.comb(k, hi))
    binomial_mass_hi = choose_hi * p_hi**hi * (1 - p_hi) ** (k - hi)
    regional_point_mass = (
        q_zero_lo**region_length * local_mass / binomial_mass_hi
    )
    return (
        arb(hi - lo + 1)
        * choose_hi
        * hamming_ball_bound_arb(n, cutoff)
        * regional_point_mass**left_degree
    )


def generate_certificate(
    *, k: int = 1_048_575, left_degree: int = 6,
    right_degree: int = 3, cutoff: int = 230_729, memory: int = 80,
    exact_limit: int = 383, low_outer_limit: int = 399_999,
    low_outer_relative_width: float = 0.1,
    high_outer_relative_width: float = 0.5,
    central_relative_width: float = 0.01,
    target_bits: int = 20, precision_bits: int = 192,
) -> dict[str, Any]:
    """Select markers and return a structurally complete certificate."""
    if (k, left_degree, right_degree) != (1_048_575, 6, 3):
        raise ValueError("this generator is specialized to the selected degree profile")
    if memory < 1:
        raise ValueError("memory must be positive")
    if exact_limit != DEFAULT_EXACT_BLOCKS[-1][1]:
        raise ValueError("exact_limit differs from the frozen block plan")
    n = left_degree * (k // right_degree)
    exact_blocks = [
        {"lo": lo, "hi": hi, "z": z}
        for lo, hi, z in DEFAULT_EXACT_BLOCKS
    ]

    low_ranges = geometric_blocks(
        exact_limit + 1, low_outer_limit, low_outer_relative_width
    )
    complement_ranges = geometric_blocks(
        1, low_outer_limit, high_outer_relative_width
    )
    high_ranges = [
        (k - hi, k - lo) for lo, hi in reversed(complement_ranges)
    ]
    if high_ranges[0][0] != k - low_outer_limit:
        raise AssertionError("high outer range does not follow the central complement")
    outer_blocks = []
    for index, (lo, hi) in enumerate(low_ranges + high_ranges, start=1):
        selected = saddle_biregular_block_logterm(
            code="ec",
            k=k,
            left_degree=left_degree,
            right_degree=right_degree,
            cutoff=cutoff,
            memory=memory,
            support_start=lo,
            support_limit=hi,
            multistart=(index == 1),
        )
        outer_blocks.append({
            "lo": lo,
            "hi": hi,
            "x": decimal_marker(selected.input_marker),
            "z": decimal_marker(selected.output_marker),
        })
        print(
            f"selected outer block {index}/{len(low_ranges) + len(high_ranges)} "
            f"{lo}..{hi}: log2={selected.log2_bound:.6f}",
            flush=True,
        )

    half = (k - 1) // 2
    central_blocks = [
        {"lo": lo, "hi": hi, "include_complement": True}
        for lo, hi in geometric_blocks(
            low_outer_limit + 1, half, central_relative_width
        )
    ]
    certificate = {
        "schema": SCHEMA,
        "parameters": {
            "k": k,
            "n": n,
            "cutoff": cutoff,
            "left_degree": left_degree,
            "right_degree": right_degree,
            "memory": memory,
            "target_bits": target_bits,
        },
        "verification": {"precision_bits": precision_bits},
        "exact_blocks": exact_blocks,
        "outer_blocks": outer_blocks,
        "central_blocks": central_blocks,
        "full_support": {
            "r": k,
            "z": full_support_marker(n=n, cutoff=cutoff, memory=memory),
        },
    }
    validate_certificate_structure(certificate)
    return certificate


def verify_certificate(certificate: dict[str, Any]) -> DegreeSixECResult:
    """Verify a frozen degree-six certificate with Arb arithmetic."""
    validate_certificate_structure(certificate)
    parameters = certificate["parameters"]
    k = int(parameters["k"])
    n = int(parameters["n"])
    cutoff = int(parameters["cutoff"])
    left_degree = int(parameters["left_degree"])
    right_degree = int(parameters["right_degree"])
    memory = int(parameters["memory"])
    target_bits = int(parameters["target_bits"])
    ctx.prec = int(certificate["verification"]["precision_bits"])

    exact_bounds = []
    for index, item in enumerate(certificate["exact_blocks"], start=1):
        bound = exact_block_layered_upper(
            k=k,
            left_degree=left_degree,
            right_degree=right_degree,
            cutoff=cutoff,
            memory=memory,
            lo=int(item["lo"]),
            hi=int(item["hi"]),
            output_marker=str(item["z"]),
        )
        exact_bounds.append(bound)
        print(
            f"verified exact block {index}/{len(certificate['exact_blocks'])} "
            f"{item['lo']}..{item['hi']}: {bound}",
            flush=True,
        )
    exact = sum(exact_bounds, arb(0))
    largest_exact_index = max(
        range(len(exact_bounds)), key=lambda index: exact_bounds[index]
    )

    outer = arb(0)
    for index, item in enumerate(certificate["outer_blocks"], start=1):
        outer += coefficient_block_arb(
            k=k,
            left_degree=left_degree,
            right_degree=right_degree,
            cutoff=cutoff,
            memory=memory,
            lo=int(item["lo"]),
            hi=int(item["hi"]),
            input_marker=str(item["x"]),
            output_marker=str(item["z"]),
        )
        if index == 1 or index % 20 == 0 or index == len(certificate["outer_blocks"]):
            print(
                f"verified outer block {index}/{len(certificate['outer_blocks'])}",
                flush=True,
            )

    central_low = sum((degree_three_dense_low_half_block_arb(
        k=k,
        left_degree=left_degree,
        cutoff=cutoff,
        lo=int(item["lo"]),
        hi=int(item["hi"]),
    ) for item in certificate["central_blocks"]), arb(0))
    central = 2 * central_low

    full_support = full_support_term_arb(
        n=n,
        cutoff=cutoff,
        memory=memory,
        output_marker=str(certificate["full_support"]["z"]),
    )
    total = exact + outer + central + full_support
    largest_item = certificate["exact_blocks"][largest_exact_index]
    return DegreeSixECResult(
        success=bool(total < arb(2) ** (-target_bits)),
        total_bound=total,
        security_bits=-total.log() / arb(2).log(),
        exact_bound=exact,
        outer_bound=outer,
        central_bound=central,
        full_support_bound=full_support,
        largest_exact_block=(int(largest_item["lo"]), int(largest_item["hi"])),
        largest_exact_block_bound=exact_bounds[largest_exact_index],
        exact_blocks=len(exact_bounds),
        outer_blocks=len(certificate["outer_blocks"]),
        central_blocks=len(certificate["central_blocks"]),
    )


def print_result(result: DegreeSixECResult) -> None:
    print(f"success: {result.success}")
    print(f"total bound: {result.total_bound}")
    print(f"security bits: {result.security_bits}")
    print(f"exact contribution: {result.exact_bound}")
    print(f"outer contribution: {result.outer_bound}")
    print(f"central contribution: {result.central_bound}")
    print(f"full-support contribution: {result.full_support_bound}")
    print(
        "largest exact block: "
        f"{result.largest_exact_block[0]}..{result.largest_exact_block[1]}, "
        f"{result.largest_exact_block_bound}"
    )
    print(f"exact blocks: {result.exact_blocks}")
    print(f"outer blocks: {result.outer_blocks}")
    print(f"central blocks: {result.central_blocks}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate")
    generate.add_argument("--output", type=Path, required=True)
    generate.add_argument("--memory", type=int, default=80)
    generate.add_argument("--target-bits", type=int, default=20)
    generate.add_argument("--precision-bits", type=int, default=192)
    verify = subparsers.add_parser("verify")
    verify.add_argument("certificate", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "generate":
        certificate = generate_certificate(
            memory=args.memory,
            target_bits=args.target_bits,
            precision_bits=args.precision_bits,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(certificate, indent=2) + "\n", encoding="utf-8"
        )
        print(f"wrote {args.output}")
        return
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    result = verify_certificate(certificate)
    print_result(result)
    if not result.success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

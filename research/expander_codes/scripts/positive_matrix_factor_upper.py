#!/usr/bin/env python3
"""Scaled nonnegative matrices with a rigorous row-sum error bound.

``FactorUpperMatrix(matrix, exponent, factor, error)`` represents the
entrywise upper bound ``2**exponent * (factor * matrix + error_matrix)`` for
some nonnegative ``error_matrix`` whose largest row sum is at most ``error``.
The relative factor handles normal rounding.  The row-sum error handles
gradual underflow without placing subnormals in structural zeros or paying a
matrix-dimension factor every time an error is propagated.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from positive_matrix_upper import MIN_SUBNORMAL, UNIT_ROUNDOFF, rational_upper


@dataclass(frozen=True)
class FactorUpperMatrix:
    matrix: np.ndarray
    exponent: int
    factor: float = 1.0
    error: float = 0.0


def multiply_factor_upper(*values: float) -> float:
    """Multiply positive upper factors with directed scalar rounding."""
    if any(value == 0.0 for value in values):
        return 0.0
    result = 1.0
    for value in values:
        result = math.nextafter(result * value, math.inf)
    if not math.isfinite(result):
        raise OverflowError("nonfinite matrix error factor")
    return result


def add_scalar_upper(*values: float) -> float:
    """Sum nonnegative error bounds with directed scalar rounding."""
    return math.nextafter(math.fsum(values), math.inf)


def operation_inflation(count: int) -> float:
    """Return a conservative upper factor for ``count`` binary64 operations."""
    if count < 1 or count * UNIT_ROUNDOFF >= 0.25:
        raise ValueError("invalid operation count")
    return math.nextafter(1.0 + 4.0 * count * UNIT_ROUNDOFF, math.inf)


def row_norm_upper(matrix: np.ndarray) -> float:
    """Upper-bound the maximum row sum of a nonnegative matrix."""
    columns = matrix.shape[1]
    nearest = float(np.max(np.sum(matrix, axis=1)))
    return multiply_factor_upper(
        add_scalar_upper(
            nearest,
            multiply_factor_upper(float(columns), MIN_SUBNORMAL),
        ),
        operation_inflation(columns),
    )


def normalize_factor(value: FactorUpperMatrix) -> FactorUpperMatrix:
    """Normalize by an exact power of two."""
    maximum = float(np.max(value.matrix))
    if maximum == 0.0:
        return value
    _, shift = math.frexp(maximum)
    matrix = np.ldexp(value.matrix, -shift)
    matrix = np.where(
        (value.matrix > 0.0) & (matrix == 0.0), MIN_SUBNORMAL, matrix
    )
    error = math.ldexp(value.error, -shift)
    if value.error > 0.0 and error == 0.0:
        error = MIN_SUBNORMAL
    return FactorUpperMatrix(
        matrix, value.exponent + shift, value.factor, error
    )


def scale_factor(
    value: FactorUpperMatrix, scalar_upper: float
) -> FactorUpperMatrix:
    """Multiply a matrix by a nonnegative upper scalar."""
    product = value.matrix * scalar_upper
    inflation = operation_inflation(1)
    factor = multiply_factor_upper(
        value.factor, inflation
    )
    columns = value.matrix.shape[1]
    scaled_error = multiply_factor_upper(value.error, scalar_upper)
    rounding_error = multiply_factor_upper(
        value.factor, float(columns), MIN_SUBNORMAL
    )
    error = add_scalar_upper(scaled_error, rounding_error)
    return normalize_factor(FactorUpperMatrix(
        product, value.exponent, factor, error
    ))


def multiply_factor(
    left: FactorUpperMatrix, right: FactorUpperMatrix
) -> FactorUpperMatrix:
    """Multiply two nonnegative matrices with one relative error factor."""
    inner = left.matrix.shape[1]
    product = left.matrix @ right.matrix
    inflation = operation_inflation(inner)
    factor = multiply_factor_upper(
        left.factor, right.factor, inflation
    )
    left_norm = row_norm_upper(left.matrix)
    right_norm = row_norm_upper(right.matrix)
    columns = right.matrix.shape[1]
    main_rounding_error = multiply_factor_upper(
        left.factor,
        right.factor,
        float(columns),
        float(inner),
        MIN_SUBNORMAL,
    )
    error = add_scalar_upper(
        multiply_factor_upper(left.factor, left_norm, right.error),
        multiply_factor_upper(left.error, right.factor, right_norm),
        multiply_factor_upper(left.error, right.error),
        main_rounding_error,
    )
    return normalize_factor(FactorUpperMatrix(
        product, left.exponent + right.exponent, factor, error
    ))


def sum_factors(values: list[FactorUpperMatrix]) -> FactorUpperMatrix:
    """Sum nonnegative scaled matrices with one common upper factor."""
    if not values:
        raise ValueError("cannot sum an empty matrix family")
    exponent = max(value.exponent for value in values)
    factor = max(value.factor for value in values)
    aligned = []
    errors = []
    for value in values:
        matrix = np.ldexp(value.matrix, value.exponent - exponent)
        aligned.append(matrix)
        error = math.ldexp(value.error, value.exponent - exponent)
        if value.error > 0.0 and error == 0.0:
            error = MIN_SUBNORMAL
        errors.append(error)
    matrix = np.sum(np.stack(aligned, axis=0), axis=0)
    inflation = operation_inflation(len(values))
    factor = multiply_factor_upper(
        factor, inflation
    )
    columns = matrix.shape[1]
    # Each alignment and the subsequent sum can lose less than one minimum
    # subnormal per entry and term.  Convert that entrywise allowance to a
    # maximum row-sum allowance.
    rounding_error = multiply_factor_upper(
        factor,
        2.0,
        float(len(values)),
        float(columns),
        MIN_SUBNORMAL,
    )
    error = add_scalar_upper(*errors, rounding_error)
    return normalize_factor(FactorUpperMatrix(
        matrix, exponent, factor, error
    ))


def power_factor(value: FactorUpperMatrix, exponent: int) -> FactorUpperMatrix:
    """Raise a square nonnegative matrix to an integer power."""
    if exponent < 0 or value.matrix.shape[0] != value.matrix.shape[1]:
        raise ValueError("invalid matrix power")
    result = FactorUpperMatrix(np.eye(value.matrix.shape[0]), 0, 1.0, 0.0)
    power = value
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply_factor(result, power)
        remaining >>= 1
        if remaining:
            power = multiply_factor(power, power)
    return result


def row_sum_factor(
    value: FactorUpperMatrix, row: int
) -> tuple[float, int, float, float]:
    """Return the components of a rigorous upper bound on one row sum."""
    count = value.matrix.shape[1]
    mantissa = float(np.sum(value.matrix[row, :]))
    inflation = operation_inflation(count)
    factor = multiply_factor_upper(
        value.factor, inflation
    )
    error = add_scalar_upper(
        value.error,
        multiply_factor_upper(
            value.factor, float(count), MIN_SUBNORMAL
        ),
    )
    return mantissa, value.exponent, factor, error


__all__ = [
    "FactorUpperMatrix",
    "multiply_factor",
    "normalize_factor",
    "power_factor",
    "rational_upper",
    "row_sum_factor",
    "scale_factor",
    "sum_factors",
]

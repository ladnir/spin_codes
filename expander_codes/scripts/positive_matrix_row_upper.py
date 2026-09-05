#!/usr/bin/env python3
"""Rigorous row-scaled arithmetic for nonnegative binary64 matrices.

For each row ``i``, ``RowUpperMatrix`` represents an entrywise upper bound

``A[i,:] <= 2**exponents[i] * (factor * matrix[i,:] + E[i,:])``,

where ``E`` is nonnegative and its row sum is at most ``errors[i]``.  A
separate power-of-two scale for every row keeps rare Markov states in the
normal binary64 range.  The residual vector accounts for gradual underflow.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from positive_matrix_upper import UNIT_ROUNDOFF, rational_upper


# Deliberately flush smaller values into the residual bound.  Keeping every
# stored operand normal avoids the severe hardware cost of subnormal BLAS.
SAFE_MIN = float.fromhex("0x1.0p-1022")


@dataclass(frozen=True)
class RowUpperMatrix:
    matrix: np.ndarray
    exponents: np.ndarray
    factor: float = 1.0
    errors: np.ndarray | None = None

    def __post_init__(self) -> None:
        rows = self.matrix.shape[0]
        if self.matrix.ndim != 2 or self.exponents.shape != (rows,):
            raise ValueError("invalid row-scaled matrix dimensions")
        if self.errors is None:
            object.__setattr__(self, "errors", np.zeros(rows, dtype=np.float64))
        elif self.errors.shape != (rows,):
            raise ValueError("invalid row-error dimensions")


def inflation(count: int) -> float:
    if count < 1 or count * UNIT_ROUNDOFF >= 0.25:
        raise ValueError("invalid operation count")
    return math.nextafter(1.0 + 4.0 * count * UNIT_ROUNDOFF, math.inf)


def product_upper(*values: float) -> float:
    if any(value == 0.0 for value in values):
        return 0.0
    result = 1.0
    for value in values:
        result = math.nextafter(result * value, math.inf)
    if not math.isfinite(result):
        raise OverflowError("nonfinite row error")
    return result


def sum_upper(*values: float) -> float:
    return math.nextafter(math.fsum(values), math.inf)


def array_upper(values: np.ndarray, count: int) -> np.ndarray:
    """Inflate a nonnegative vector expression for rounding and underflow."""
    result = values * inflation(count) + count * SAFE_MIN
    return np.nextafter(result, math.inf)


def row_norms_upper(matrix: np.ndarray) -> np.ndarray:
    """Upper-bound every row sum of a stored nonnegative matrix."""
    columns = matrix.shape[1]
    values = np.sum(matrix, axis=1)
    values = values * inflation(columns) + columns * SAFE_MIN
    return np.nextafter(values, math.inf)


def normalize_rows(value: RowUpperMatrix) -> RowUpperMatrix:
    """Normalize every row by an exact power of two."""
    assert value.errors is not None
    maxima = np.maximum(np.max(value.matrix, axis=1), value.errors)
    _, shifts = np.frexp(maxima)
    shifts = shifts.astype(np.int64)
    matrix = np.ldexp(value.matrix, -shifts[:, None])
    errors = np.ldexp(value.errors, -shifts)
    small_matrix = (matrix > 0.0) & (matrix < SAFE_MIN)
    matrix = np.where(small_matrix, 0.0, matrix)
    lost = np.count_nonzero(
        ((value.matrix > 0.0) & (matrix == 0.0)) | small_matrix, axis=1
    )
    lost_error = (value.errors > 0.0) & (errors < SAFE_MIN)
    errors = np.where(lost_error, 0.0, errors)
    errors = array_upper(
        errors
        + lost_error.astype(np.float64) * SAFE_MIN
        + value.factor * lost.astype(np.float64) * SAFE_MIN,
        5,
    )
    return RowUpperMatrix(
        matrix,
        value.exponents + shifts,
        value.factor,
        errors,
    )


def scale_rows(value: RowUpperMatrix, scalar: float) -> RowUpperMatrix:
    """Multiply by a nonnegative upper scalar."""
    if scalar < 0.0 or not math.isfinite(scalar):
        raise ValueError("invalid nonnegative scalar")
    assert value.errors is not None
    matrix = value.matrix * scalar
    factor = product_upper(value.factor, inflation(1))
    columns = value.matrix.shape[1]
    errors = array_upper(
        value.errors * scalar
        + value.factor * float(columns) * SAFE_MIN,
        5,
    )
    return normalize_rows(RowUpperMatrix(
        matrix, value.exponents.copy(), factor, errors
    ))


def sum_rows(values: list[RowUpperMatrix]) -> RowUpperMatrix:
    """Sum row-scaled nonnegative matrices."""
    if not values:
        raise ValueError("cannot sum an empty matrix family")
    rows, columns = values[0].matrix.shape
    exponents = np.maximum.reduce([value.exponents for value in values])
    factor0 = max(value.factor for value in values)
    aligned = []
    aligned_errors = []
    lost_counts = np.zeros(rows, dtype=np.int64)
    for value in values:
        assert value.errors is not None
        shifts = value.exponents - exponents
        matrix = np.ldexp(value.matrix, shifts[:, None])
        small_matrix = (matrix > 0.0) & (matrix < SAFE_MIN)
        matrix = np.where(small_matrix, 0.0, matrix)
        lost_counts += np.count_nonzero(
            ((value.matrix > 0.0) & (matrix == 0.0)) | small_matrix,
            axis=1,
        )
        aligned.append(matrix)
        errors = np.ldexp(value.errors, shifts)
        lost_error = (value.errors > 0.0) & (errors < SAFE_MIN)
        errors = np.where(lost_error, SAFE_MIN, errors)
        aligned_errors.append(errors)
    matrix = np.sum(np.stack(aligned, axis=0), axis=0)
    factor = product_upper(factor0, inflation(len(values)))
    errors = array_upper(
        np.sum(np.stack(aligned_errors, axis=0), axis=0)
        + factor
        * (2.0 * len(values) * columns + lost_counts.astype(np.float64))
        * SAFE_MIN,
        2 * len(values) + 5,
    )
    return normalize_rows(RowUpperMatrix(matrix, exponents, factor, errors))


def multiply_rows(left: RowUpperMatrix, right: RowUpperMatrix) -> RowUpperMatrix:
    """Multiply two row-scaled nonnegative matrices."""
    if left.matrix.shape[1] != right.matrix.shape[0]:
        raise ValueError("matrix dimensions do not match")
    assert left.errors is not None and right.errors is not None
    rows, inner = left.matrix.shape
    columns = right.matrix.shape[1]
    right_base = int(np.max(right.exponents))
    shifts = right.exponents - right_base
    scale = np.ldexp(np.ones(inner, dtype=np.float64), shifts)
    scaled_left = left.matrix * scale[None, :]
    small_left = (scaled_left > 0.0) & (scaled_left < SAFE_MIN)
    scaled_left = np.where(small_left, 0.0, scaled_left)
    lost_counts = np.count_nonzero(
        ((left.matrix > 0.0) & (scaled_left == 0.0)) | small_left, axis=1
    )
    matrix = scaled_left @ right.matrix
    exponents = left.exponents + right_base
    factor = product_upper(left.factor, right.factor, inflation(inner))

    right_norms = row_norms_upper(right.matrix)
    scaled_right_errors = np.ldexp(right.errors, shifts)
    scaled_right_errors = np.where(
        (right.errors > 0.0) & (scaled_right_errors < SAFE_MIN),
        SAFE_MIN,
        scaled_right_errors,
    )
    scaled_right_norms = np.ldexp(right_norms, shifts)
    scaled_right_norms = np.where(
        (right_norms > 0.0) & (scaled_right_norms < SAFE_MIN),
        SAFE_MIN,
        scaled_right_norms,
    )
    right_total = array_upper(
        right.factor * scaled_right_norms + scaled_right_errors, 3
    )
    maximum_right = math.nextafter(float(np.max(right_total)), math.inf)

    cross = scaled_left @ scaled_right_errors
    cross = cross * inflation(inner) + inner * SAFE_MIN
    cross = np.nextafter(cross, math.inf)
    errors = array_upper(
        left.factor * cross
        + left.errors * maximum_right
        + left.factor
        * lost_counts.astype(np.float64)
        * SAFE_MIN
        * maximum_right
        + left.factor
        * right.factor
        * float(columns)
        * float(inner)
        * SAFE_MIN,
        12,
    )
    return normalize_rows(RowUpperMatrix(matrix, exponents, factor, errors))


def power_rows(value: RowUpperMatrix, exponent: int) -> RowUpperMatrix:
    """Raise a square row-scaled matrix to a nonnegative integer power."""
    if exponent < 0 or value.matrix.shape[0] != value.matrix.shape[1]:
        raise ValueError("invalid matrix power")
    size = value.matrix.shape[0]
    result = RowUpperMatrix(
        np.eye(size, dtype=np.float64), np.zeros(size, dtype=np.int64)
    )
    power = value
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply_rows(result, power)
        remaining >>= 1
        if remaining:
            power = multiply_rows(power, power)
    return result


def row_sum_rows(
    value: RowUpperMatrix, row: int
) -> tuple[float, int, float, float]:
    """Return components bounding one row sum."""
    assert value.errors is not None
    columns = value.matrix.shape[1]
    mantissa = float(np.sum(value.matrix[row, :]))
    factor = product_upper(value.factor, inflation(columns))
    error = sum_upper(
        float(value.errors[row]),
        product_upper(value.factor, float(columns), SAFE_MIN),
    )
    return mantissa, int(value.exponents[row]), factor, error


__all__ = [
    "RowUpperMatrix",
    "multiply_rows",
    "normalize_rows",
    "power_rows",
    "rational_upper",
    "row_sum_rows",
    "scale_rows",
    "sum_rows",
]

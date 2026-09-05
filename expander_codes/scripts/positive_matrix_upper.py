#!/usr/bin/env python3
"""Directed upper arithmetic for scaled nonnegative binary64 matrices.

Every ``ScaledUpperMatrix`` represents an entrywise upper bound
``matrix * 2**exponent``.  Matrix products use the standard dot-product error
bound and add an absolute allowance for gradual underflow.  All other
operations round their nonnegative result toward positive infinity.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


UNIT_ROUNDOFF = 2.0**-53
MIN_SUBNORMAL = float.fromhex("0x0.0000000000001p-1022")


@dataclass(frozen=True)
class ScaledUpperMatrix:
    matrix: np.ndarray
    exponent: int


def upward(value: np.ndarray) -> np.ndarray:
    """Move every positive finite entry to the next binary64 value."""
    result = value.copy()
    np.nextafter(result, math.inf, out=result, where=result > 0.0)
    return result


def rational_upper(numerator: int, denominator: int) -> float:
    """Return a binary64 upper bound on a nonnegative rational number."""
    if numerator < 0 or denominator <= 0:
        raise ValueError("rational must be nonnegative with positive denominator")
    if numerator == 0:
        return 0.0
    nearest = numerator / denominator
    if nearest == 0.0 and numerator:
        return MIN_SUBNORMAL
    return math.nextafter(nearest, math.inf)


def scalar_multiply_upper(matrix: np.ndarray, scalar: float) -> np.ndarray:
    """Upper-bound a nonnegative matrix times a nonnegative binary64 scalar."""
    if scalar < 0.0 or not math.isfinite(scalar):
        raise ValueError("scalar must be finite and nonnegative")
    product = matrix * scalar
    inflation = math.nextafter(1.0 + 4.0 * UNIT_ROUNDOFF, math.inf)
    allowance = np.where(matrix > 0.0, MIN_SUBNORMAL, 0.0)
    return upward(product * inflation + allowance)


def add_upper(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Upper-bound the entrywise sum of two nonnegative matrices."""
    return upward(left + right)


def matmul_upper(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Upper-bound a nonnegative binary64 matrix product."""
    if left.shape[1] != right.shape[0]:
        raise ValueError("matrix dimensions do not match")
    inner = left.shape[1]
    if inner * UNIT_ROUNDOFF >= 0.25:
        raise ValueError("inner dimension is too large for the rounding bound")
    inflation = math.nextafter(
        1.0 + 4.0 * inner * UNIT_ROUNDOFF, math.inf
    )
    product = left @ right
    # At most ``inner`` products and additions can underflow in each dot
    # product.  Impossible transitions remain exactly zero; otherwise a
    # dense allowance would create spurious paths through the state machine.
    reachable = (
        left.astype(np.bool_).astype(np.uint8)
        @ right.astype(np.bool_).astype(np.uint8)
    ) > 0
    allowance = np.where(reachable, inner * MIN_SUBNORMAL, 0.0)
    return upward(product * inflation + allowance)


def normalize_upper(value: ScaledUpperMatrix) -> ScaledUpperMatrix:
    """Normalize a scaled upper matrix using an exact power-of-two factor."""
    maximum = float(np.max(value.matrix))
    if maximum == 0.0:
        return value
    _, shift = math.frexp(maximum)
    scaled = np.ldexp(value.matrix, -shift)
    scaled = np.where(
        (value.matrix > 0.0) & (scaled == 0.0), MIN_SUBNORMAL, scaled
    )
    return ScaledUpperMatrix(scaled, value.exponent + shift)


def align_upper(value: ScaledUpperMatrix, exponent: int) -> np.ndarray:
    """Express ``value`` at a greater common power-of-two exponent."""
    if exponent < value.exponent:
        raise ValueError("common exponent must not be smaller")
    shifted = np.ldexp(value.matrix, value.exponent - exponent)
    # Preserve a positive upper bound when the exact power-of-two shift falls
    # below the subnormal range.
    shifted = np.where(
        (value.matrix > 0.0) & (shifted == 0.0), MIN_SUBNORMAL, shifted
    )
    return shifted


def add_scaled_upper(
    left: ScaledUpperMatrix | None, right: ScaledUpperMatrix
) -> ScaledUpperMatrix:
    """Add two scaled upper matrices and normalize the result."""
    if left is None:
        return normalize_upper(right)
    exponent = max(left.exponent, right.exponent)
    matrix = add_upper(
        align_upper(left, exponent), align_upper(right, exponent)
    )
    return normalize_upper(ScaledUpperMatrix(matrix, exponent))


def multiply_scaled_upper(
    left: ScaledUpperMatrix, right: ScaledUpperMatrix
) -> ScaledUpperMatrix:
    """Multiply two scaled upper matrices and normalize the result."""
    return normalize_upper(ScaledUpperMatrix(
        matmul_upper(left.matrix, right.matrix),
        left.exponent + right.exponent,
    ))


def scale_scaled_upper(
    value: ScaledUpperMatrix, scalar: float
) -> ScaledUpperMatrix:
    """Multiply a scaled upper matrix by an upper scalar bound."""
    return normalize_upper(ScaledUpperMatrix(
        scalar_multiply_upper(value.matrix, scalar), value.exponent
    ))


def power_scaled_upper(value: ScaledUpperMatrix, exponent: int) -> ScaledUpperMatrix:
    """Raise a square scaled upper matrix to a nonnegative integer power."""
    if exponent < 0 or value.matrix.shape[0] != value.matrix.shape[1]:
        raise ValueError("invalid matrix power")
    size = value.matrix.shape[0]
    result = ScaledUpperMatrix(np.eye(size, dtype=np.float64), 0)
    power = value
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply_scaled_upper(result, power)
        remaining >>= 1
        if remaining:
            power = multiply_scaled_upper(power, power)
    return result


def row_sum_upper(value: ScaledUpperMatrix, row: int) -> tuple[float, int]:
    """Return an upper mantissa and exponent for one nonnegative row sum."""
    count = value.matrix.shape[1]
    if count * UNIT_ROUNDOFF >= 0.25:
        raise ValueError("row is too long for the rounding bound")
    inflation = math.nextafter(
        1.0 + 4.0 * count * UNIT_ROUNDOFF, math.inf
    )
    nearest = float(np.sum(value.matrix[row, :]))
    nonzero = int(np.count_nonzero(value.matrix[row, :]))
    upper = math.nextafter(
        nearest * inflation + nonzero * MIN_SUBNORMAL,
        math.inf,
    )
    if not math.isfinite(upper):
        raise OverflowError("nonfinite row sum")
    return upper, value.exponent

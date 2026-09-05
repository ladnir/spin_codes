#!/usr/bin/env python3
"""Rigorous layered upper bounds for nonnegative binary64 matrices.

Each logical matrix is a sum of dense binary64 layers.  A layer stores
``matrix * 2**exponent``.  Canonicalization assigns every nonzero entry to a
symmetric power-of-two band, so one logical entry can have a scale independent
of every other entry.  Nearby entries still share dense BLAS operations.

The default band width is 1,000 bits.  Stored mantissas lie roughly between
``2**-500`` and ``2**500``.  Their products remain normal and finite for the
matrix dimensions used here, so multiplication needs a relative rounding
inflation but no absolute underflow allowance.  Structural zeros remain exact.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from positive_matrix_upper import ScaledUpperMatrix, UNIT_ROUNDOFF


DEFAULT_BAND_BITS = 1000


@dataclass(frozen=True)
class LayeredUpperMatrix:
    """A nonnegative matrix upper bound represented as a sum of layers."""

    layers: tuple[ScaledUpperMatrix, ...]
    shape: tuple[int, int]
    band_bits: int = DEFAULT_BAND_BITS
    canonical: bool = True

    def __post_init__(self) -> None:
        if self.band_bits < 1 or self.band_bits >= 1022:
            raise ValueError("band width does not keep scalar products normal")
        for layer in self.layers:
            if layer.matrix.shape != self.shape:
                raise ValueError("layer shape mismatch")


def _inflation(count: int) -> float:
    if count < 1 or count * UNIT_ROUNDOFF >= 0.25:
        raise ValueError("invalid operation count")
    return math.nextafter(1.0 + 4.0 * count * UNIT_ROUNDOFF, math.inf)


def _upward_positive(values: np.ndarray) -> np.ndarray:
    result = values.copy()
    np.nextafter(result, math.inf, out=result, where=result > 0.0)
    return result


def _matmul_arrays_upper(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Multiply canonical layer mantissas without an underflow residual."""
    if left.shape[1] != right.shape[0]:
        raise ValueError("matrix dimensions do not match")
    inner = left.shape[1]
    product = left @ right
    product *= _inflation(inner)
    np.nextafter(product, math.inf, out=product, where=product > 0.0)
    # With canonical inputs, every nonzero scalar product is normal.  Hence
    # a zero output is a structural zero, not an arithmetic underflow.  Do
    # not recompute boolean reachability here: that second dense product
    # would double the cost of the hot kernel.
    return product


def _canonicalize(
    values: list[ScaledUpperMatrix],
    shape: tuple[int, int],
    band_bits: int,
) -> LayeredUpperMatrix:
    """Split entries into exponent bands and merge equal-band layers."""
    pending = [
        value for value in values if np.any(value.matrix > 0.0)
    ]
    if not pending:
        return LayeredUpperMatrix((), shape, band_bits)
    half_band = band_bits // 2

    while True:
        bucket_terms: dict[int, list[np.ndarray]] = {}
        for value in pending:
            if value.matrix.shape != shape:
                raise ValueError("matrix shape mismatch")
            positive = value.matrix > 0.0
            if not np.any(positive):
                continue
            _, local_exponents = np.frexp(value.matrix)
            absolute_exponents = local_exponents.astype(np.int64) + value.exponent
            layer_exponents = (
                ((absolute_exponents + half_band) // band_bits) * band_bits
            )
            selected_exponents = layer_exponents[positive]
            minimum_exponent = int(np.min(selected_exponents))
            maximum_exponent = int(np.max(selected_exponents))
            if minimum_exponent == maximum_exponent:
                exponent_int = minimum_exponent
                matrix = np.ldexp(
                    value.matrix, value.exponent - exponent_int
                )
                if np.any(positive & ((matrix == 0.0) | ~np.isfinite(matrix))):
                    raise ArithmeticError("layer canonicalization lost an entry")
                bucket_terms.setdefault(exponent_int, []).append(matrix)
                continue
            for exponent in np.unique(layer_exponents[positive]):
                exponent_int = int(exponent)
                mask = positive & (layer_exponents == exponent)
                matrix = np.zeros(shape, dtype=np.float64)
                matrix[mask] = np.ldexp(
                    value.matrix[mask], value.exponent - exponent_int
                )
                if np.any((matrix[mask] == 0.0) | ~np.isfinite(matrix[mask])):
                    raise ArithmeticError("layer canonicalization lost an entry")
                bucket_terms.setdefault(exponent_int, []).append(matrix)

        buckets: dict[int, np.ndarray] = {}
        for exponent, terms in bucket_terms.items():
            if len(terms) == 1:
                buckets[exponent] = terms[0]
                continue
            matrix = np.sum(np.stack(terms, axis=0), axis=0)
            matrix *= _inflation(len(terms))
            np.nextafter(matrix, math.inf, out=matrix, where=matrix > 0.0)
            buckets[exponent] = matrix

        # Merging contributions can carry an entry into the next band.  A
        # second pass is required only when such a carry occurs.
        stable = True
        for exponent, matrix in buckets.items():
            positive = matrix > 0.0
            if not np.any(positive):
                continue
            _, local_exponents = np.frexp(matrix[positive])
            absolute = local_exponents.astype(np.int64) + exponent
            expected = ((absolute + half_band) // band_bits) * band_bits
            if np.any(expected != exponent):
                stable = False
                break
        if stable:
            layers = tuple(
                ScaledUpperMatrix(buckets[exponent], exponent)
                for exponent in sorted(buckets, reverse=True)
            )
            return LayeredUpperMatrix(layers, shape, band_bits)
        pending = [
            ScaledUpperMatrix(matrix, exponent)
            for exponent, matrix in buckets.items()
        ]


def layered_from_matrix(
    matrix: np.ndarray, exponent: int = 0,
    band_bits: int = DEFAULT_BAND_BITS,
) -> LayeredUpperMatrix:
    """Create a canonical layered matrix from one nonnegative component."""
    if matrix.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    if np.any(matrix < 0.0) or not np.all(np.isfinite(matrix)):
        raise ValueError("matrix must be finite and nonnegative")
    return _canonicalize(
        [ScaledUpperMatrix(matrix.copy(), exponent)], matrix.shape, band_bits
    )


def layered_identity(
    size: int, band_bits: int = DEFAULT_BAND_BITS
) -> LayeredUpperMatrix:
    return layered_from_matrix(np.eye(size, dtype=np.float64), 0, band_bits)


def sum_layered(values: list[LayeredUpperMatrix]) -> LayeredUpperMatrix:
    """Return a directed upper bound on a sum of layered matrices."""
    if not values:
        raise ValueError("cannot sum an empty matrix family")
    shape = values[0].shape
    band_bits = values[0].band_bits
    if any(
        value.shape != shape or value.band_bits != band_bits
        for value in values
    ):
        raise ValueError("incompatible layered matrices")
    return _canonicalize(
        [layer for value in values for layer in value.layers],
        shape,
        band_bits,
    )


def multiply_layered(
    left: LayeredUpperMatrix, right: LayeredUpperMatrix
) -> LayeredUpperMatrix:
    """Return a directed upper bound on a matrix product."""
    if left.shape[1] != right.shape[0] or left.band_bits != right.band_bits:
        raise ValueError("incompatible layered matrices")
    if not left.canonical or not right.canonical:
        raise ValueError("matrix products require canonical exponent layers")
    headroom = math.ceil(math.log2(max(1, left.shape[1])))
    if left.band_bits + headroom >= 1023:
        raise ValueError("band width can overflow this matrix product")
    shape = (left.shape[0], right.shape[1])
    products = [
        ScaledUpperMatrix(
            _matmul_arrays_upper(a.matrix, b.matrix),
            a.exponent + b.exponent,
        )
        for a in left.layers
        for b in right.layers
    ]
    return _canonicalize(products, shape, left.band_bits)


def multiply_scale_layered_rational(
    left: LayeredUpperMatrix,
    right: LayeredUpperMatrix,
    numerator: int,
    denominator: int,
) -> LayeredUpperMatrix:
    """Fuse a matrix product with a rational scale for slice convolution.

    The returned components are intentionally not canonicalized.  The slice
    kernel immediately passes all such summands to ``sum_layered``, which
    canonicalizes the completed coefficient once.
    """
    if left.shape[1] != right.shape[0] or left.band_bits != right.band_bits:
        raise ValueError("incompatible layered matrices")
    headroom = math.ceil(math.log2(max(1, left.shape[1])))
    if left.band_bits + headroom >= 1023:
        raise ValueError("band width can overflow this matrix product")
    mantissa, scalar_exponent = rational_mantissa_exponent_upper(
        numerator, denominator
    )
    shape = (left.shape[0], right.shape[1])
    if mantissa == 0.0 or not left.layers or not right.layers:
        return LayeredUpperMatrix((), shape, left.band_bits)
    products = []
    for a in left.layers:
        for b in right.layers:
            product = _matmul_arrays_upper(a.matrix, b.matrix)
            product *= mantissa
            np.nextafter(product, math.inf, out=product, where=product > 0.0)
            products.append(ScaledUpperMatrix(
                product,
                a.exponent + b.exponent + scalar_exponent,
            ))
    return LayeredUpperMatrix(
        tuple(products), shape, left.band_bits, canonical=False
    )


def rational_mantissa_exponent_upper(
    numerator: int, denominator: int
) -> tuple[float, int]:
    """Upper-bound a nonnegative rational as ``mantissa * 2**exponent``."""
    if numerator < 0 or denominator <= 0:
        raise ValueError("rational must be nonnegative with positive denominator")
    if numerator == 0:
        return 0.0, 0
    exponent = numerator.bit_length() - denominator.bit_length() + 1

    def upper_significand(selected_exponent: int) -> int:
        if selected_exponent >= 0:
            scaled_numerator = numerator << 53
            scaled_denominator = denominator << selected_exponent
        else:
            scaled_numerator = numerator << (53 - selected_exponent)
            scaled_denominator = denominator
        return (
            scaled_numerator + scaled_denominator - 1
        ) // scaled_denominator

    significand = upper_significand(exponent)
    if significand < 1 << 52:
        exponent -= 1
        significand = upper_significand(exponent)
    elif significand > 1 << 53:
        exponent += 1
        significand = upper_significand(exponent)
    mantissa = math.ldexp(float(significand), -53)
    if not (0.5 <= mantissa <= 1.0):
        raise ArithmeticError("failed to normalize rational upper bound")
    return mantissa, exponent


def scale_layered_rational(
    value: LayeredUpperMatrix, numerator: int, denominator: int
) -> LayeredUpperMatrix:
    """Multiply by a nonnegative rational using a scaled upper scalar."""
    mantissa, exponent = rational_mantissa_exponent_upper(
        numerator, denominator
    )
    if mantissa == 0.0 or not value.layers:
        return LayeredUpperMatrix((), value.shape, value.band_bits)
    scaled = []
    for layer in value.layers:
        product = _upward_positive(layer.matrix * mantissa)
        positive = layer.matrix > 0.0
        if np.any(positive & (product == 0.0)):
            raise ArithmeticError("scaled layer unexpectedly underflowed")
        scaled.append(ScaledUpperMatrix(product, layer.exponent + exponent))
    return _canonicalize(scaled, value.shape, value.band_bits)


def power_layered(
    value: LayeredUpperMatrix, exponent: int
) -> LayeredUpperMatrix:
    """Raise a square layered matrix to a nonnegative integer power."""
    if exponent < 0 or value.shape[0] != value.shape[1]:
        raise ValueError("invalid matrix power")
    if not value.canonical:
        raise ValueError("matrix powers require canonical exponent layers")
    result = layered_identity(value.shape[0], value.band_bits)
    power = value
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply_layered(result, power)
        remaining >>= 1
        if remaining:
            power = multiply_layered(power, power)
    return result


def row_sum_layers(value: LayeredUpperMatrix, row: int) -> list[tuple[float, int]]:
    """Return scaled binary64 upper components for one logical row sum."""
    if not 0 <= row < value.shape[0]:
        raise ValueError("row index out of range")
    result: list[tuple[float, int]] = []
    columns = value.shape[1]
    for layer in value.layers:
        nearest = float(np.sum(layer.matrix[row, :]))
        if nearest == 0.0:
            continue
        upper = math.nextafter(nearest * _inflation(columns), math.inf)
        result.append((upper, layer.exponent))
    return result


__all__ = [
    "DEFAULT_BAND_BITS",
    "LayeredUpperMatrix",
    "layered_from_matrix",
    "layered_identity",
    "multiply_layered",
    "multiply_scale_layered_rational",
    "power_layered",
    "rational_mantissa_exponent_upper",
    "row_sum_layers",
    "scale_layered_rational",
    "sum_layered",
]

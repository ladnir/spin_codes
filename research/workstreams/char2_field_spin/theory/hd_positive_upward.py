"""Directed positive binary64 arithmetic with explicit power-of-two scales.

Every value represents an UPPER bound mantissa * 2**exponent. All
arithmetic uses nonnegative values; there is no subtraction, logarithm,
exponential, dot product, or unspecified library reduction. Each primitive
addition/multiplication/scaling rounds toward +infinity by nextafter.
Exact zero operands preserve structural zero; nonzero underflow becomes
the least positive subnormal. Overflow and failed IEEE checks are errors.

This module supplies arithmetic, not a theorem about its caller's model.
Inputs are already binary64 upper mantissas, normally produced by the
exact-checked Arb exporter. Converting arbitrary integers, decimals or
rationals here is not directed and is outside this interface.
"""
from dataclasses import dataclass

import numpy as np


def _array(value):
    value = np.asarray(value)
    if value.dtype != np.dtype(np.float64):
        raise ValueError('Inputs must already be binary64 upper mantissas; use a checked exporter')
    if np.any(~np.isfinite(value)) or np.any(value < 0):
        raise ValueError('Finite nonnegative binary64 mantissas required')
    return value


def _outward(raw, positive):
    with np.errstate(over='raise', invalid='raise', under='ignore'):
        result = np.where(positive, np.nextafter(raw, np.inf), 0.)
    if np.any(~np.isfinite(result)):
        raise ArithmeticError('Positive directed arithmetic overflowed')
    return result


def upper_add(left, right):
    left, right = _array(left), _array(right)
    with np.errstate(over='raise', invalid='raise', under='ignore'):
        value = np.add(left, right)
    return _outward(value, (left != 0) | (right != 0))


def upper_multiply(left, right):
    left, right = _array(left), _array(right)
    with np.errstate(over='raise', invalid='raise', under='ignore'):
        value = np.multiply(left, right)
    return _outward(value, (left != 0) & (right != 0))


def upper_ldexp(value, shift):
    value = _array(value)
    if type(shift) is not int or not -(1 << 30) < shift < (1 << 30):
        raise ValueError('Explicit bounded integer binary scaling required')
    if shift == 0:
        return value.copy()
    with np.errstate(over='raise', invalid='raise', under='ignore'):
        answer = np.ldexp(value, shift)
    return _outward(answer, value != 0)


@dataclass(frozen=True, slots=True)
class Scaled:
    mantissas: np.ndarray
    exponent: int

    def __post_init__(self):
        if type(self.exponent) is not int:
            raise ValueError('An exact integer scale exponent is required')
        object.__setattr__(self, 'mantissas', _array(self.mantissas))

    def __getitem__(self, index):
        return Scaled(self.mantissas[index], self.exponent)


def normalize(value, exponent=0):
    if type(exponent) is not int:
        raise ValueError('An exact integer scale exponent is required')
    value = _array(value)
    maximum = float(np.max(value)) if value.size else 0.
    if not maximum:
        return Scaled(value.copy(), 0)
    shift = int(np.frexp(maximum)[1])
    return Scaled(upper_ldexp(value, -shift), exponent + shift)


def from_export(export):
    return normalize(export['mantissas'], export['exponent'])


def zeros(shape):
    return Scaled(np.zeros(shape, dtype=np.float64), 0)


def identity(shape=()):
    value = np.zeros((*shape, 2, 2), dtype=np.float64)
    value[..., 0, 0] = value[..., 1, 1] = 1.
    return Scaled(value, 0)


def add(left, right):
    exponent = _common_exponent((left, right))
    a = upper_ldexp(left.mantissas, left.exponent - exponent)
    b = upper_ldexp(right.mantissas, right.exponent - exponent)
    return normalize(upper_add(a, b), exponent)


def multiply(left, right):
    return normalize(upper_multiply(left.mantissas, right.mantissas), left.exponent + right.exponent)


def concatenate(values, axis=0):
    values = tuple(values)
    if not values:
        raise ValueError('A nonempty list of scaled arrays required')
    exponent = _common_exponent(values)
    return normalize(np.concatenate([upper_ldexp(value.mantissas, value.exponent - exponent)
                                    for value in values], axis=axis), exponent)


def matmul(left, right):
    if left.mantissas.shape[-2:] != (2, 2) or right.mantissas.shape[-2:] != (2, 2):
        raise ValueError('Two-state matrices required')
    first = multiply(left[..., :, 0, None], right[..., None, 0, :])
    second = multiply(left[..., :, 1, None], right[..., None, 1, :])
    return add(first, second)


def matrix_power(matrix, exponent):
    if type(exponent) is not int or exponent < 0:
        raise ValueError('Nonnegative integer matrix exponent required')
    result = identity(matrix.mantissas.shape[:-2])
    while exponent:
        if exponent & 1:
            result = matmul(result, matrix)
        exponent >>= 1
        if exponent:
            matrix = matmul(matrix, matrix)
    return result


def sum_axis0(value):
    """Explicit directed pairwise reduction; never np.sum/einsum/dot."""
    current = value.mantissas
    if current.ndim < 1 or not len(current):
        raise ValueError('Nonempty leading reduction axis required')
    while len(current) > 1:
        pairs = len(current) // 2
        reduced = upper_add(current[:2*pairs:2], current[1:2*pairs:2])
        if len(current) % 2:
            reduced = np.concatenate((reduced, current[-1:]), axis=0)
        current = reduced
    return normalize(current[0], value.exponent)


def injective_add(destination, targets, values):
    """Positive scatter-add with checked injective leading-axis indices."""
    targets = np.asarray(targets)
    if (targets.ndim != 1 or targets.dtype.kind not in 'iu'
            or len(targets) != len(values.mantissas) or not len(targets)
            or np.min(targets) < 0 or np.max(targets) >= len(destination.mantissas)
            or len(np.unique(targets)) != len(targets)
            or destination.mantissas[targets].shape != values.mantissas.shape):
        raise ValueError('Injective in-range targets and matching values required')
    exponent = _common_exponent((destination, values))
    result = upper_ldexp(destination.mantissas, destination.exponent - exponent)
    addition = upper_ldexp(values.mantissas, values.exponent - exponent)
    result[targets] = upper_add(result[targets], addition)
    return normalize(result, exponent)


def _common_exponent(values):
    return max((value.exponent for value in values if np.any(value.mantissas)), default=0)


def ieee_checks():
    """Reject FTZ/DAZ or unexpected float format before a directed replay."""
    info = np.finfo(np.float64)
    if info.nmant != 52 or info.maxexp != 1024 or info.minexp != -1022:
        raise ArithmeticError('The replay requires IEEE binary64')
    tiny = np.nextafter(np.float64(0), np.float64(1))
    if not tiny > 0 or float(tiny).as_integer_ratio() != (1, 1 << 1074):
        raise ArithmeticError('Subnormal generation is unavailable')
    with np.errstate(under='ignore'):
        # Expected subnormals are integer bit patterns. Computing them by
        # ldexp under FTZ could flush BOTH sides and falsely accept zero.
        sum_bits = int(np.asarray(np.add(tiny, tiny), dtype=np.float64).view(np.uint64))
        half_bits = int(np.asarray(np.multiply(np.float64(2.)**-1022, np.float64(0.5)), dtype=np.float64).view(np.uint64))
        if (np.multiply(tiny, np.float64(2.)**1023) != np.float64(2.)**-51
                or sum_bits != 2 or half_bits != 0x0008000000000000):
            raise ArithmeticError('Subnormal input/output was flushed; fail closed')
    return dict(format='IEEE binary64', mantissa_bits=53, gradual_underflow=True,
                primitives='nextafter upward after each positive add/mul/ldexp',
                reductions='explicit pairwise directed additions')

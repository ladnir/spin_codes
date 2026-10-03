"""Positive dyadic upper arrays with an unbounded power-of-two exponent.

Every object denotes 2**exponent * value. All arithmetic is upward; nonzero
entries lost while aligning scales receive a minimum-subnormal allowance.
The exponent is an integer, not a floating logarithm.
"""
from __future__ import annotations

from fractions import Fraction
from math import frexp
import numpy as np
import outward_positive as op


class Scaled:
    __slots__ = ('value', 'exponent')

    def __init__(self, value, exponent=0):
        if type(exponent) is not int:
            raise ValueError('integer scale required')
        value = op._positive(value)
        largest = float(value.max(initial=0.))
        if not largest:
            self.value, self.exponent = value.copy(), 0
            return
        shift = frexp(largest)[1]
        self.value = op.ldexp_upper(value, -shift)
        self.exponent = exponent+shift


def from_fraction(value):
    value = Fraction(value)
    if value < 0: raise ValueError('nonnegative scalar required')
    if not value: return Scaled(0.)
    shift = value.numerator.bit_length()-value.denominator.bit_length()
    mantissa = value/Fraction(2)**shift
    return Scaled(op.upper_float(mantissa), shift)


def add(left, right):
    if not np.any(left.value): return Scaled(right.value, right.exponent)
    if not np.any(right.value): return Scaled(left.value, left.exponent)
    exponent = max(left.exponent, right.exponent)
    # The finite proof geometries use offsets much smaller than C-int range.
    # Clamp larger negative shifts: 2^-1075 times a normalized value is below
    # minsubnormal, so replacing it by minsubnormal is still an upper bound.
    def align(item):
        shift = max(-1075, item.exponent-exponent)
        return op.ldexp_upper(item.value, shift)
    return Scaled(op.add_upper(align(left), align(right)), exponent)


def matmul(left, right):
    return Scaled(op.matmul_upper(left.value, right.value), left.exponent+right.exponent)


def multiply(left, right):
    if not isinstance(right, Scaled): right = from_fraction(Fraction(right))
    return Scaled(op.multiply_upper(left.value, right.value), left.exponent+right.exponent)


def power(matrix, exponent):
    if type(exponent) is not int or exponent < 0:
        raise ValueError('nonnegative integer power required')
    if matrix.value.ndim != 2 or matrix.value.shape[0] != matrix.value.shape[1]:
        raise ValueError('square matrix required')
    result = Scaled(np.eye(matrix.value.shape[0]))
    while exponent:
        if exponent&1: result = matmul(result, matrix)
        exponent //= 2
        if exponent: matrix = matmul(matrix, matrix)
    return result


def integer_power(value, exponent):
    if not isinstance(value, Scaled): value = from_fraction(value)
    if value.value.size != 1 or type(exponent) is not int or exponent < 0:
        raise ValueError('scalar and nonnegative integer power required')
    result = Scaled(1.)
    while exponent:
        if exponent&1: result = multiply(result, value)
        exponent //= 2
        if exponent: value = multiply(value, value)
    return result


def fractional_power(value, alpha):
    alpha = Fraction(alpha)
    if value.value.size != 1 or not 0 < alpha <= 1:
        raise ValueError('scalar and fractional exponent in(0,1] required')
    if alpha == 1: return Scaled(value.value, value.exponent)
    exponent, remainder = divmod(value.exponent*alpha.numerator, alpha.denominator)
    factor = op.tangent_table(alpha)[2][remainder]
    return Scaled(op.multiply_upper(op.fractional_power_upper(value.value, alpha), factor), exponent)


def scalar_fraction(value):
    if value.value.size != 1: raise ValueError('scalar endpoint required')
    return Fraction(float(value.value.item()))*Fraction(2)**value.exponent


def scalar_record(value):
    if value.value.size != 1: raise ValueError('scalar endpoint required')
    return dict(mantissa_hex=float(value.value.item()).hex(), exponent=value.exponent)


def scalar_from_record(record):
    value = float.fromhex(record['mantissa_hex'])
    exponent = record['exponent']
    if type(exponent) is not int or not np.isfinite(value) or value < 0:
        raise ValueError('finite positive dyadic endpoint required')
    return Fraction(value)*Fraction(2)**exponent

"""Positive binary64 enclosures and certified fractional-power tangents.

Assumes round-to-nearest basic operations with gradual underflow. The runtime
check rejects FTZ/DAZ. Dot products may reorder positive terms; the gamma bound
allows two rounded operations per term. No libm power is trusted as an endpoint.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import isfinite, nextafter
import numpy as np

UNIT = Fraction(1, 2**53)
ETA = np.nextafter(0., 1.)


def upper_float(value):
    value = Fraction(value)
    result = float(value)
    if not isfinite(result):
        raise OverflowError('finite binary64 endpoint required')
    while Fraction.from_float(result) < value:
        result = nextafter(result, float('inf'))
    return result


def check_runtime():
    if np.finfo(float).nmant != 52 or ETA != float.fromhex('0x0.0000000000001p-1022'):
        raise RuntimeError('IEEE binary64 required')
    # Exercise both elementwise operations and the matrix-product backend.
    a = np.array([[np.finfo(float).tiny, ETA]])
    if (a[0,0]*.5 != float.fromhex('0x0.8p-1022')
            or a[0,1]*2 != 2*ETA
            or (a@np.array([[.5],[0.]])).item() != float.fromhex('0x0.8p-1022')
            or (np.array([[ETA]])@np.array([[2.]])).item() != 2*ETA):
        raise RuntimeError('gradual underflow required; FTZ/DAZ is unsupported')


def _positive(value):
    array = np.asarray(value, dtype=float)
    if not np.isfinite(array).all() or np.any(array < 0):
        raise ValueError('finite nonnegative endpoints required')
    return array


def _next_upper(value):
    result = np.nextafter(value, np.inf)
    if not np.isfinite(result).all():
        raise OverflowError('upward endpoint overflow')
    return result


def add_upper(left, right):
    left, right = _positive(left), _positive(right)
    value = left+right
    if not np.isfinite(value).all(): raise OverflowError('sum overflow')
    return np.where((left>0)|(right>0), _next_upper(value), 0.)


def multiply_upper(left, right):
    left, right = _positive(left), _positive(right)
    with np.errstate(under='ignore'):
        value = left*right
    if not np.isfinite(value).all(): raise OverflowError('product overflow')
    return np.where((left>0)&(right>0), _next_upper(value), 0.)


def ldexp_upper(value, exponent):
    value = _positive(value)
    with np.errstate(under='ignore'):
        shifted = np.ldexp(value, exponent)
    if not np.isfinite(shifted).all(): raise OverflowError('scale overflow')
    return np.where(value>0,_next_upper(shifted),0.)


@lru_cache(None)
def dot_inflation(length):
    if type(length) is not int or not 1 <= length < 2**50:
        raise ValueError('bounded positive dot length required')
    gamma = 2*length*UNIT/(1-2*length*UNIT)
    return upper_float(1/(1-gamma))


def matmul_upper(left, right):
    left, right = _positive(left), _positive(right)
    if left.shape[-1] != right.shape[-2]: raise ValueError('shape mismatch')
    length = left.shape[-1]
    with np.errstate(under='ignore'):
        value = left@right
    # Use a wide-enough integer for Boolean supports; our matrices are 10x10.
    support = (left>0).astype(np.int64)@(right>0).astype(np.int64)>0
    correction = upper_float(length*Fraction.from_float(float(ETA)))
    value = multiply_upper(add_upper(value,correction),dot_inflation(length))
    return np.where(support,value,0.)


def sum_upper(value, axis=0):
    value = _positive(value)
    length = value.shape[axis]
    if not length: raise ValueError('nonempty positive sum required')
    # The dot bound is conservative for summation, which has no multiplications.
    summed = value.sum(axis=axis)
    correction = upper_float(length*Fraction.from_float(float(ETA)))
    summed = multiply_upper(add_upper(summed,correction),dot_inflation(length))
    return np.where(np.any(value>0,axis=axis),summed,0.)


def root_power_upper(value, numerator, denominator, *, precision=80):
    """Exact dyadic u with u^denominator >= value^numerator."""
    value = Fraction(value)
    if value < 0 or numerator < 1 or denominator < 1 or precision < 1:
        raise ValueError('positive root parameters required')
    if not value: return Fraction(0)
    target = (value.numerator**numerator) << (precision*denominator)
    divisor = value.denominator**numerator
    bits = max(1,(target.bit_length()-divisor.bit_length()+denominator)//denominator+1)
    low,high = 0,1<<bits
    while high**denominator*divisor < target: high <<= 1
    while high-low>1:
        middle=(low+high)//2
        if middle**denominator*divisor >= target: high=middle
        else: low=middle
    result=Fraction(high,1<<precision)
    if result**denominator < value**numerator:
        raise ArithmeticError('root endpoint integer check failed')
    return result


@lru_cache(None)
def tangent_table(alpha, bins=256):
    alpha=Fraction(alpha)
    if not 0<alpha<1 or bins<1 or bins&(bins-1):
        raise ValueError('alpha in(0,1), power-of-two bin count required')
    p,q=alpha.numerator,alpha.denominator
    intercept,slope=[],[]
    for i in range(bins):
        center=Fraction(2*bins+2*i+1,4*bins)
        upper=root_power_upper(center,p,q)
        intercept.append(upper_float((1-alpha)*upper))
        slope.append(upper_float(alpha*upper/center))
    factors=[upper_float(root_power_upper(Fraction(2**r),1,q)) for r in range(q)]
    return np.asarray(intercept),np.asarray(slope),np.asarray(factors)


def fractional_power_upper(value, alpha, *, bins=256):
    """Bound x^alpha using concavity, exact root inequalities, and positive ops."""
    value,alpha=_positive(value),Fraction(alpha)
    if not 0<alpha<=1: raise ValueError('alpha in(0,1] required')
    if alpha==1: return value.copy()
    intercept,slope,factors=tangent_table(alpha,bins)
    mantissa,exponent=np.frexp(value)
    # These bin operations use exact powers of two and a dyadic subtraction.
    index=np.clip(((mantissa-.5)*(2*bins)).astype(np.int64),0,bins-1)
    tangent=add_upper(multiply_upper(slope[index],mantissa),intercept[index])
    scaled_exponent=exponent.astype(np.int64)*alpha.numerator
    remainder=scaled_exponent%alpha.denominator
    shift=scaled_exponent//alpha.denominator
    upper=ldexp_upper(multiply_upper(tangent,factors[remainder]),shift)
    return np.where(value>0,upper,0.)

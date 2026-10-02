"""Exact row-pair shell bounds from the BCH dual-distance premise.

Two independent uniform codewords give independent uniform two-bit symbols
on any fewer than d_dual coordinates. Their union weight therefore matches
Bin(n,3/4) moments through degree d_dual-1. A squared reproducing polynomial
gives pointwise shell caps; unlike CDF differences, these really bound shells.
"""
from fractions import Fraction as Q
from math import comb


def krawtchouk4(n, u, degree):
    if (any(type(v) is not int for v in (n,u,degree))
            or not 0 <= u <= n or not 0 <= degree <= n):
        raise ValueError('integer weight and degree within the block required')
    values = [1]
    if degree:
        values.append(3*n-4*u)
    for j in range(1,degree):
        numerator = (3*n-4*u-2*j)*values[j]-3*(n-j+1)*values[j-1]
        value,remainder = divmod(numerator,j+1)
        assert remainder == 0
        values.append(value)
    return values


def shell_caps(n, k, dual_distance):
    if (any(type(v) is not int for v in (n,k,dual_distance))
            or not 0 <= k <= n or not 1 <= dual_distance <= n+1):
        raise ValueError('valid length, dimension and dual-distance premise required')
    degree = (dual_distance-1)//2
    result = []
    for u in range(n+1):
        values = krawtchouk4(n,u,degree)
        kernel = sum((Q(v*v,3**j*comb(n,j)) for j,v in enumerate(values)),Q(0))
        cap = Q(1 << (2*k))/kernel
        result.append(cap.numerator//cap.denominator)
    return result

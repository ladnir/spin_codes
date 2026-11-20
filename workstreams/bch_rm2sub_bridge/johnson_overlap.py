"""Normalized Hahn kernels for constant-weight pair distributions.

Formula: Chailloux--Debris-Alazard, arXiv:2405.07666v2, section4.2,
equation22, divided by multiplicity. Johnson distance is half Hamming.
For any constant-weight set, the pair average of each kernel is nonnegative.
"""
from fractions import Fraction as F
from functools import lru_cache


@lru_cache(maxsize=None)
def hahn(n, w, degree, distance):
    assert 0 <= degree <= w <= n//2 and 0 <= distance <= w
    term = total = F(1)
    for j in range(min(degree, distance)):
        term *= -F((degree-j)*(n+1-degree-j)*(distance-j),
                   (w-j)*(n-w-j)*(j+1))
        total += term
    return total

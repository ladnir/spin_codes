"""Exact first two moments of additive polynomials of fixed-slice row counts.

Only one- and two-column moments are needed. A truncated bivariate EGF
avoids enumerating row matrices or a count distribution of dimension n.
"""
import math
from fractions import Fraction as F


def moments(n, w, q, degree):
    assert 0 < w < n and q >= 1 and degree >= 1
    p = F(w, n)
    both = F(w*(w-1), n*(n-1))
    outcomes = [(F(0)-p, F(0)-p, 1-2*p+both),
                (1-p, F(0)-p, p-both), (F(0)-p, 1-p, p-both), (1-p, 1-p, both)]
    keys = [(a, b) for a in range(2*degree+1) for b in range(2*degree+1)
            if a+b <= 2*degree and ((a <= degree and b <= degree) or a == 0 or b == 0)]
    index = {key:i for i, key in enumerate(keys)}
    schedule = [(i, j, index[a+c, b+d]) for i, (a, b) in enumerate(keys)
                for j, (c, d) in enumerate(keys) if (a+c, b+d) in index]
    def multiply(left, right):
        out = [F(0)]*len(keys)
        for i, j, k in schedule:
            if left[i] and right[j]:
                out[k] += left[i]*right[j]
        return out
    base = [sum(prob*x**a*y**b for x, y, prob in outcomes)/(math.factorial(a)*math.factorial(b))
            for a, b in keys]
    result = [F(0)]*len(keys)
    result[index[0, 0]] = F(1)
    remaining = q
    while remaining:
        if remaining & 1:
            result = multiply(result, base)
        remaining >>= 1
        if remaining:
            base = multiply(base, base)
    return {key:value*math.factorial(key[0])*math.factorial(key[1])
            for key, value in zip(keys, result)}


def additive_mean_variance(n, w, q, coefficients):
    degree = len(coefficients)-1
    table = moments(n, w, q, degree)
    one = sum(c*table[j, 0] for j, c in enumerate(coefficients))
    one_square = sum(a*b*table[i+j, 0] for i, a in enumerate(coefficients) for j, b in enumerate(coefficients))
    two = sum(a*b*table[i, j] for i, a in enumerate(coefficients) for j, b in enumerate(coefficients))
    mean = n*one
    variance = n*one_square+n*(n-1)*two-mean*mean
    assert variance >= 0
    return mean, variance

"""Scaled positive matrix-polynomial products for numerical IMT regions."""
import math
import numpy as np


def normalize(coefficients,logs):
    scales = coefficients.max(axis=(1,2))
    if np.any(scales <= 0) or not np.all(np.isfinite(scales)):
        raise ArithmeticError('Vanished matrix-polynomial coefficient')
    return coefficients/scales[:,None,None],logs+np.log(scales)


def product(left,right,maximum):
    a,sa = left
    b,sb = right
    limit = min(maximum,len(a)+len(b)-2)
    out = np.empty((limit+1,a.shape[1],a.shape[2]))
    logs = np.empty(limit+1)
    for degree in range(limit+1):
        lo,hi = max(0,degree-len(b)+1),min(degree,len(a)-1)
        indices = np.arange(lo,hi+1)
        exponents = sa[indices]+sb[degree-indices]
        pivot = max(exponents)
        weighted = a[indices]*np.exp(exponents-pivot)[:,None,None]
        out[degree] = (weighted @ b[degree-indices]).sum(axis=0)
        logs[degree] = pivot
    return normalize(out,logs)


def regions(epoch,t,length,maximum):
    assert length >= t and length % t == 0 and 0 <= maximum <= length
    assert len(epoch) == min(t,maximum)+1
    logs = epoch.max(axis=(1,2))
    if not np.all(np.isfinite(logs)):
        raise ArithmeticError('Invalid epoch coefficients')
    power = (np.exp(epoch-logs[:,None,None]),
             logs+np.array([math.log(math.comb(t,j)) for j in range(len(epoch))]))
    current = np.eye(epoch.shape[1])[None],np.array([0.])
    remaining = length//t
    while remaining:
        if remaining & 1:
            current = product(current,power,maximum)
        remaining >>= 1
        if remaining:
            power = product(power,power,maximum)
    coefficients,logs = current
    with np.errstate(divide='ignore'):
        result = np.log(coefficients)+logs[:,None,None]
    return result-np.array([math.log(math.comb(length,j)) for j in range(maximum+1)])[:,None,None]

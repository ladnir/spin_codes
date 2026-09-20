"""Scaled positive polynomial products with a log fallback for rare tails.

Direct positive convolutions avoid FFT cancellation. Per-entry scales and
a common degree tilt keep ordinary coefficients representable. Tiny
outputs are recomputed in log space rather than discarded.
"""
import math

import numpy as np


def polynomial_product(left, right, maximum):
    limit = min(maximum, len(left)+len(right)-2)
    degrees = len(left)+len(right)-2
    slope = ((float(np.max(left[-1]))-float(np.max(left[0])))
             +(float(np.max(right[-1]))-float(np.max(right[0]))))/max(1, degrees)
    eta = -slope
    a = left+eta*np.arange(len(left))[:, None, None]
    b = right+eta*np.arange(len(right))[:, None, None]
    a_scale, b_scale = np.max(a, axis=0), np.max(b, axis=0)
    with np.errstate(invalid='ignore'):
        aa = np.exp(a-np.where(np.isfinite(a_scale), a_scale, 0))
        bb = np.exp(b-np.where(np.isfinite(b_scale), b_scale, 0))
    result = np.full((limit+1, 4, 4), -np.inf)
    correction = eta*np.arange(limit+1)
    for i in range(4):
        for k in range(4):
            if not math.isfinite(a_scale[i, k]):
                continue
            for j in range(4):
                if not math.isfinite(b_scale[k, j]):
                    continue
                raw = np.convolve(aa[:, i, k], bb[:, k, j])[:limit+1]
                with np.errstate(divide='ignore'):
                    logs = np.log(raw)+a_scale[i, k]+b_scale[k, j]-correction
                # Scaling error cannot hide an underflowed summand of
                # material relative size in the ordinary branch. Use the
                # original logarithms for every tiny or zero coefficient.
                for degree in np.flatnonzero(raw < 1e-130):
                    lo, hi = max(0, degree-len(right)+1), min(degree, len(left)-1)
                    logs[degree] = np.logaddexp.reduce(left[lo:hi+1, i, k]+right[degree-hi:degree-lo+1, k, j][::-1])
                np.logaddexp(result[:, i, j], logs, out=result[:, i, j])
    return result


def region_logs(epoch, t, length, maximum):
    if length < t or length % t or maximum > length or len(epoch) != min(t, maximum)+1:
        raise ValueError('incomplete or non-native region')
    power = epoch+np.array([math.log(math.comb(t, j)) for j in range(len(epoch))])[:, None, None]
    current = None; exponent = length//t
    while exponent:
        if exponent & 1:
            current = power.copy() if current is None else polynomial_product(current, power, maximum)
        exponent >>= 1
        if exponent:
            power = polynomial_product(power, power, maximum)
    return current-np.array([math.log(math.comb(length, j)) for j in range(maximum+1)])[:, None, None]

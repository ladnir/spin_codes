"""Outward positive placement and full continuous-state occupancy bounds."""
from __future__ import annotations

from fractions import Fraction
from math import comb, log2
import numpy as np
import scaled_positive as sp


def placement(local, limit, *, epochs=8, windows=32):
    local = np.asarray(local, dtype=float)
    if (local.ndim != 3 or local.shape[0] != windows+1
            or local.shape[1] != local.shape[2]
            or type(limit) is not int or not 0 <= limit <= epochs*windows):
        raise ValueError('consistent placement geometry required')
    size = local.shape[1]
    matrices = [sp.Scaled(m) for m in local]
    previous = [sp.Scaled(np.eye(size))]
    for epoch in range(1, epochs+1):
        current = []
        for occupied in range(min(limit, epoch*windows)+1):
            denominator = comb(epoch*windows, occupied)
            terms = []
            for k in range(max(0, occupied-(epoch-1)*windows), min(windows, occupied)+1):
                coefficient = Fraction(comb(windows,k)*comb((epoch-1)*windows,occupied-k), denominator)
                terms.append(sp.multiply(sp.matmul(previous[occupied-k], matrices[k]), coefficient))
            total = terms[0]
            for term in terms[1:]: total = sp.add(total, term)
            current.append(total)
        previous = current
    return previous


def continuous_moment(regional, regions=64):
    full = sp.power(regional, regions)
    row = sp.Scaled(full.value[0:1], full.exponent)
    return sp.matmul(row, sp.Scaled(np.ones((full.value.shape[0],1))))


def occupancy_upper(regional, q, z, alpha, *, groups=256, regions=64,
                    cutoff=13107, beta=None, a=None):
    z, alpha = Fraction(z), Fraction(alpha)
    if not 0 < z < 1 or not 0 < alpha <= 1 or not 1 <= q <= groups:
        raise ValueError('valid occupancy and witness parameters required')
    beta = Fraction(2**512, (2**32-1)**8) if beta is None else Fraction(beta)
    a = ((1+z)/2)**8 if a is None else Fraction(a)
    # This exactly mirrors [beta^q z^-cutoff a^(regions*q)]^alpha.
    prefactor = sp.multiply(sp.integer_power(beta,q), sp.integer_power(1/z,cutoff))
    prefactor = sp.multiply(prefactor, sp.integer_power(a,regions*q))
    prefactor = sp.fractional_power(prefactor,alpha)
    result = sp.multiply(continuous_moment(regional,regions), prefactor)
    return sp.multiply(result, comb(groups,q))


def display_margin(value):
    """Diagnostic only: acceptance must use exact endpoint comparisons."""
    if not value.value.item(): return float('inf')
    return -log2(float(value.value.item()))-value.exponent

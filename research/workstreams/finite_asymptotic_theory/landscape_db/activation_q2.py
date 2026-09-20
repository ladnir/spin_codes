"""Activation-aware Q2 transfers and exact pair-support averaging.

Q2 requires distinct nonzero B columns, verified by the map loader. Regions
with two bits use uniform two-subsets, not independent draws with replacement.
All arithmetic is binary64; these screens are not outward certificates.
"""
import ctypes
import math
from pathlib import Path

import numpy as np
import activation_q1 as q1

HERE = Path(__file__).resolve().parent


def region_logs(t, s, spectrum, lambdas, length):
    if length % t or length < t:
        raise ValueError('complete epochs required')
    lam = np.asarray(lambdas, dtype=float)
    zero, one = q1.epoch_logs(t, s, spectrum, lam)
    den = (1 << s) - 1
    m2 = np.full(lam.shape, -np.inf)
    for w, count in spectrum.items():
        for overlap in range(3):
            if overlap <= w and 2-overlap <= t-w:
                mass = count * math.comb(w, overlap) * math.comb(t-w, 2-overlap)
                term = math.log(mass) - math.log(den) - math.log(math.comb(t, 2)) - lam*(w+2-2*overlap)
                np.logaddexp(m2, term, out=m2)
    d = min(spectrum)
    if d < 2:
        raise ValueError('Q2 arbitrary-live formula requires A distance at least two')
    two = np.full_like(zero, -np.inf)
    two[:, 0, 1] = -2*lam
    two[:, 1, 0] = -lam*(d-2)-math.log(den)
    two[:, 1, 2] = -lam*(d-2)
    two[:, 2, 0] = math.log1p(1/(den-1))+m2-math.log(den)
    two[:, 2, 2] = math.log1p(1/(den-1))+m2
    base = (zero, one+math.log(t), two+math.log(math.comb(t,2)))
    identity = np.full_like(zero, -np.inf)
    identity[:,0,0] = identity[:,1,1] = identity[:,2,2] = 0
    result = (identity, np.full_like(zero,-np.inf), np.full_like(zero,-np.inf))
    def multiply(a, b):
        return (q1.matrix_product(a[0],b[0]),
                np.logaddexp(q1.matrix_product(a[0],b[1]),q1.matrix_product(a[1],b[0])),
                np.logaddexp(np.logaddexp(q1.matrix_product(a[0],b[2]),q1.matrix_product(a[1],b[1])),q1.matrix_product(a[2],b[0])))
    n = length//t
    while n:
        if n & 1:
            result = multiply(result, base)
        n >>= 1
        if n:
            base = multiply(base, base)
    return np.stack((result[0], result[1]-math.log(length), result[2]-math.log(math.comb(length,2))), axis=1)


class PairKernel:
    def __init__(self, path=None):
        self.library = ctypes.CDLL(str(path or HERE/'native_build'/'activation_q2_kernel.dll'))
        self.function = self.library.pair_log_coefficients
        pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        self.function.argtypes = [pointer, ctypes.c_int, pointer]
        self.function.restype = ctypes.c_int

    def coefficients(self, regions, length):
        regions = np.ascontiguousarray(regions, dtype=np.float64)
        if regions.shape != (3,3,3) or np.isnan(regions).any() or np.isposinf(regions).any():
            raise ValueError('invalid region log matrices')
        output = np.empty((length+1,length+1), dtype=np.float64)
        status = self.function(regions, length, output)
        if status:
            raise ArithmeticError(f'native pair recurrence failed: {status}')
        normalizers = np.array([math.log(math.comb(length,w)) for w in range(length+1)])
        output -= normalizers[:,None]
        output -= normalizers[None,:]
        return output

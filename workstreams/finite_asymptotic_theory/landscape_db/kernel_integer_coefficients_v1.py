"""Exact truncated coefficients of an even kernel polynomial raised to E."""
import math

import numpy as np


def coefficients(kernel, epochs, maximum):
    if epochs < 1 or maximum < 0 or kernel[0] != 1 or any(kernel[1::2]):
        raise ValueError('positive epoch count and an even kernel with constant one required')
    degree = min(maximum//2, (len(kernel)-1)*epochs//2)
    terms = [(j, int(a)) for j, a in enumerate(kernel[::2]) if j and a]
    result = [1]+[0]*degree
    # F G' = E F' G for G=F^E gives a short exact recurrence. Integer
    # arithmetic retains cancellation exactly; no polynomial FFT is used.
    for n in range(1, degree+1):
        total = 0
        for j, a in terms:
            if j > n:
                break
            total += ((epochs+1)*j-n)*a*result[n-j]
        value, remainder = divmod(total, n)
        if remainder or value < 0:
            raise ArithmeticError('kernel coefficient recurrence lost exactness')
        result[n] = value
    return result


def regional_logs(kernel, epochs, maximum):
    values = coefficients(kernel, epochs, maximum)
    length = (len(kernel)-1)*epochs
    return np.array([math.log(value)-math.log(math.comb(length, 2*j)) if value else -math.inf
                     for j, value in enumerate(values)])

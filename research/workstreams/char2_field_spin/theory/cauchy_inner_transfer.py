"""Exact local weight/return enumerator for a superregular packet inner.

One fixed Cauchy matrix T maps (X,R) to (Y,Q), of dimensions t+s.
Before this step R is zero or uniform nonzero in F_q^s. On a prescribed
weight-w support, X has independent uniform nonzero amplitudes. Fresh
per-step GL state refresh and independent input diagonals realize this law
for each fixed source word. No whole-code rank/existence claim is made.

Counts use integers. Floating evaluation is explicitly only a screen.
"""

from fractions import Fraction
from functools import lru_cache
import math


def validate(width, state, q):
    if (not all(isinstance(x, int) for x in (width, state, q))
            or width < 1 or state < 1 or q & (q-1)
            or q < 2*(width+state)):
        raise ValueError('require positive dimensions and a large enough binary field')


@lru_cache(maxsize=None)
def forced_zeros(weight, state, q, equations, incoming):
    """Count fixed-support inputs with prescribed output coordinates zero."""
    if incoming not in (0, 1) or not 0 <= weight or equations < 0:
        raise ValueError('invalid weight, zero count, or incoming state flag')
    result = 0
    for removed in range(weight+1):
        free = weight-removed
        count = q**max(0, free-equations)
        if incoming:
            count = q**max(0, free+state-equations)-count
        term = math.comb(weight, removed)*count
        result += term if removed % 2 == 0 else -term
    if result < 0:
        raise ArithmeticError('negative superregular solution count')
    return result


@lru_cache(maxsize=None)
def count_table(width=16, state=4, q=256):
    """counts[w][incoming][next_nonzero][output_weight], fixed input support."""
    validate(width, state, q)
    table = []
    for weight in range(width+1):
        rows = []
        for incoming in (0, 1):
            polynomials = [[0]*(width+1) for _ in range(2)]
            for output in range(width+1):
                exact = []
                for next_zero in (False, True):
                    total = 0
                    for extra in range(output+1):
                        equations = width-output+extra+(state if next_zero else 0)
                        term = math.comb(output, extra)*forced_zeros(weight, state, q, equations, incoming)
                        total += term if extra % 2 == 0 else -term
                    exact.append(math.comb(width, output)*total)
                polynomials[0][output] = exact[1]
                polynomials[1][output] = exact[0]-exact[1]
            if any(value < 0 for row in polynomials for value in row):
                raise ArithmeticError('negative joint weight coefficient')
            expected = (q-1)**weight*(q**state-1 if incoming else 1)
            if sum(map(sum, polynomials)) != expected:
                raise ArithmeticError('local mass mismatch')
            rows.append(tuple(map(tuple, polynomials)))
        table.append(tuple(rows))
    return tuple(table)


def transfer(z, width=16, state=4, q=256):
    """Exact Fraction matrices for rational z in [0,1]."""
    z = Fraction(z)
    if not 0 <= z <= 1:
        raise ValueError('require 0 <= z <= 1')
    powers = [z**h for h in range(width+1)]
    return tuple(tuple(tuple(sum((value*powers[h] for h, value in enumerate(poly)), Fraction(0))
                             /((q-1)**w*(q**state-1 if incoming else 1))
                             for poly in row) for incoming, row in enumerate(matrix))
                 for w, matrix in enumerate(count_table(width, state, q)))


@lru_cache(maxsize=None)
def normalized_coefficients(width=16, state=4, q=256):
    import numpy as np
    return np.array([[[[value/((q-1)**w*(q**state-1 if incoming else 1)) for value in poly]
                       for poly in row] for incoming, row in enumerate(matrix)]
                     for w, matrix in enumerate(count_table(width, state, q))])


def floating_transfer(z, width=16, state=4, q=256):
    import numpy as np
    if not math.isfinite(z) or not 0 <= z <= 1:
        raise ValueError('require finite 0 <= z <= 1')
    return normalized_coefficients(width, state, q) @ np.array([z**h for h in range(width+1)])

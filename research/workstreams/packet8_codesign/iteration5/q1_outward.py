"""Outward q=1 contribution for the wider RS16/8 byte-packet construction.

Input matrices are certified ACTIVE-label upper operators in normalized
comparison measures, not thinned or fractionally powered matrices. The
initial coordinate is zero; every terminal coordinate has mass one.
All endpoint Fractions are exact dyadics. This module does not authenticate
the local producer by itself and does not claim a whole-code certificate.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import hashlib
from math import comb
from pathlib import Path
import sys

import numpy as np
import outward_positive as op
import scaled_positive as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]/'k16_design'))
import rs_outer


def dyadic_upper(value, precision=80):
    """Round a nonnegative rational upward to precision significant bits."""
    value = Fraction(value)
    if value < 0 or type(precision) is not int or precision < 2:
        raise ValueError('nonnegative rational and precision >=2 required')
    if not value:
        return Fraction(0)
    n, d = value.numerator, value.denominator
    exponent = n.bit_length()-d.bit_length()
    if exponent >= 0:
        exponent -= int(n < d << exponent)
    else:
        exponent -= int(n << -exponent < d)
    shift = exponent-precision+1
    if shift >= 0:
        denominator = d << shift
        integer = (n+denominator-1)//denominator
        result = Fraction(integer << shift)
    else:
        numerator = n << -shift
        integer = (numerator+d-1)//d
        result = Fraction(integer, 1 << -shift)
    if result < value:
        raise ArithmeticError('exact dyadic ceiling check failed')
    return result


def endpoint_record(value):
    """Small JSON-compatible exact dyadic record, even for tiny values."""
    value = Fraction(value)
    denominator = value.denominator
    if value < 0 or denominator & (denominator-1):
        raise ValueError('nonnegative dyadic endpoint required')
    return dict(numerator_hex=hex(value.numerator), exponent=1-denominator.bit_length())


def _validate(local, z, regions, epochs, cutoff):
    local = np.asarray(local, dtype=float)
    z = Fraction(z)
    if (local.ndim != 3 or local.shape[0] < 2 or local.shape[1] != local.shape[2]
            or local.shape[1] < 1 or not np.isfinite(local).all() or np.any(local < 0)):
        raise ValueError('finite nonnegative active square operator family required')
    if not 0 < z < 1:
        raise ValueError('exact weight z must lie strictly between0 and1')
    if (type(regions) is not int or regions < 1 or type(epochs) is not int or epochs < 1
            or type(cutoff) is not int or cutoff < 0):
        raise ValueError('positive geometry and nonnegative integer cutoff required')
    return local, z


def regional_operators(active, epochs):
    """Return R0 and R1; R1 averages the single active step over epochs.

    The active one-packet operator already averages its physical packet slot.
    Therefore the number of packet slots cancels from the regional average.
    """
    size = active.shape[1]
    t0, t1 = sp.Scaled(active[0]), sp.Scaled(active[1])
    zero = sp.Scaled(np.zeros((size, size)))
    r0, unnormalized_r1 = sp.Scaled(np.eye(size)), zero
    for _ in range(epochs):
        unnormalized_r1 = sp.add(sp.matmul(unnormalized_r1, t0), sp.matmul(r0, t1))
        r0 = sp.matmul(r0, t0)
    return r0, sp.multiply(unnormalized_r1, Fraction(1, epochs))


def q1_support_upper(active, z, *, regions=64, epochs=32, cutoff=13107, precision=80):
    """Return probability uppers for support sizes0..regions.

    Defaults are the wider actual24 geometry. The arithmetic accepts generic
    dimensions for independent tiny tests; callers must establish the
    normalized comparison-family contract and authenticate their local input.
    """
    active, z = _validate(active, z, regions, epochs, cutoff)
    op.check_runtime()
    r0, r1 = regional_operators(active, epochs)
    size = active.shape[1]
    initial = np.zeros((1, size))
    initial[0, 0] = 1.
    coefficients = [sp.Scaled(initial)]
    zero = sp.Scaled(np.zeros((1, size)))
    for _ in range(regions):
        following = [zero for _ in range(len(coefficients)+1)]
        for support, row in enumerate(coefficients):
            following[support] = sp.add(following[support], sp.matmul(row, r0))
            following[support+1] = sp.add(following[support+1], sp.matmul(row, r1))
        coefficients = following

    reciprocal = sp.multiply(sp.Scaled(np.ones((1, 1))), 1/z)
    cutoff_factor = sp.scalar_fraction(sp.integer_power(reciprocal, cutoff))
    terminal = sp.Scaled(np.ones((size, 1)))
    result = []
    for support, row in enumerate(coefficients):
        moment = sp.matmul(row, terminal)
        bound = sp.multiply(moment, cutoff_factor/Fraction(comb(regions, support)))
        endpoint = min(Fraction(1), sp.scalar_fraction(bound))
        result.append(dyadic_upper(endpoint, precision))
    return tuple(result)


@lru_cache(None)
def wider_shell_counts():
    """Exact expected nonzero-message shells for eight parallel GF16 RS16/8 rows."""
    counts = rs_outer.expected_group_support_counts(16, 8, 8, 4)
    if len(counts) != 65 or counts[0] != 0 or sum(counts) != (1 << 256)-1:
        raise ArithmeticError('wider outer shell dimensions or total mass mismatch')
    return counts


def min_supports(*witnesses):
    """Select any already-certified witness separately for each support."""
    if not witnesses or not witnesses[0] or any(len(w) != len(witnesses[0]) for w in witnesses):
        raise ValueError('nonempty equally sized witness vectors required')
    converted = [tuple(map(Fraction, witness)) for witness in witnesses]
    if any(x < 0 or x > 1 for witness in converted for x in witness):
        raise ValueError('probability endpoints in[0,1] required')
    return tuple(map(min, zip(*converted)))


def combine_supports(supports, *, counts=None, groups=256, precision=80):
    """Return an exact dyadic upper on the single-nonzero-group union term."""
    supports = tuple(map(Fraction, supports))
    counts = wider_shell_counts() if counts is None else tuple(map(Fraction, counts))
    if (len(supports) != len(counts) or not supports
            or any(x < 0 or x > 1 for x in supports) or any(x < 0 for x in counts)
            or type(groups) is not int or groups < 1):
        raise ValueError('matching nonnegative shells, probability endpoints, and positive groups required')
    return dyadic_upper(groups*sum(x*y for x, y in zip(supports, counts)), precision)


def q1_upper(active, z, *, precision=80):
    """Convenience interface for the fixed wider construction."""
    support = q1_support_upper(active, z, precision=precision)
    return dict(support_upper=support, bound=combine_supports(support, precision=precision))


def source_pins():
    files = [Path(__file__).resolve(), HERE/'test_q1_outward.py',
             Path(op.__file__).resolve(), Path(sp.__file__).resolve(), Path(rs_outer.__file__).resolve()]
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in files}

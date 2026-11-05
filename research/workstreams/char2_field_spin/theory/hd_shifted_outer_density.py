"""Explicit outer-measure bounds for the existing GF256 Toeplitz4 ensemble.

Four RS[64,32] constituents form one group. Every column receives an
independent uniformly sampled invertible4x4 Toeplitz matrix. Its atom cap
on any fixed nonzero input is rho=1/((q-1)*q**3), where q=256.

Fix a nonzero target with h nonzero columns. It is impossible for h<33.
Otherwise condition on the mixers at any h-32 nonzero target columns.
The shortened RS information set determines at most one source group.
The remaining32 column mixers are independent of this conditioning. Each
must map a fixed nonzero input to its target, or the event is impossible.
Thus the expected output measure is at most rho**32 at every point.
Relative to uniform GF256**256 its density is q**256*rho**32.

The inherited diagonal-Cauchy bound is weaker but remains the default.
Both bounds concern the same implemented Toeplitz4 ensemble. Changing the
bound changes neither setup nor packet operators, and the projective
factor1/255 is unchanged. No arbitrary caller-supplied density is accepted.
"""
from fractions import Fraction as F
from math import log2

from flint import arb


DEFAULT_OUTER_MEASURE = 'inherited-diagonal'
TOEPLITZ_OUTER_MEASURE = 'invertible-toeplitz'


def checked_outer_measure(mode):
    if mode not in (DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE):
        raise ValueError('An established outer bound for the GF256 Toeplitz4 ensemble is required')
    return mode


def density_exact(mode=DEFAULT_OUTER_MEASURE):
    """Exact density relative to uniform256-symbol output, not an atom cap."""
    checked_outer_measure(mode)
    if mode == DEFAULT_OUTER_MEASURE:
        return F(256**256, 255**128)
    return F(256**160, 255**32)


def log_density_float(mode=DEFAULT_OUTER_MEASURE):
    """FLOAT proposal value only, preserving the inherited expression."""
    checked_outer_measure(mode)
    if mode == DEFAULT_OUTER_MEASURE:
        return 2048-128*log2(255)
    return 1280-32*log2(255)


def log_density_arb(mode=DEFAULT_OUTER_MEASURE):
    """Directed enclosure evaluated at the caller's current Arb precision."""
    checked_outer_measure(mode)
    logarithm = arb(255).log()/arb(2).log()
    if mode == DEFAULT_OUTER_MEASURE:
        return 2048-128*logarithm
    return 1280-32*logarithm


def bound_description(mode=DEFAULT_OUTER_MEASURE):
    checked_outer_measure(mode)
    return dict(outer_measure=mode, field_order=256, constituents=4,
                constituent_length=64, constituent_dimension=32,
                density_formula=('256^256/255^128' if mode == DEFAULT_OUTER_MEASURE
                                 else '256^160/255^32'),
                setup_changed=False)

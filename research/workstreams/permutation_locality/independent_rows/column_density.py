"""One-window upper bounds for multi-packet lazy density columns.

Fix linear maps E from states to output words and B in the reverse
direction. For a fixed packet shape b_1,...,b_j, sample distinct windows
uniformly and sample a uniform weight-b_i mask within each window.
Write X for their XOR and z = exp(-lambda), with lambda >= 0.

For a fixed target state t, the full tilted kernel column sum is

    K(t) = E_X z**wt(E(t + B X) + X).

The sum includes source state zero. For one packet of weight b, let D_b
upper-bound this quantity for every t, including t = 0. Fix all packets
except one of weight b, and write Y for their XOR and R for their total
weight. The remaining packet has N-j+1 available windows. Its conditional
law is at most N/(N-j+1) times the unrestricted one-packet law pointwise.
Furthermore,

    z**wt(E(t+B Y+B x)+x+Y)
        <= z**(-R) * z**wt(E(t+B Y+B x)+x).

Here wt(Y)=R because its packet windows are disjoint. Maximizing the
one-packet average over all states absorbs the translated state t+B Y.
Each packet therefore gives K(t) <= N/(N-j+1) * z**(-R) * D_b.
Taking their minimum, and the trivial bound 1, remains valid.

For a mature source measure bounded pointwise by C, the lazy target
density is at most alpha*C*max_t K(t). A uniform source class of size h
instead contributes at most alpha*U*K(t)/h. Restricting the full column
to that source class only decreases it. These statements do not apply
the lazy factor to source state zero, whose separate transfer is kept.

This opt-in candidate changes only C->C and optionally U_i->C in the
existing unsplit eleven-coordinate envelope. shape_inner enables it only
with the column_density option. Operators must use output tilt only: no extra input-weight
penalties or density-coordinate rescalings are supported.
"""
from fractions import Fraction as Q
from numbers import Integral, Rational
from pathlib import Path
import sys

from flint import arb

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from group_rank_one_verify import up
from occupancy_memory import M, C, U
from mature_tail import L48, L56


def _rational(value):
    if not isinstance(value, (Rational, str)) or isinstance(value, bool):
        raise ValueError('an exact rational or rational string is required')
    return Q(value)


def _shape(weights, windows, packet_bits):
    if not isinstance(windows, Integral) or isinstance(windows, bool) or windows < 1:
        raise ValueError('positive integer window count required')
    if not isinstance(packet_bits, Integral) or isinstance(packet_bits, bool) or packet_bits < 1:
        raise ValueError('positive integer packet width required')
    weights = tuple(weights)
    if len(weights) > windows or any(not isinstance(b, Integral) or isinstance(b, bool)
                                     or not 1 <= b <= packet_bits for b in weights):
        raise ValueError('shape must fit distinct windows with nonzero packet weights')
    return weights


def density_caps(window_averages):
    """Read exact one-window records (mass, nonzero density, zero, denominator).

    These are the records returned by occupancy_window_average.averages().
    Both density numerators are retained; their maximum covers every target.
    """
    result = {}
    for b, record in window_averages.items():
        if len(record) != 4 or any(not isinstance(v, Integral) or isinstance(v, bool) for v in record):
            raise ValueError('four exact integer entries per window record required')
        mass, nonzero, zero, denominator = record
        if denominator <= 0 or min(mass, nonzero, zero) < 0 or max(mass, nonzero, zero) > denominator:
            raise ValueError('invalid one-window moment upper')
        result[b] = Q(max(nonzero, zero), denominator)
    return result


def _selected_caps(weights, densities):
    result = {}
    for b in set(weights):
        if b not in densities:
            raise ValueError('missing one-window density bound')
        value = _rational(densities[b])
        if not 0 <= value <= 1:
            raise ValueError('one-window density bounds must belong to [0,1]')
        result[b] = value
    return result


def exact_bound(weights, densities, z, *, windows=32, packet_bits=4):
    """Exact rational form for tests or rational output weights z in (0,1]."""
    weights = _shape(weights, windows, packet_bits)
    if len(weights) < 2:
        raise ValueError('at least two packets required')
    z = _rational(z)
    if not 0 < z <= 1:
        raise ValueError('require 0 < z <= 1')
    caps = _selected_caps(weights, densities)
    factor, total = Q(windows, windows-len(weights)+1), sum(weights)
    return min(Q(1), *(factor*z**(-(total-b))*cap for b, cap in caps.items()))


def bound(weights, window_averages, tilt, *, windows=32, packet_bits=4):
    """Outward full-column bound from exact window records and rational tilt.

    Records must upper-bound the same E, B, window geometry, and tilt.
    Their zero-target entry is required even when only nonzero target
    densities will ultimately be used.
    """
    weights = _shape(weights, windows, packet_bits)
    if len(weights) < 2:
        raise ValueError('at least two packets required')
    tilt = _rational(tilt)
    if tilt < 0:
        raise ValueError('nonnegative output tilt required')
    caps = _selected_caps(weights, density_caps(window_averages))
    lam = arb(tilt.numerator)/tilt.denominator
    factor, total = Q(windows, windows-len(weights)+1), sum(weights)
    candidates = [arb(1)]
    for b, cap in caps.items():
        scale = factor*cap
        candidates.append(up(arb(scale.numerator)/scale.denominator*(lam*(total-b)).exp()))
    return min(candidates)


def refine(t, weights, window_averages, tilt, *, spectrum=None, rounds=2):
    """Copy an unpenalized eleven-coordinate operator and tighten its density.

    Call after conversion to the indicated number of transvections. The
    lazy probability is 2**(-rounds). If spectrum is supplied, its five
    positive class sizes must match the operator's sorted U coordinates.
    All other coefficients, including transitions into zero, are unchanged.
    """
    if t.nrows() != 11 or t.ncols() != 11:
        raise ValueError('the existing eleven-coordinate envelope is required')
    if not isinstance(rounds, Integral) or isinstance(rounds, bool) or rounds < 1:
        raise ValueError('positive integer transvection count required')
    weights = _shape(weights, 32, 4)
    if any(t[source, C] != 0 for source in (M, L48, L56)):
        raise ValueError('the mature density coordinate must be unsplit')
    sizes = None
    if spectrum is not None:
        sizes = [spectrum[v] for v in sorted(spectrum)]
        if len(sizes) != 5 or any(not isinstance(h, Integral) or isinstance(h, bool) or h < 1 for h in sizes):
            raise ValueError('five positive integer uniform-class sizes required')
    result = t*1
    if len(weights) < 2:
        return result
    candidate = up(bound(weights, window_averages, tilt)*arb(2)**(-rounds))
    result[C, C] = min(t[C, C], candidate)
    if sizes is not None:
        for i, size in enumerate(sizes):
            result[U+i, C] = min(t[U+i, C], up(candidate/size))
    return result

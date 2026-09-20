"""Sparse case selection for the frozen directed all-one-split recurrence.

Skip cases already retained at -80; in replay, evaluate only cases assigned
to the current witness. Counting uses exact integer binomial recurrences.
"""
import math

import numpy as np
from flint import arb

import constant_split_grid as frozen

core, base, density, sparse = frozen.core, frozen.base, frozen.density, frozen.sparse
UNSET = 10**9


def table(lo, hi):
    return np.full((hi-lo+1, hi+1), UNSET, dtype=np.int64)


def unpack(powers, lo, hi):
    frozen.bounds(powers, lo, hi)
    result = table(lo, hi)
    for q, row in zip(range(lo, hi+1), powers):
        result[q-lo, :q+1] = row
    return result


def pack(powers, lo, hi):
    result = [powers[q-lo, :q+1].tolist() for q in range(lo, hi+1)]
    core.require(all(all(-80 <= p < UNSET for p in row) for row in result), 'Unvisited cases')
    return result


def counts(rows, d, qs):
    """Yield C(rows,q) C(q,d) 12^d for sorted selected q values."""
    factor = 12**d
    previous, choose = -1, 0
    for q0 in qs:
        q = int(q0)
        choose = choose*q//(q-d) if q == previous+1 and q > d else math.comb(q, d)
        yield math.comb(rows, q)*choose*factor
        previous = q


def evaluate(region, ps, spec, lo, hi, tilt, best, owners, witness, replay=False):
    bands = frozen.ordinary_bands()
    core.require(best.shape == owners.shape == (hi-lo+1, hi+1) and
        len(region) == hi+1 and len(ps) == len(bands) and 1 <= lo <= hi <= spec['rows'], 'Invalid selected grid')
    probs = [base.decode(p) for p in ps]
    core.require(all(0 < p < 1 for p in probs), 'Invalid probabilities')
    costs = sparse.costs_for(bands, probs, core.inputs.caps_module.caps())
    roots = [((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
    probabilities = [arb(p.numerator)/p.denominator for p in probs]
    left = frozen.up(np.array([float((r*(1-p)).upper()) for r, p in zip(roots, probabilities)]))
    right = frozen.up(np.array([float((r*p).upper()) for r, p in zip(roots, probabilities)]))
    keep = sparse.hull.indices(left, right)
    correction = (spec['cutoff']*(arb(tilt)/10).exp()).exp().upper()
    man, exp = correction.man_exp()
    correction_exp = int(exp)+int(man).bit_length()
    correction_man = float(np.nextafter(float(correction*arb(2)**(-correction_exp)), np.inf))
    mantissas, exponents = density.initial(region)
    locations = [math.comb(spec['rows'], q) for q in range(hi+1)]
    ordinary_count, visited = 1, 0
    for d, matrices, powers in frozen.folds(mantissas, exponents, left[keep], right[keep]):
        qs = np.arange(max(lo, d), hi+1, dtype=np.int64)
        hs = qs-d
        mask = owners[qs-lo, hs] == witness if replay else best[qs-lo, hs] > -80
        qs, hs = qs[mask], hs[mask]
        if len(qs):
            terminal, terminal_exp = frozen.terminal_256(matrices[hs], powers[hs])
            terminal = frozen.up(terminal*correction_man)
            previous_q, choose = -1, 0
            for q0, h0, value, exponent in zip(qs, hs, terminal, terminal_exp):
                q, h = int(q0), int(h0)
                choose = choose*q//(q-d) if q == previous_q+1 and q > d else math.comb(q, d)
                count = locations[q]*choose*ordinary_count
                power = frozen.dyadic_ceiling(value, int(exponent)+correction_exp, count)
                if replay:
                    core.require(power <= int(best[q-lo, h]), 'Replay exceeded bound at Q=%d,h=%d' % (q, h))
                elif power < best[q-lo, h]:
                    best[q-lo, h], owners[q-lo, h] = power, witness
                previous_q = q
            visited += len(qs)
        ordinary_count *= 12
    return visited

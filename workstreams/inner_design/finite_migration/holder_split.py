"""Sum band assignments using a directed entrywise l_256 recurrence."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from flint import arb
import activation_density
import split_imt as base

model, up = base.model, base.up
POLICY = 'imt-seven-state;activation-density;l256-band-sum;all-one-count-exact-v1'


def norm256(values):
    """Upper l_256 norm along axis zero, including all bands."""
    maximum = values.max(axis=0)
    assert (maximum > 0).all() and np.isfinite(maximum).all()
    powers = up(values / maximum)
    for _ in range(8):
        powers = up(powers * powers)
    total = powers[0].copy()
    for row in powers[1:]:
        total = up(total + row)
    for _ in range(8):
        total = up(np.sqrt(total))
    return up(maximum * total)


def folds(mantissas, exponents, left, right):
    current, powers = mantissas.copy(), exponents.copy()
    yield 0, current, powers
    for d in range(1, len(current)):
        common = np.maximum(powers[:-1], powers[1:])
        a = up(np.ldexp(current[:-1], (powers[:-1] - common)[:, None, None]))
        b = up(np.ldexp(current[1:], (powers[1:] - common)[:, None, None]))
        values = up(up(left[:, None, None, None] * a[None]) +
                    up(right[:, None, None, None] * b[None]))
        updated = norm256(values)
        maximum = updated.max(axis=(1, 2))
        assert (maximum > 0).all() and np.isfinite(maximum).all()
        _, shift = np.frexp(maximum)
        current = up(np.ldexp(updated, -shift[:, None, None]))
        powers = common + shift
        yield d, current, powers


def evaluate(engine, region, ps, lo, hi, tilt, best, owners, witness, replay=False):
    assert engine.n == 7 and len(region) == hi + 1
    assert best.shape == owners.shape == (hi - lo + 1, hi + 1)
    assert 1 <= lo <= hi <= engine.length and len(ps) == 12 and all(0 < p < 1 for p in ps)
    assert model.BANDS[-1] == (256,) and len(model.BANDS) == 13
    roots = [model.up((c.log() / 256).exp()) for c in engine.costs(list(ps) + [F(1)])[:-1]]
    probabilities = list(map(model.number, ps))
    left = up(np.array([float(model.up(r * (1 - p))) for r, p in zip(roots, probabilities)]))
    right = up(np.array([float(model.up(r * p)) for r, p in zip(roots, probabilities)]))
    correction = model.up((engine.cutoff * model.number(tilt).exp()).exp())
    man, exp = correction.man_exp()
    correction_exp = int(exp) + int(man).bit_length()
    correction_man = float(up(float(correction * arb(2)**(-correction_exp))))
    locations = [math.comb(engine.length, q) for q in range(hi + 1)]
    visited = 0
    for d, matrices, powers in folds(*base.sparse_ranges.initial(region, 7), left, right):
        qs = np.arange(max(lo, d), hi + 1, dtype=np.int64)
        hs = qs - d
        mask = owners[qs - lo, hs] == witness if replay else best[qs - lo, hs] > -80
        qs, hs = qs[mask], hs[mask]
        if len(qs):
            terminal, terminal_exp = base.terminal_256(matrices[hs], powers[hs])
            terminal = up(terminal * correction_man)
            previous_q, binomial = -1, 0
            for q0, h0, value, exponent in zip(qs, hs, terminal, terminal_exp):
                q, h = int(q0), int(h0)
                binomial = binomial * q // (q - d) if q == previous_q + 1 and q > d else math.comb(q, d)
                # The norm recurrence has already summed every ordinary-band
                # assignment. There is no additional 12^d factor.
                count = locations[q] * binomial
                power = base.fold.dyadic_ceiling(value, int(exponent) + correction_exp, count)
                if replay:
                    assert power <= int(best[q - lo, h]), (q, h, power, int(best[q - lo, h]))
                elif power < best[q - lo, h]:
                    best[q - lo, h], owners[q - lo, h] = power, witness
                previous_q = q
            visited += len(qs)
    return visited


def run(*args):
    original = base.ladder, base.evaluate, base.POLICY
    base.ladder = SimpleNamespace(Engine=activation_density.Engine,
                                  candidate=activation_density.ladder.candidate)
    base.evaluate, base.POLICY = evaluate, POLICY
    try:
        return base.run(*args)
    finally:
        base.ladder, base.evaluate, base.POLICY = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=(16, 18), default=16)
    p.add_argument('--first', type=int, default=218)
    p.add_argument('--last', type=int, default=218)
    p.add_argument('--witness', action='append')
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    if not a.verify and not a.witness:
        p.error('Provide witness proposals')
    jobs = []
    for text in a.witness or []:
        tilt, q, h = text.split(':')
        jobs.append(dict(tilt=str(F(tilt) / 10), anchor=[int(q), int(h)]))
    run(a.output.resolve(), a.m, a.first, a.last, jobs, a.verify)

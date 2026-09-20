"""Singleton outer weights inside the IMT l_256 band-sum bound."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
from flint import arb
import holder_split as holder
import activation_density

base, model, up = holder.base, holder.model, holder.up
POLICY = 'imt-seven-state;activation-density;l256-singleton-sum;all-one-count-exact-v1'


class Engine(activation_density.Engine):
    def __init__(self, exponent):
        super().__init__(exponent)
        self.weights = tuple(w for w in model.base.WEIGHTS if w != 256)
        assert set(self.weights) == set(model.outer.inputs.caps_module.caps()) - {256}

    def costs(self, ps):
        assert len(ps) == len(self.weights) and all(0 < p < 1 for p in ps)
        caps = model.outer.inputs.caps_module.caps()
        return [model.up(model.number(F(caps[w], math.comb(256, w))) /
                         model.number(p)**w / (1 - model.number(p))**(256 - w))
                for w, p in zip(self.weights, ps)]


def choose(engine, region, ordinary):
    assert ordinary > 0 and len(region) == ordinary + 1
    logs = np.array([[float(v.log()) if v > 0 else -math.inf for v in row]
                     for row in region]).reshape(-1, 7, 7)
    floating = np.exp(logs - float(logs.max()))
    j = np.arange(ordinary + 1)
    log_counts = gammaln(ordinary + 1) - gammaln(j + 1) - gammaln(ordinary - j + 1)
    result = []
    for weight in engine.weights:
        def objective(theta):
            p = 1 / (1 + math.exp(-theta))
            lp, lq = math.log(p), math.log1p(-p)
            probabilities = np.exp(log_counts + j * lp + (ordinary - j) * lq)
            matrix = np.einsum('i,ijk->jk', probabilities, floating)
            return base.sparse_bch.terminal_log(matrix) - ordinary * (weight * lp + (256 - weight) * lq)
        fit = minimize_scalar(objective, bounds=(-8, 12), method='bounded')
        assert fit.success and math.isfinite(fit.fun)
        result.append(F.from_float(1 / (1 + math.exp(-float(fit.x)))))
    return result


def evaluate(engine, region, ps, lo, hi, tilt, best, owners, witness, replay=False):
    assert engine.n == 7 and len(region) == hi + 1
    assert best.shape == owners.shape == (hi - lo + 1, hi + 1)
    assert 1 <= lo <= hi <= engine.length
    roots = [model.up((c.log() / 256).exp()) for c in engine.costs(ps)]
    probabilities = list(map(model.number, ps))
    left = up(np.array([float(model.up(r * (1 - p))) for r, p in zip(roots, probabilities)]))
    right = up(np.array([float(model.up(r * p)) for r, p in zip(roots, probabilities)]))
    correction = model.up((engine.cutoff * model.number(tilt).exp()).exp())
    man, exp = correction.man_exp()
    correction_exp = int(exp) + int(man).bit_length()
    correction_man = float(up(float(correction * arb(2)**(-correction_exp))))
    locations = [math.comb(engine.length, q) for q in range(hi + 1)]
    visited = 0
    for d, matrices, powers in holder.folds(*base.sparse_ranges.initial(region, 7), left, right):
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
    original = base.ladder, base.evaluate, base.choose, base.POLICY
    base.ladder = SimpleNamespace(Engine=Engine, candidate=activation_density.ladder.candidate)
    base.evaluate, base.choose, base.POLICY = evaluate, choose, POLICY
    try:
        return base.run(*args)
    finally:
        base.ladder, base.evaluate, base.choose, base.POLICY = original


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

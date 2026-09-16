"""Independent-map moments inside the existing outer-only two-tilt reduction.

The inherited scalar bound does not assume an inner spectrum or mixer. This
adapter replaces every matrix and terminal moment; no old inner transfer runs.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path

import bch_model as model
import numpy as np
from flint import arb, ctx
import ladder_dense_fixed_reference as fixed
import poisson_density_factor as density


class Checker(fixed.Checker):
    def __init__(self, exponent, ps=None):
        self.engine = model.Engine(exponent)
        self.rows, self.cutoff = self.engine.length, self.engine.cutoff
        self.t, self.s, self.power = 128, 19, exponent-6
        assert 1 << self.power == self.engine.output_bits//128
        self.ps = tuple(fixed.labels.row_probabilities(.5) if ps is None else map(F, ps))
        self.bands = model.BANDS[:-1]
        self.caps = model.outer.inputs.caps_module.caps()
        assert len(self.ps) == len(self.bands) and all(0 < p < 1 for p in self.ps)
        assert sorted(w for b in self.bands for w in b) == sorted(w for w in self.caps if w != 256)
        self.pmin = min(self.ps)
        self.float_p = np.array([float(p) for p in self.ps]+[1.])
        self.shells = [np.array([math.log(self.caps[w])-math.log(math.comb(256, w))
                       -w*math.log(float(p))-(256-w)*math.log1p(-float(p)) for w in band])
                       for band, p in zip(self.bands, self.ps)]
        self.shell_logs = [[model.number(F(self.caps[w], math.comb(256, w)) /
                           (p**w*(1-p)**(256-w))).log() for w in band]
                           for band, p in zip(self.bands, self.ps)]
        self.epoch = lru_cache(maxsize=256)(self._epoch)
        self.fixed_moment = lru_cache(maxsize=1024)(self._fixed_moment)

    def _epoch(self, index):
        assert type(index) is int and -400 <= index <= 160
        lam = (arb(index)/40).exp()
        coefficients = self.engine.epoch(F(index, 40))
        floating = np.array([[float(v) for v in row] for row in coefficients]).reshape(-1, self.engine.n, self.engine.n)
        return lam, None, None, floating

    def _fixed_moment(self, index, r):
        assert type(index) is int and -400 <= index <= 160 and 0 < r < 1
        lam = (arb(index)/40).exp()
        coefficients = self.engine.epoch(F(index, 40))
        weights = fixed.labels.intervals.bernstein_weights(128, model.number(r))
        matrix = tuple(model.up(sum((w*row[k] for w, row in zip(weights, coefficients)), arb(0)))
                       for k in range(self.engine.n*self.engine.n))
        ordinary = model.independent.terminal(matrix, self.engine.n, 1 << self.power)
        direct = model.independent.terminal(self.engine.bernoulli(model.number(r), lam), self.engine.n, 1 << self.power)
        # Each complete moment is an upper bound; never mix representation entries.
        return model.up(self.cutoff*lam+min(model.up(ordinary), model.up(direct)).log())

    def bound(self, lo, hi, a, b, witness):
        value = super().bound(lo, hi, a, b, witness)
        replacement = 256*(arb(self.rows+1).log()-arb(density.density_factor(self.rows)).log())
        return model.up(value-replacement)


def screen(output):
    assert not output.exists()
    ctx.prec = 256
    results = []
    for exponent in (16, 18, 20):
        checker = Checker(exponent)
        for q in sorted(set([checker.rows//16, checker.rows//4, checker.rows//2, checker.rows])):
            for v in (F(0), F(1, 4), F(1, 2), F(3, 4), F(999, 1000)):
                witness = checker.witness(q, q, v, v)
                upper = checker.bound(q, q, v, v, witness)
                margin = float(-upper/arb(2).log())
                results.append(dict(message_exponent=exponent, occupation=q, coordinate=str(v),
                                    margin_bits=margin, witness=witness))
                print('dense point', exponent, q, str(v), margin, flush=True)
    model.base.write_new(output, dict(status='OUTWARD_POINT_DIAGNOSTICS_NOT_DENSE_COVER',
                                     full_distance_proved=False, points=results, source_sha256=model.sources()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=model.HERE/'DENSE_POINTS.json')
    args = parser.parse_args()
    screen(args.output.resolve())

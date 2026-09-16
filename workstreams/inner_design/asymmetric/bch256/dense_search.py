"""Retune dense witnesses using the independent-map Fourier moment itself."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp, expit
from flint import arb, ctx

import dense_bch as initial
import screen_dense as direct
model = initial.model


class Checker(initial.Checker):
    def witness(self, lo, hi, a, b):
        nlo, nhi = self.density(a, b)
        nu, q = float((nlo+nhi)/2), (lo+hi)/2
        theta = q/self.rows*nu
        assert 0 < theta < 1
        log_odds = math.log(theta)-math.log1p(-theta)
        def evaluate(z, u):
            r = float(expit(log_odds-u))
            if not 1e-13 < r < 1-1e-13 or not -10 <= z <= 4:
                return math.inf, 0., r
            costs = np.array([float(np.max(g+np.array(band)*u))
                              for g, band in zip(self.shells, self.bands)]+[256*u])
            def derivative(e):
                values = costs-e*self.float_p
                return nu-float(np.exp(values-logsumexp(values))@self.float_p)
            eta = -10000. if derivative(-10000.) >= 0 else 10000. if derivative(10000.) <= 0 else brentq(derivative, -10000., 10000.)
            scalar = q*(float(logsumexp(costs-eta*self.float_p))+eta*nu)
            scalar += 256*self.rows*(math.log1p(-theta)-math.log1p(-r))
            matrix = direct.bernoulli(self.engine.spectrum, self.engine.b_spectrum, self.engine.kernel, r, math.exp(z))
            moment = model.independent.g.terminal(matrix, 1 << self.power)
            return scalar+moment+self.cutoff*math.exp(z), eta, r
        starts = []
        prediction = math.log(q/self.rows)
        for z in np.linspace(max(-10, prediction-2), min(3, prediction+3), 9):
            for u in np.linspace(-3, 3, 9):
                value, eta, r = evaluate(float(z), float(u))
                starts.append((value, float(z), float(u), eta, r))
        best = min(starts)
        for _, z, u, _, _ in sorted(starts)[:2]:
            optimum = minimize(lambda v: evaluate(*v)[0]/self.engine.output_bits, [z, u],
                               method='Nelder-Mead', bounds=[(-10, 4), (-8, 8)],
                               options=dict(maxiter=120, xatol=1e-5, fatol=1e-10))
            z, u = map(float, optimum.x)
            value, eta, r = evaluate(z, u)
            if value < best[0]:
                best = value, z, u, eta, r
        _, z, u, _, _ = best
        # Quantization changes only the proposed witness, never an accepted bound.
        candidates = []
        for index in (math.floor(40*z), math.ceil(40*z)):
            value, eta, r = evaluate(index/40, u)
            candidates.append((value, dict(tilt=index, input_tilt=round(40*u),
                               eta=model.base.encode(F(round(4*eta), 4)),
                               fixed_r=model.base.encode(F.from_float(r)))))
        return min(candidates, key=lambda v: v[0])[1]


def run(output):
    assert not output.exists()
    ctx.prec = 256
    prior = model.base.read(model.HERE/'DENSE_POINTS.json')
    results = []
    for exponent in (16, 18, 20):
        checker = Checker(exponent)
        points = [p for p in prior['points'] if p['message_exponent'] == exponent and p['margin_bits'] < 40]
        for row in points:
            q, v = row['occupation'], F(row['coordinate'])
            witness = checker.witness(q, q, v, v)
            upper = checker.bound(q, q, v, v, witness)
            margin = float(-upper/arb(2).log())
            results.append(dict(message_exponent=exponent, occupation=q, coordinate=str(v),
                                old_margin_bits=row['margin_bits'], margin_bits=margin, witness=witness))
            print('retuned', exponent, q, str(v), margin, flush=True)
    model.base.write_new(output, dict(status='OUTWARD_RETUNED_POINTS_NOT_DENSE_COVER',
                                     full_distance_proved=False, points=results, source_sha256=model.sources()))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=model.HERE/'DENSE_RETUNED_POINTS.json')
    args = p.parse_args()
    run(args.output.resolve())

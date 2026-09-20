"""Select witnesses against the actual fixed-input bound, not an old surrogate."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit, logsumexp
from flint import arb, ctx
import short_dense

model = short_dense.mixed_dense.ladder.model


def proposal(checker, q, coordinate, index, u):
    nu = float(checker.density(coordinate, coordinate)[0])
    theta = q / checker.rows * nu
    r = float(expit(math.log(theta) - math.log1p(-theta) - u))
    costs = np.array([float(np.max(g + np.array(band)*u))
                      for g, band in zip(checker.shells, checker.bands)] + [256*u])
    def derivative(eta):
        values = costs - eta*checker.float_p
        return nu - float(np.exp(values - logsumexp(values)) @ checker.float_p)
    eta = (-10000. if derivative(-10000.) >= 0 else
           10000. if derivative(10000.) <= 0 else brentq(derivative, -10000., 10000.))
    return dict(tilt=index, input_tilt=round(40*u),
                eta=model.base.encode(F(round(64*eta), 64)),
                fixed_r=model.base.encode(F.from_float(r)))


def run(seed, output, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'IMT_DIRECT_FIXED_INPUT_POINT'
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    prior = model.base.read(seed)
    model.authenticate(prior)
    point = prior['points'][1]
    q, coordinate = point['q'], F(point['coordinate'])
    checker = short_dense.Checker(16)
    assert checker.engine.identity() == prior['instance']
    if saved:
        assert saved['q'] == q and saved['coordinate'] == str(coordinate)
        witness = saved['witness']
        upper = checker.fixed_input_bound(q, q, coordinate, coordinate, witness)
        assert upper <= model.number(model.base.decode(saved['log_upper']))
        model.base.write_new(output.with_name(output.stem+'_replay.json'), dict(
            status='IMT_DIRECT_FIXED_INPUT_POINT_512_REPLAY_PASSED',
            producer_sha256=model.base.sha(output), full_distance_proved=False))
        print('replay passed', float(-upper/arb(2).log()), flush=True)
        return
    best = None
    for index in range(-24, 1, 2):
        def evaluate(u):
            nonlocal best
            witness = proposal(checker, q, coordinate, index, float(u))
            upper = checker.fixed_input_bound(q, q, coordinate, coordinate, witness)
            score = float(upper)
            if best is None or score < best[0]:
                best = score, witness, model.base.encode(model.exact(upper))
            return score
        grid = np.linspace(-1.5, 1.5, 13)
        scores = [evaluate(u) for u in grid]
        i = int(np.argmin(scores))
        minimize_scalar(evaluate, bounds=(grid[max(0, i-1)], grid[min(12, i+1)]),
                        method='bounded', options=dict(xatol=1e-5, maxiter=35))
        print('direct tilt', index, 'best margin', -best[0]/math.log(2), flush=True)
    sources = short_dense.mixed_dense.ladder.candidate.sources()
    sources[seed.relative_to(model.ROOT).as_posix()] = model.base.sha(seed)
    model.base.write_new(output, dict(status='IMT_DIRECT_FIXED_INPUT_POINT',
        instance=checker.engine.identity(), q=q, coordinate=str(coordinate),
        witness=best[1], log_upper=best[2], margin_bits=-best[0]/math.log(2),
        full_distance_proved=False, source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.seed.resolve(), a.output.resolve(), a.verify)

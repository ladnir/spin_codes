"""Selected sparse occupancies with outward evaluation and fixed-witness replay."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb, ctx
import bch_model as model


def terminal_log(matrix):
    matrix = matrix.copy()
    scale = float(matrix.max())
    if not scale > 0:
        return math.inf
    matrix /= scale
    logarithm = math.log(scale)
    for _ in range(8):
        matrix = matrix@matrix
        scale = float(matrix.max())
        if not scale > 0:
            return math.inf
        matrix /= scale
        logarithm = 2*logarithm+math.log(scale)
    total = float(matrix[0].sum())
    return math.log(total)+logarithm if total > 0 else math.inf


def choose(engine, region, q):
    n = engine.n
    logs = np.array([[float(v.log()) if v > 0 else -math.inf for v in row] for row in region]).reshape(-1, n, n)
    offset = float(logs.max())
    floating = np.exp(logs-offset)
    counts = model.outer.inputs.caps_module.caps()
    probabilities = []
    for band in model.BANDS[:-1]:
        weights = np.array(band)
        costs = np.array([math.log(counts[w])-math.log(math.comb(256, w)) for w in band])
        def objective(theta):
            p = 1/(1+math.exp(-theta))
            density = float(np.max(costs-weights*math.log(p)-(256-weights)*math.log1p(-p)))
            binomial = np.array([math.comb(q, j)*p**j*(1-p)**(q-j) for j in range(q+1)])
            return terminal_log(np.einsum('i,ijk->jk', binomial, floating))+q*density
        fit = minimize_scalar(objective, bounds=(-8, 12), method='bounded')
        assert fit.success and math.isfinite(fit.fun)
        probabilities.append(F.from_float(1/(1+math.exp(-float(fit.x)))))
    probabilities.append(F(1))
    return probabilities


def run(output, occupations, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'OUTWARD_SELECTED_BCH256_ASYMMETRIC_OCCUPANCIES'
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    results = []
    for exponent in (16, 18, 20):
        engine = model.Engine(exponent)
        selected = [r for r in saved['results'] if r['instance']['message_bits'] == 1 << exponent] if saved else None
        for index, q in enumerate(occupations):
            old = selected[index] if selected else None
            if old:
                assert old['occupation'] == q and old['instance'] == engine.identity()
                trials = [(F(old['log_tilt']), [model.base.decode(p) for p in old['probabilities']])]
            else:
                prediction = (-76+6*math.log2(q)-10*(exponent-20)*math.log(2))/10
                trials = [(F.from_float(prediction+d), None) for d in (-.5, 0, .5)]
            best = None
            for tilt, probabilities in trials:
                region = engine.region(tilt, q)
                if probabilities is None:
                    probabilities = choose(engine, region, q)
                upper = model.up(engine.adaptive(region, probabilities, q)*(engine.cutoff*model.number(tilt).exp()).exp())
                if best is None or upper < best[0]:
                    best = upper, tilt, probabilities
            value, tilt, probabilities = best
            if old:
                assert value <= model.unpack(old['upper'])
                bound = old['upper']
            else:
                bound = model.pack(value)
            result = dict(instance=engine.identity(), occupation=q, log_tilt=str(tilt),
                          probabilities=[model.base.encode(p) for p in probabilities], upper=bound,
                          margin_bits=float(-model.unpack(bound).log()/arb(2).log()))
            if old:
                assert result == old
            results.append(result)
            print(ctx.prec, 'sparse', exponent, q, result['margin_bits'], flush=True)
        engine.region.cache_clear()
        engine.epoch.cache_clear()
    if saved:
        model.base.write_new(output.with_name(output.stem+'_replay.json'),
                             dict(status='BCH256_ASYMMETRIC_SPARSE_512_BIT_REPLAY_PASSED',
                                  producer_sha256=model.base.sha(output), full_distance_proved=False))
    else:
        model.base.write_new(output, dict(status='OUTWARD_SELECTED_BCH256_ASYMMETRIC_OCCUPANCIES',
                                         occupations=occupations, full_distance_proved=False,
                                         results=results, source_sha256=model.sources()))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=model.HERE/'SPARSE_SELECTED.json')
    p.add_argument('--occupations', type=int, nargs='+', default=[2, 3, 4, 8, 16, 32])
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    occupations = model.base.read(args.output)['occupations'] if args.verify else sorted(set(args.occupations))
    assert occupations and all(2 <= q <= 512 for q in occupations)
    run(args.output.resolve(), occupations, args.verify)

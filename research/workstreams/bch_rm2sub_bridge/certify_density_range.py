"""Outward fixed-weight ranges using the four-state activation-density transfer."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import numpy as np
from flint import arb, ctx
from scipy.optimize import minimize_scalar

import activation_density_arb as engine
import bch256_dense_types as dependencies
import frontier_sparse as sparse
import screen_constant_row_split as pure

core, base = engine.core, engine.core.base


def initial(region):
    values, exponents = [], []
    for row in region:
        m, e = max(row).man_exp()
        exponent = int(e)+int(m).bit_length() if m else 0
        values.append([float(v*arb(2)**(-exponent)) for v in row])
        exponents.append(exponent)
    return np.nextafter(np.array(values).reshape(-1, 4, 4), np.inf), np.array(exponents, dtype=np.int64)


def choose(region, q, caps):
    selected = np.array([[float(v.log()) if v > 0 else -math.inf for v in row]
                         for row in region[:q+1]]).reshape(-1, 4, 4)
    ps = []
    for band in sparse.BANDS:
        if band == (256,):
            ps.append(F(1))
            continue
        w = np.array(band)
        offsets = np.array([math.log(caps[x])-math.log(math.comb(256, x)) for x in band])
        def objective(theta):
            p = 1/(1+math.exp(-theta))
            g = float(np.max(offsets-w*math.log(p)-(256-w)*math.log1p(-p)))
            return pure.moment(selected, q, p)+q*g
        fit = minimize_scalar(objective, bounds=(-6, 14), method='bounded', options={'xatol': 1e-6})
        core.require(fit.success, 'Witness optimizer failed')
        ps.append(F.from_float(1/(1+math.exp(-float(fit.x)))))
    return [base.encode(p) for p in ps]


def evaluate(region, probabilities, lo, hi, spec, tilt, caps):
    ps = [base.decode(p) for p in probabilities]
    costs = sparse.costs_for(sparse.BANDS, ps, caps)
    roots = [((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
    probs = [arb(p.numerator)/p.denominator for p in ps]
    left = np.nextafter(np.array([float((r*(1-p)).upper()) for r, p in zip(roots, probs)]), np.inf)
    right = np.nextafter(np.array([float((r*p).upper()) for r, p in zip(roots, probs)]), np.inf)
    keep = sparse.hull.indices(left, right)
    mantissas, exponents = initial(region)
    correction = (spec['cutoff']*(arb(tilt)/10).exp()).exp()
    results = []
    for q, mat, exponent in sparse.scaled.matrices(mantissas, exponents, left[keep], right[keep]):
        if q < lo:
            continue
        value = engine.power(tuple(arb(float(v)) for v in mat.flat), 256)
        upper = sparse.rational((sum(value[:4], arb(0))*correction*arb(2)**(256*exponent)).upper())
        upper *= math.comb(spec['rows'], q)*len(sparse.BANDS)**q
        results.append(dict(occupation=q, upper_power=sparse.dyadic_power(upper),
                            diagnostic_margin_bits=math.log2(upper.denominator)-math.log2(upper.numerator)))
    core.require([r['occupation'] for r in results] == list(range(lo, hi+1)), 'Coverage gap')
    return results


def run(output, m, lo, hi, tilts, verify=False):
    saved = base.read(output) if verify else None
    if saved:
        core.authenticate(saved, base.ROOT)
        m = saved['instance']['message_exponent']
        lo, hi = saved['interval']
        tilts = [s['tilt'] for s in saved['shards']]
    else:
        core.require(not output.exists(), 'Use a fresh output')
    spec = core.instance('t128_s19', m)
    if saved:
        core.require(spec == saved['instance'], 'Wrong replay instance')
    core.require(2 <= lo <= hi <= spec['rows'] and tilts and
                 all(type(t) is int and -120 <= t <= 30 for t in tilts), 'Invalid interval/tilts')
    t, s, ac, kernel = core.inputs.load(spec['configuration'])
    caps = core.inputs.caps_module.caps()
    ctx.prec = 512 if verify else 256
    best, shards = {}, []
    for index, tilt in enumerate(tilts):
        region = engine.regions(t, s, ac, kernel, (-(arb(tilt)/10).exp()).exp(), hi, spec['rows'])
        print('four-state region ready', lo, hi, tilt, ctx.prec, flush=True)
        ps = saved['shards'][index]['p'] if saved else choose(region, lo, caps)
        rows = evaluate(region, ps, lo, hi, spec, tilt, caps)
        if saved:
            previous = saved['shards'][index]['rows']
            core.require([r['occupation'] for r in previous] == list(range(lo, hi+1)), 'Saved coverage gap')
            core.require(all(a['upper_power'] <= b['upper_power'] for a, b in zip(rows, previous)), 'Outward replay failed')
            rows = previous
        for row in rows:
            q, power = row['occupation'], row['upper_power']
            best[q] = min(best.get(q, power), power)
        shards.append(dict(tilt=tilt, p=ps, rows=rows))
        print('anchor margin', rows[0]['diagnostic_margin_bits'], 'last margin', rows[-1]['diagnostic_margin_bits'], flush=True)
    passing = core.intervals(q for q, p in best.items() if p <= -50)
    if saved:
        core.require([best[q] for q in range(lo, hi+1)] == saved['upper_powers'] and
                     passing == saved['per_occupancy_50_bit_intervals'], 'Saved aggregation mismatch')
        base.write_new(output.with_name(output.stem+'_replay.json'), dict(
            status='ACTIVATION_DENSITY_512_BIT_REPLAY_PASSED', producer_sha256=base.sha(output),
            instance=spec, interval=[lo, hi], per_occupancy_50_bit_intervals=passing,
            source_sha256=dependencies.sources()))
    else:
        base.write_new(output, dict(status='ACTIVATION_DENSITY_RANGE_OUTWARD_PRODUCER', instance=spec,
            interval=[lo, hi], shards=shards, upper_powers=[best[q] for q in range(lo, hi+1)],
            per_occupancy_50_bit_intervals=passing, source_sha256=dependencies.sources(), full_distance_proved=False))
    print('50-bit per-occupancy coverage', passing, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, default=16)
    p.add_argument('--first', type=int, default=150)
    p.add_argument('--last', type=int, default=512)
    p.add_argument('--tilts', nargs='+', type=int, default=[-5, 0, 5])
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output, a.m, a.first, a.last, a.tilts, a.verify)

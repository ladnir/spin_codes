"""Certify one complete occupancy by fixing the number of all-one rows.

For Q occupied rows, h all-one rows force h ones per region; d=Q-h
ordinary rows use the 12 nonconstant band majorants. All h=0..Q are summed.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import numpy as np
from flint import arb, ctx

import certify_density_range as density

core, base, sparse = density.core, density.base, density.sparse


def constant_matrices(mantissas, exponents, left, right):
    """The frozen directed fold, retaining the final offset h=Q-d instead of 0."""
    current, powers = mantissas.copy(), exponents.copy()
    yield 0, current[-1].copy(), int(powers[-1])
    def up(x):
        return np.nextafter(x, np.inf)
    for d in range(1, len(current)):
        common = np.maximum(powers[:-1], powers[1:])
        a = up(np.ldexp(current[:-1], (powers[:-1]-common)[:, None, None]))
        b = up(np.ldexp(current[1:], (powers[1:]-common)[:, None, None]))
        updated = up(up(left[0]*a)+up(right[0]*b))
        for x, y in zip(left[1:], right[1:]):
            np.maximum(updated, up(up(x*a)+up(y*b)), out=updated)
        maximum = updated.max(axis=(1, 2))
        core.require(bool((maximum > 0).all() and np.isfinite(maximum).all()), 'Invalid scaled recurrence')
        _, shift = np.frexp(maximum)
        current = up(np.ldexp(updated, -shift[:, None, None]))
        powers = common+shift
        yield d, current[-1].copy(), int(powers[-1])


def evaluate(region, ps, spec, q, tilt):
    bands = sparse.BANDS[:-1]
    core.require(sparse.BANDS[-1] == (256,) and len(ps) == len(bands), 'Wrong ordinary bands')
    costs = sparse.costs_for(bands, [base.decode(p) for p in ps], core.inputs.caps_module.caps())
    roots = [((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
    probabilities = [arb(p.numerator)/p.denominator for p in map(base.decode, ps)]
    left = np.nextafter(np.array([float((r*(1-p)).upper()) for r, p in zip(roots, probabilities)]), np.inf)
    right = np.nextafter(np.array([float((r*p).upper()) for r, p in zip(roots, probabilities)]), np.inf)
    keep = sparse.hull.indices(left, right)
    correction = (spec['cutoff']*(arb(tilt)/10).exp()).exp()
    mantissas, exponents = density.initial(region)
    results = []
    for d, mat, exponent in constant_matrices(mantissas, exponents, left[keep], right[keep]):
        h = q-d
        value = density.engine.power(tuple(arb(float(v)) for v in mat.flat), 256)
        upper = sparse.rational((sum(value[:4], arb(0))*correction*arb(2)**(256*exponent)).upper())
        upper *= math.comb(spec['rows'], q)*math.comb(q, h)*len(bands)**d
        results.append(dict(all_one_rows=h, upper_power=sparse.dyadic_power(upper),
            diagnostic_margin_bits=math.log2(upper.denominator)-math.log2(upper.numerator)))
    core.require([r['all_one_rows'] for r in results] == list(range(q, -1, -1)), 'Missing constant-row case')
    return results


def run(output, m, q, witnesses, verify=False):
    saved = base.read(output) if verify else None
    if saved:
        core.authenticate(saved, base.ROOT)
        m, q = saved['instance']['message_exponent'], saved['occupation']
        witnesses = [(sh['tilt'], sh['anchor_all_one_rows']) for sh in saved['shards']]
    else:
        core.require(not output.exists(), 'Use a fresh output')
    spec = core.instance('t128_s19', m)
    core.require(2 <= q <= spec['rows'] and witnesses and
        all(type(tilt) is int and -120 <= tilt <= 30 and type(h) is int and 0 <= h < q for tilt, h in witnesses),
        'Invalid occupation/witness')
    if saved:
        core.require(spec == saved['instance'], 'Wrong instance')
    t, s, ac, kernel = core.inputs.load(spec['configuration'])
    caps = core.inputs.caps_module.caps()
    ctx.prec = 512 if verify else 256
    best, shards = {}, []
    for index, (tilt, h) in enumerate(witnesses):
        region = density.engine.regions(t, s, ac, kernel, (-(arb(tilt)/10).exp()).exp(), q, spec['rows'])
        ps = saved['shards'][index]['p'] if saved else density.choose(region[h:], q-h, caps)[:-1]
        rows = evaluate(region, ps, spec, q, tilt)
        if saved:
            old = saved['shards'][index]['rows']
            core.require([r['all_one_rows'] for r in old] == list(range(q, -1, -1)), 'Saved case gap')
            core.require(all(a['upper_power'] <= b['upper_power'] for a, b in zip(rows, old)), '512-bit replay failed')
            rows = old
        for row in rows:
            h1, power = row['all_one_rows'], row['upper_power']
            best[h1] = min(best.get(h1, power), power)
        shards.append(dict(tilt=tilt, anchor_all_one_rows=h, p=ps, rows=rows))
        print('Q', q, 'tilt/anchor', tilt, h, 'cases at least 70 bits', sum(p <= -70 for p in best.values()), '/', q+1, flush=True)
    total = sum((F(2)**best[h] for h in range(q+1)), F(0))
    status = 'COMPLETE_OCCUPANCY_CERTIFIED_50_BITS' if total <= F(1, 1 << 50) else 'COMPLETE_OCCUPANCY_WEAK_UPPER_BOUND'
    if saved:
        core.require([best[h] for h in range(q+1)] == saved['upper_powers'] and total == base.decode(saved['union_upper']) and
                     status == saved['status'], 'Aggregation mismatch')
        base.write_new(output.with_name(output.stem+'_replay.json'), dict(status='CONSTANT_SPLIT_512_BIT_REPLAY_PASSED',
            producer_sha256=base.sha(output), instance=spec, occupation=q, source_sha256=density.dependencies.sources()))
    else:
        base.write_new(output, dict(status=status, instance=spec, occupation=q, shards=shards,
            upper_powers=[best[h] for h in range(q+1)], union_upper=base.encode(total),
            source_sha256=density.dependencies.sources(), full_distance_proved=False))
    print(status, 'Q', q, 'margin', math.log2(total.denominator)-math.log2(total.numerator), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, default=16)
    p.add_argument('--occupation', type=int, default=512)
    p.add_argument('--witnesses', nargs='+', default=['6:0', '0:256', '-10:480'])
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output, a.m, a.occupation, [tuple(map(int, x.split(':'))) for x in a.witnesses], a.verify)

"""Batch sparse occupancies with directed positive folds and Arb terminals."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import numpy as np
from flint import arb, ctx
import sparse_bch as discovery
import scaled_adaptive as folds
model = discovery.model


def initial(region, n):
    mantissas, exponents = [], []
    for row in region:
        maximum = max(row)
        m, e = map(int, maximum.man_exp())
        exponent = e+m.bit_length() if m else 0
        mantissas.append([float(v*arb(2)**(-exponent)) for v in row])
        exponents.append(exponent)
    return np.nextafter(np.array(mantissas).reshape(-1, n, n), np.inf), np.array(exponents, dtype=np.int64)


def evaluate(engine, region, probabilities, tilt, lower, upper):
    assert len(region) == upper+1 and 2 <= lower <= upper <= engine.length
    costs = engine.costs(probabilities)
    roots = [model.up((cost.log()/256).exp()) for cost in costs]
    ps = list(map(model.number, probabilities))
    left = np.nextafter(np.array([float(model.up(r*(1-p))) for r, p in zip(roots, ps)]), np.inf)
    right = np.nextafter(np.array([float(model.up(r*p)) for r, p in zip(roots, ps)]), np.inf)
    correction = (engine.cutoff*model.number(tilt).exp()).exp()
    values = []
    for q, matrix, exponent in folds.matrices(*initial(region, engine.n), left, right):
        if q < lower:
            continue
        moment = model.independent.terminal(tuple(arb(float(v)) for v in matrix.flat), engine.n, 256)
        value = moment*arb(2)**(256*exponent)*correction*math.comb(engine.length, q)*len(probabilities)**q
        bits = model.exact((value.log()/arb(2).log()).upper())
        power = max(-200, -(-bits.numerator//bits.denominator))
        values.append(dict(occupation=q, power=power))
    assert [r['occupation'] for r in values] == list(range(lower, upper+1))
    return values


def run(output, exponent, last, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        exponent, last = saved['message_exponent'], saved['last']
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    engine = model.Engine(exponent)
    assert 2 <= last <= engine.length
    best, witnesses = {}, []
    if saved:
        assert saved['instance'] == engine.identity()
        for job in saved['witnesses']:
            tilt, (lo, hi) = F(job['log_tilt']), job['interval']
            ps = [model.base.decode(p) for p in job['probabilities']]
            rows = evaluate(engine, engine.region(tilt, hi), ps, tilt, lo, hi)
            assert len(rows) == len(job['rows'])
            for actual, old in zip(rows, job['rows']):
                assert actual['occupation'] == old['occupation'] and actual['power'] <= old['power']
                if old['power'] <= -60:
                    q = old['occupation']
                    best[q] = min(best.get(q, old['power']), old['power'])
            print('range replay', exponent, lo, hi, flush=True)
    else:
        while len(best) < last-1:
            anchor = next(q for q in range(2, last+1) if q not in best)
            hi = min(last, anchor+63, max(anchor+7, int(anchor*1.3)))
            prediction = (-76+6*math.log2(anchor)-10*(exponent-20)*math.log(2))/10
            for delta in (0., -.5, .5, -1., 1.):
                tilt = F.from_float(prediction+delta)
                region = engine.region(tilt, hi)
                ps = discovery.choose(engine, region[:anchor+1], anchor)
                rows = evaluate(engine, region, ps, tilt, anchor, hi)
                witnesses.append(dict(interval=[anchor, hi], log_tilt=str(tilt),
                                      probabilities=[model.base.encode(p) for p in ps], rows=rows))
                for row in rows:
                    if row['power'] <= -60:
                        q = row['occupation']
                        best[q] = min(best.get(q, row['power']), row['power'])
                print('range', exponent, anchor, hi, 'covered', len(best), 'anchor power', rows[0]['power'], flush=True)
                if anchor in best:
                    break
            if anchor not in best:
                break
            engine.region.cache_clear()
            engine.epoch.cache_clear()
    total = sum((F(2)**v for v in best.values()), F(0))
    if saved:
        assert best == {int(q):p for q, p in saved['best_powers'].items()}
        assert total == model.base.decode(saved['upper'])
        model.base.write_new(output.with_name(output.stem+'_replay.json'),
                             dict(status='ASYMMETRIC_SPARSE_RANGE_512_BIT_REPLAY_PASSED',
                                  producer_sha256=model.base.sha(output), full_distance_proved=False))
    else:
        complete = len(best) == last-1
        model.base.write_new(output, dict(status='OUTWARD_ASYMMETRIC_SPARSE_RANGE' if complete else 'OUTWARD_SPARSE_RANGE_UNRESOLVED',
                                         instance=engine.identity(), message_exponent=exponent, last=last,
                                         witnesses=witnesses, best_powers=best, complete_sparse_range=complete,
                                         upper=model.base.encode(total), source_sha256=model.sources(), full_distance_proved=False))
    print('ranges complete', exponent, len(best), '/', last-1, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=(16,18,20), default=20)
    p.add_argument('--last', type=int, default=511)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output.resolve(), a.m, a.last, a.verify)

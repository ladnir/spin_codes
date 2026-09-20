"""Recompute an IMT sparse range from length-scaled witness proposals."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import ladder
import sparse_ranges as transfer
from flint import ctx

model = ladder.model


def run(output, exponent, seed_path, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'IMT_LADDER_SPARSE_RANGE'
        exponent = saved['exponent']
    elif output.exists():
        raise FileExistsError('Use a fresh output path')
    ctx.prec = 512 if verify else 256
    engine = ladder.Engine(exponent)
    if saved:
        assert saved['instance'] == engine.identity()
        jobs, last = saved['witnesses'], saved['last']
    else:
        prior = model.base.read(seed_path)
        model.authenticate(prior)
        assert prior['candidate'] == 'weight5_seed0' and prior['complete_sparse_range']
        last = prior['last']
        # Only the proposal changes. Every transfer, count, and retained bound
        # is recomputed at the target length with outward arithmetic.
        shift = F.from_float((exponent - prior['exponent']) * math.log(2))
        jobs = [dict(interval=job['interval'], tilt=str(F(job['tilt']) - shift),
                     probabilities=job['probabilities']) for job in prior['witnesses']]
    assert 2 <= last <= engine.length
    best, witnesses = {}, []
    for job in jobs:
        lo, hi = job['interval']
        assert 2 <= lo <= hi <= last
        tilt = F(job['tilt'])
        probabilities = list(map(model.base.decode, job['probabilities']))
        rows = transfer.evaluate(engine, engine.region(tilt, hi), probabilities, tilt, lo, hi)
        if saved:
            assert len(rows) == len(job['rows'])
            for current, old in zip(rows, job['rows']):
                assert current['occupation'] == old['occupation'] and current['power'] <= old['power']
            rows = job['rows']
        witnesses.append(dict(interval=[lo, hi], tilt=str(tilt),
                              probabilities=job['probabilities'], rows=rows))
        for row in rows:
            q, power = row['occupation'], row['power']
            if power <= -60:
                best[q] = min(best.get(q, power), power)
        engine.region.cache_clear()
        engine.epoch.cache_clear()
        print('sparse', ctx.prec, exponent, lo, hi, 'covered', len(best), '/', last - 1, flush=True)
    complete = set(best) == set(range(2, last + 1))
    total = sum((F(2)**p for p in best.values()), F(0))
    if saved:
        assert complete == saved['complete_sparse_range']
        assert best == {int(q): p for q, p in saved['best'].items()}
        assert total == model.base.decode(saved['upper'])
        model.base.write_new(output.with_name(output.stem + '_replay.json'), dict(
            status='IMT_LADDER_SPARSE_512_BIT_REPLAY_PASSED',
            producer_sha256=model.base.sha(output), complete_sparse_range=complete,
            full_distance_proved=False))
    else:
        sources = ladder.candidate.sources()
        sources[seed_path.relative_to(model.ROOT).as_posix()] = model.base.sha(seed_path)
        model.base.write_new(output, dict(status='IMT_LADDER_SPARSE_RANGE',
            exponent=exponent, instance=engine.identity(), last=last,
            witnesses=witnesses, best=best, upper=model.base.encode(total),
            complete_sparse_range=complete, full_distance_proved=False,
            source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=ladder.EXPONENTS, default=22)
    p.add_argument('--seed', type=Path)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    if not a.verify and not a.seed:
        p.error('Producer requires a seed bank')
    run(a.output.resolve(), a.m, a.seed.resolve() if a.seed else None, a.verify)

"""Size-specific split search with compact receipts and assigned-case replay.

Only witness locations are scaled from m16. Every coefficient, probability,
case bound, and exact union is recomputed for the requested instance.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import time

import numpy as np
from flint import arb, ctx

import scalable_split_grid as fast
import fine_split_search as prior

core, base, density = fast.core, fast.base, fast.density
PRODUCER = 'SCALED_SPLIT_SELECTED_PRODUCER_V1'
REPLAY = 'SCALED_SPLIT_SELECTED_512_BIT_REPLAY_PASSED'


def encode_table(table, lo, hi):
    result = []
    for q in range(lo, hi+1):
        row = table[q-lo, :q+1]
        starts = np.concatenate(([0], np.flatnonzero(row[1:] != row[:-1])+1, [q+1]))
        result.append([[int(row[a]), int(b-a)] for a, b in zip(starts[:-1], starts[1:])])
    return result


def decode_table(saved, lo, hi, low, high):
    core.require(len(saved) == hi-lo+1, 'Wrong RLE row count')
    table = fast.table(lo, hi)
    for q, runs in zip(range(lo, hi+1), saved):
        offset = 0
        for run in runs:
            core.require(len(run) == 2, 'Invalid run')
            value, count = run
            core.require(type(value) is int and type(count) is int and low <= value <= high and
                1 <= count <= q+1 and offset+count <= q+1, 'Invalid RLE value/count')
            table[q-lo, offset:offset+count] = value
            offset += count
        core.require(offset == q+1, 'Missing cases in RLE row')
    return table


def passing(best, lo, hi):
    values = {}
    for q in range(lo, hi+1):
        row = best[q-lo, :q+1]
        if int(row.max()) > -50:
            continue
        counts = np.bincount(row+80, minlength=31)
        value = F(sum(int(n) << i for i, n in enumerate(counts)), 1 << 80)
        if value <= F(2)**-50:
            values[q] = value
    return values


def proposal(best, lo, hi, shards, factor):
    good = passing(best, lo, hi)
    gaps = core.intervals(q for q in range(lo, hi+1) if q not in good)
    tried = {(F(v['tilt_numerator'], v['tilt_denominator']), v['anchor_occupation'], v['anchor_all_one_rows']) for v in shards}
    anchors = [(150, -6), (200, -4), (256, -2), (384, 1), (512, 8)]
    for left, right in sorted(gaps, key=lambda x: (-(x[1]-x[0]), x[0])):
        q = (left+right)//2
        h = min(int(best[q-lo, :q+1].argmax()), q-1)
        x = F(q, factor)
        lower = max((a for a in anchors if a[0] <= x), default=anchors[0])
        upper = min((a for a in anchors if a[0] >= x), default=anchors[-1])
        center = F(lower[1]) if lower[0] == upper[0] else lower[1]+(x-lower[0])*F(upper[1]-lower[1], upper[0]-lower[0])
        center = round(4*center)
        for delta in (0, -1, 1, -2, 2, -4, 4, -8, 8, -12, 12):
            candidate = (F(center+delta, 4), q, h)
            if candidate not in tried:
                return candidate
    return None


def run(output, m=18, lo=480, hi=2048, jobs=28, seconds=600, resume=None):
    core.require(not output.exists() and jobs > 0 and seconds > 0, 'Invalid output/budget')
    spec = core.instance('t128_s19', m)
    core.require(16 <= m <= 20 and 1 <= lo <= hi <= spec['rows'] and
        (hi-lo+1)*(hi+1) <= 10**7, 'Invalid or oversized grid; use a smaller interval')
    factor = 1 << (m-16)
    if resume:
        saved, best, owners = load(resume, require_replay=True)
        core.require(saved['instance'] == spec and saved['interval'] == [lo, hi], 'Wrong resume instance')
        shards = saved['shards'].copy()
        candidates = []
    else:
        best, owners = fast.table(lo, hi), np.full((hi-lo+1, hi+1), -1, dtype=np.int64)
        shards = []
        seed_path = base.HERE/'generated/t128_s19_m16_fine_split_grid_v1.json'
        seed = prior.load_verified(seed_path)
        candidates = [(F(-8), min(hi, 125*factor), 0)]
        candidates += [(F(v['tilt_numerator'], v['tilt_denominator']), min(hi, v['anchor_occupation']*factor),
                        min(hi-1, v['anchor_all_one_rows']*factor)) for v in seed['shards']]
        candidates = list(dict.fromkeys(candidates))
    t, s, ac, kernel = core.inputs.load(spec['configuration'])
    caps = core.inputs.caps_module.caps()
    ctx.prec = 256
    started = time.monotonic()
    for job in range(jobs):
        good = passing(best, lo, hi)
        if len(good) == hi-lo+1 or time.monotonic()-started >= seconds:
            break
        candidate = candidates.pop(0) if candidates else proposal(best, lo, hi, shards, factor)
        if candidate is None:
            break
        index, q, h = candidate
        tilt = arb(index.numerator)/index.denominator
        print('start witness', len(shards)+1, (str(index), q, h), flush=True)
        region = density.engine.regions(t, s, ac, kernel, (-(tilt/10).exp()).exp(), hi, spec['rows'])
        ps = density.choose(region[h:q+1], q-h, caps)[:-1]
        count = fast.evaluate(region, ps, spec, lo, hi, tilt, best, owners, len(shards))
        shards.append(dict(tilt_numerator=index.numerator, tilt_denominator=index.denominator,
            anchor_occupation=q, anchor_all_one_rows=h, p=ps))
        print('evaluated', count, 'cases; coverage', core.intervals(passing(best, lo, hi)),
            'seconds', round(time.monotonic()-started, 2), flush=True)
    core.require(shards and all(int(best[q-lo, :q+1].max()) < fast.UNSET for q in range(lo, hi+1)), 'No complete producer sweep')
    result = dict(status=PRODUCER, instance=spec, interval=[lo, hi], policy=fast.frozen.POLICY,
        shards=shards, powers_rle=encode_table(best, lo, hi), owners_rle=encode_table(owners, lo, hi),
        per_occupancy_50_bit_intervals=core.intervals(passing(best, lo, hi)),
        source_sha256=density.dependencies.sources(), full_distance_proved=False)
    base.write_new(output, result)
    print(PRODUCER, 'coverage', result['per_occupancy_50_bit_intervals'], flush=True)
    return result


def load(path, require_replay=False):
    saved = base.read(path)
    core.authenticate(saved, base.ROOT)
    core.require(saved['status'] == PRODUCER and saved['policy'] == fast.frozen.POLICY and saved['shards'], 'Wrong producer')
    lo, hi = saved['interval']
    core.require(1 <= lo <= hi <= saved['instance']['rows'] and (hi-lo+1)*(hi+1) <= 10**7, 'Invalid grid')
    best = decode_table(saved['powers_rle'], lo, hi, -80, 10**8)
    owners = decode_table(saved['owners_rle'], lo, hi, 0, len(saved['shards'])-1)
    core.require(core.intervals(passing(best, lo, hi)) == saved['per_occupancy_50_bit_intervals'], 'Wrong coverage')
    if require_replay:
        receipt = base.read(path.with_name(path.stem+'_replay.json'))
        core.authenticate(receipt, base.ROOT)
        core.require(receipt['status'] == REPLAY and receipt['producer_sha256'] == base.sha(path) and
            receipt['instance'] == saved['instance'] and receipt['interval'] == saved['interval'] and
            receipt['cases_checked'] == sum(q+1 for q in range(lo, hi+1)), 'Wrong replay')
    return saved, best, owners


def verify(path):
    saved, best, owners = load(path)
    spec = core.instance(saved['instance']['configuration'], saved['instance']['message_exponent'])
    core.require(spec == saved['instance'], 'Wrong instance')
    lo, hi = saved['interval']
    t, s, ac, kernel = core.inputs.load(spec['configuration'])
    ctx.prec = 512
    visited = 0
    for i, v in enumerate(saved['shards']):
        n, d = v['tilt_numerator'], v['tilt_denominator']
        core.require(type(n) is int and type(d) is int and 1 <= d <= 10000 and -120 <= F(n, d) <= 30,
                     'Invalid tilt')
        if not bool((owners == i).any()):
            continue
        tilt = arb(n)/d
        print('replay witness', i+1, '/', len(saved['shards']), flush=True)
        region = density.engine.regions(t, s, ac, kernel, (-(tilt/10).exp()).exp(), hi, spec['rows'])
        count = fast.evaluate(region, v['p'], spec, lo, hi, tilt, best, owners, i, replay=True)
        visited += count
        print('checked', count, 'cases; total', visited, flush=True)
    core.require(visited == sum(q+1 for q in range(lo, hi+1)), 'Missing or repeated replay cases')
    result = dict(status=REPLAY, producer_sha256=base.sha(path), instance=spec, interval=[lo, hi],
        cases_checked=visited, source_sha256=density.dependencies.sources())
    base.write_new(path.with_name(path.stem+'_replay.json'), result)
    print(REPLAY, flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, default=18)
    p.add_argument('--first', type=int, default=480)
    p.add_argument('--last', type=int, default=2048)
    p.add_argument('--max-jobs', type=int, default=28)
    p.add_argument('--seconds', type=float, default=600)
    p.add_argument('--resume', type=Path)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    if a.verify:
        verify(a.output)
    else:
        run(a.output, a.m, a.first, a.last, a.max_jobs, a.seconds, a.resume)

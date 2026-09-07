"""Bounded gap-driven continuation of the default all-one-split search.

Resume from an authenticated, replayed grid, propose witnesses at actual gaps,
and emit a new immutable bank. Replay still recomputes every retained witness.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import time

from flint import arb, ctx

import search_split_certificates as search

grid, core, base, density = search.grid, search.core, search.base, search.density


def propose(powers, lo, hi, shards):
    gaps = core.intervals(q for q, b in zip(range(lo, hi+1), grid.bounds(powers, lo, hi)) if b > F(2)**-50)
    if not gaps:
        return None
    attempted = {(v['tilt'], v['anchor_occupation'], v['anchor_all_one_rows']) for v in shards}
    # Purely a proposal rule: interpolate successful initial h=0 anchor tilts.
    # No bounds or coverage are interpolated, and every new p is recomputed.
    anchors = sorted((q, tilt) for tilt, q, h in search.DEFAULT_WITNESSES if h == 0)
    for left, right in sorted(gaps, key=lambda x: (-(x[1]-x[0]), x[0])):
        q = (left+right)//2
        h = min(max(range(q+1), key=lambda h1: powers[q-lo][h1]), q-1)
        low = max((a for a in anchors if a[0] <= q), default=anchors[0])
        high = min((a for a in anchors if a[0] >= q), default=anchors[-1])
        center = low[1] if low[0] == high[0] else round(low[1]+(q-low[0])*(high[1]-low[1])/(high[0]-low[0]))
        for shift in (0, -1, 1, -2, 2, -3, 3):
            candidate = (center+shift, q, h)
            if candidate not in attempted:
                return candidate
    return None


def run(seed, output, max_jobs=16, seconds=180):
    core.require(not output.exists() and type(max_jobs) is int and 1 <= max_jobs <= 100 and seconds > 0,
                 'Invalid budget or existing output')
    saved = base.read(seed)
    receipt_path = seed.with_name(seed.stem+'_replay.json')
    receipt = base.read(receipt_path)
    core.authenticate(saved, base.ROOT)
    core.authenticate(receipt, base.ROOT)
    core.require(saved['status'] == 'CONSTANT_SPLIT_GRID_OUTWARD_PRODUCER' and
        receipt['status'] == 'CONSTANT_SPLIT_GRID_512_BIT_REPLAY_PASSED' and
        receipt['producer_sha256'] == base.sha(seed) and
        saved['policy'] == grid.POLICY == receipt['policy'] and saved['instance'] == receipt['instance'] and
        saved['interval'] == receipt['interval'], 'Unreplayed or mismatched seed')
    spec = core.instance(saved['instance']['configuration'], saved['instance']['message_exponent'])
    core.require(spec == saved['instance'] and spec['message_exponent'] == 16,
                 'This proposal policy is currently tuned for m16 only')
    lo, hi = saved['interval']
    core.require(receipt['cases_checked'] == sum(q+1 for q in range(lo, hi+1)) and
        search.passing(saved['upper_powers'], lo, hi) == saved['per_occupancy_50_bit_intervals'], 'Wrong seed coverage')
    best = [row.copy() for row in saved['upper_powers']]
    shards = saved['shards'].copy()
    t, s, ac, kernel = core.inputs.load(spec['configuration'])
    caps = core.inputs.caps_module.caps()
    ctx.prec = 256
    started = time.monotonic()
    for job in range(max_jobs):
        candidate = propose(best, lo, hi, shards)
        if candidate is None or time.monotonic()-started >= seconds:
            break
        tilt, q, h = candidate
        region = density.engine.regions(t, s, ac, kernel, (-(arb(tilt)/10).exp()).exp(), hi, spec['rows'])
        ps = density.choose(region[h:q+1], q-h, caps)[:-1]
        grid.merge(best, grid.evaluate(region, ps, spec, lo, hi, tilt))
        shards.append(dict(tilt=tilt, anchor_occupation=q, anchor_all_one_rows=h, p=ps))
        print('gap witness', job+1, candidate, '50-bit occupancies', search.passing(best, lo, hi),
              'seconds', round(time.monotonic()-started, 2), flush=True)
    sources = density.dependencies.sources()
    for path in (seed, receipt_path):
        sources[path.resolve().relative_to(base.ROOT).as_posix()] = base.sha(path)
    result = dict(status='CONSTANT_SPLIT_GRID_OUTWARD_PRODUCER', instance=spec, interval=[lo, hi],
        policy=grid.POLICY, shards=shards, upper_powers=best,
        per_occupancy_50_bit_intervals=search.passing(best, lo, hi), source_sha256=sources,
        full_distance_proved=False, continuation_of=seed.resolve().relative_to(base.ROOT).as_posix())
    base.write_new(output, result)
    print('Search complete; output still requires 512-bit replay:', output, flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--max-jobs', type=int, default=16)
    p.add_argument('--seconds', type=float, default=180)
    a = p.parse_args()
    run(a.seed, a.output, a.max_jobs, a.seconds)

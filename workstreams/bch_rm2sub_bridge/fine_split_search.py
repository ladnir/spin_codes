"""Gap-driven all-one-split search with exact rational tilt indices.

The original index step was 1 (a 0.1 step in log(lambda)). Quarter-index
proposals avoid missing narrow passing windows. Every resulting (Q,h) bound
is evaluated outward; neither witness interpolation nor the proposal is proof.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import time

from flint import arb, ctx

import search_split_certificates as coarse

grid, core, base, density = coarse.grid, coarse.core, coarse.base, coarse.density
PRODUCER = 'FINE_CONSTANT_SPLIT_GRID_OUTWARD_PRODUCER'
REPLAY = 'FINE_CONSTANT_SPLIT_GRID_512_BIT_REPLAY_PASSED'


def load_verified(path):
    saved = base.read(path)
    replay_path = path.with_name(path.stem+'_replay.json')
    replay = base.read(replay_path)
    core.authenticate(saved, base.ROOT)
    core.authenticate(replay, base.ROOT)
    expected = {'CONSTANT_SPLIT_GRID_OUTWARD_PRODUCER': 'CONSTANT_SPLIT_GRID_512_BIT_REPLAY_PASSED', PRODUCER: REPLAY}
    lo, hi = saved['interval']
    core.require(saved['status'] in expected and replay['status'] == expected[saved['status']] and
        replay['producer_sha256'] == base.sha(path) and saved['instance'] == replay['instance'] and
        saved['policy'] == grid.POLICY == replay['policy'] and saved['interval'] == replay['interval'] and
        replay['cases_checked'] == sum(q+1 for q in range(lo, hi+1)) and
        coarse.passing(saved['upper_powers'], lo, hi) == saved['per_occupancy_50_bit_intervals'], 'Invalid seed/replay')
    return saved


def normal_shards(saved):
    if saved['status'] == PRODUCER:
        return saved['shards'].copy()
    return [dict(tilt_numerator=v['tilt'], tilt_denominator=1, anchor_occupation=v['anchor_occupation'],
                 anchor_all_one_rows=v['anchor_all_one_rows'], p=v['p']) for v in saved['shards']]


def propose(powers, lo, hi, shards):
    gaps = core.intervals(q for q, b in zip(range(lo, hi+1), grid.bounds(powers, lo, hi)) if b > F(2)**-50)
    attempted = {(F(v['tilt_numerator'], v['tilt_denominator']), v['anchor_occupation'], v['anchor_all_one_rows']) for v in shards}
    anchors = sorted((q, tilt) for tilt, q, h in coarse.DEFAULT_WITNESSES if h == 0)
    for left, right in sorted(gaps, key=lambda x: (-(x[1]-x[0]), x[0])):
        q = (left+right)//2
        h = min(max(range(q+1), key=lambda h1: powers[q-lo][h1]), q-1)
        low = max((a for a in anchors if a[0] <= q), default=anchors[0])
        high = min((a for a in anchors if a[0] >= q), default=anchors[-1])
        center = F(low[1]) if low[0] == high[0] else low[1]+F((q-low[0])*(high[1]-low[1]), high[0]-low[0])
        center = round(4*center)
        for shift in (0, -1, 1, -2, 2, -3, 3, -4, 4, -6, 6, -8, 8):
            candidate = (F(center+shift, 4), q, h)
            if candidate not in attempted:
                return candidate
    return None


def run(seed, output, max_jobs=16, seconds=180):
    core.require(not output.exists() and 1 <= max_jobs <= 100 and seconds > 0, 'Invalid output/budget')
    saved = load_verified(seed)
    spec = core.instance(saved['instance']['configuration'], saved['instance']['message_exponent'])
    core.require(spec == saved['instance'] and spec['message_exponent'] == 16 and spec['configuration'] == 't128_s19',
                 'Automatic proposal policy is currently tuned for t128_s19/m16 only')
    lo, hi = saved['interval']
    best, shards = [row.copy() for row in saved['upper_powers']], normal_shards(saved)
    t, s, ac, kernel = core.inputs.load(spec['configuration'])
    caps = core.inputs.caps_module.caps()
    ctx.prec = 256
    started = time.monotonic()
    for job in range(max_jobs):
        candidate = propose(best, lo, hi, shards)
        if candidate is None or time.monotonic()-started >= seconds:
            break
        index, q, h = candidate
        tilt = arb(index.numerator)/index.denominator
        region = density.engine.regions(t, s, ac, kernel, (-(tilt/10).exp()).exp(), hi, spec['rows'])
        ps = density.choose(region[h:q+1], q-h, caps)[:-1]
        grid.merge(best, grid.evaluate(region, ps, spec, lo, hi, tilt))
        shards.append(dict(tilt_numerator=index.numerator, tilt_denominator=index.denominator,
                           anchor_occupation=q, anchor_all_one_rows=h, p=ps))
        print('fine gap witness', job+1, (str(index), q, h), '50-bit occupancies', coarse.passing(best, lo, hi),
              'seconds', round(time.monotonic()-started, 2), flush=True)
    sources = density.dependencies.sources()
    for path in (seed, seed.with_name(seed.stem+'_replay.json')):
        sources[path.resolve().relative_to(base.ROOT).as_posix()] = base.sha(path)
    result = dict(status=PRODUCER, instance=spec, interval=[lo, hi], policy=grid.POLICY, shards=shards,
        upper_powers=best, per_occupancy_50_bit_intervals=coarse.passing(best, lo, hi), source_sha256=sources,
        full_distance_proved=False, continuation_of=seed.resolve().relative_to(base.ROOT).as_posix())
    base.write_new(output, result)
    print('Producer complete; replay required:', output, flush=True)
    return result


def verify(path):
    saved = base.read(path)
    core.authenticate(saved, base.ROOT)
    core.require(saved['status'] == PRODUCER and saved['policy'] == grid.POLICY, 'Wrong producer/policy')
    spec = core.instance(saved['instance']['configuration'], saved['instance']['message_exponent'])
    core.require(spec == saved['instance'] and saved['shards'], 'Wrong instance or empty bank')
    lo, hi = saved['interval']
    core.require(1 <= lo <= hi <= spec['rows'], 'Wrong interval')
    grid.bounds(saved['upper_powers'], lo, hi)
    t, s, ac, kernel = core.inputs.load(spec['configuration'])
    best = [[None]*(q+1) for q in range(lo, hi+1)]
    ctx.prec = 512
    for i, v in enumerate(saved['shards']):
        n, d = v['tilt_numerator'], v['tilt_denominator']
        core.require(type(n) is int and type(d) is int and 1 <= d <= 10000 and -120 <= F(n, d) <= 30 and
            type(v['anchor_occupation']) is int and type(v['anchor_all_one_rows']) is int and
            0 <= v['anchor_all_one_rows'] < v['anchor_occupation'] <= hi, 'Invalid rational witness')
        tilt = arb(n)/d
        region = density.engine.regions(t, s, ac, kernel, (-(tilt/10).exp()).exp(), hi, spec['rows'])
        grid.merge(best, grid.evaluate(region, v['p'], spec, lo, hi, tilt))
        print('512-bit fine replay', i+1, '/', len(saved['shards']), coarse.passing(best, lo, hi), flush=True)
    core.require(all(all(a <= b for a, b in zip(new, old)) for new, old in zip(best, saved['upper_powers'])),
                 '512-bit replay exceeded retained case bound')
    core.require(coarse.passing(saved['upper_powers'], lo, hi) == saved['per_occupancy_50_bit_intervals'], 'Wrong coverage')
    result = dict(status=REPLAY, producer_sha256=base.sha(path), instance=spec, interval=[lo, hi],
        policy=grid.POLICY, cases_checked=sum(q+1 for q in range(lo, hi+1)), source_sha256=density.dependencies.sources())
    base.write_new(path.with_name(path.stem+'_replay.json'), result)
    print(REPLAY, flush=True)
    return result


def combine(output, report_path, paths, bits=40):
    core.require(not output.exists() and paths, 'Existing output or no grid')
    coarse.verifier.verify_report(report_path)
    report = base.read(report_path)
    spec, best, sources = report['instance'], {}, {}
    def record(path):
        sources[path.resolve().relative_to(base.ROOT).as_posix()] = base.sha(path)
    record(report_path)
    for item in report['accepted_certificates']:
        cert, replay = report_path.parent/item['certificate'], report_path.parent/item['replay']
        values = coarse.verifier.load_verified(cert, replay, spec)
        core.merge_rows(best, [v for v in values if v['occupation'] in item['occupancies']], spec['rows'])
        record(cert)
        record(replay)
    for path in paths:
        saved = load_verified(path)
        core.require(saved['instance'] == spec, 'Wrong grid instance')
        lo, hi = saved['interval']
        core.merge_rows(best, [dict(occupation=q, upper=base.encode(b)) for q, b in
            zip(range(lo, hi+1), grid.bounds(saved['upper_powers'], lo, hi)) if b <= F(2)**-50], spec['rows'])
        record(path)
        record(path.with_name(path.stem+'_replay.json'))
    sources.update(density.dependencies.sources())
    result = dict(status='EXACT_FINE_SPLIT_COVERAGE_LEDGER', instance=spec, coverage=core.summarize(best, spec['rows'], bits),
        per_occupancy_bounds=[dict(occupation=q, upper=base.encode(b)) for q, b in sorted(best.items())],
        source_sha256=sources, numerical_replays_rerun_by_this_ledger=False)
    base.write_new(output, result)
    print(result['coverage'], flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--verify', action='store_true')
    p.add_argument('--max-jobs', type=int, default=16)
    p.add_argument('--seconds', type=float, default=180)
    p.add_argument('--sparse-report', type=Path)
    p.add_argument('--grid', type=Path, action='append')
    p.add_argument('--target-bits', type=int, default=40)
    a = p.parse_args()
    if a.sparse_report:
        core.require(not a.verify and not a.seed, 'Choose one action')
        combine(a.output, a.sparse_report, a.grid, a.target_bits)
    elif a.verify:
        core.require(not a.seed, 'Choose one action')
        verify(a.output)
    else:
        core.require(a.seed, 'Continuation requires a replayed seed')
        run(a.seed, a.output, a.max_jobs, a.seconds)

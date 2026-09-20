"""Default dense BCH-256 search: all-one rows are always counted separately.

Produce a bounded witness-bank run, replay at 512 bits, and combine with an
authenticated sparse report. Existing frozen producers remain unchanged.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import time

from flint import arb, ctx

import constant_split_grid as grid
import verify_certificate_search as verifier

core, base, density = grid.core, grid.base, grid.density
DEFAULT_WITNESSES = [(-6, 150, 0), (-4, 200, 0), (-2, 256, 0), (1, 384, 0),
                     (8, 512, 0), (6, 512, 64), (4, 512, 128), (2, 512, 160), (0, 512, 256)]


def passing(powers, lo, hi):
    values = grid.bounds(powers, lo, hi)
    return core.intervals(q for q, value in zip(range(lo, hi+1), values) if value <= F(2)**-50)


def run(output, m=16, lo=150, hi=512, witnesses=None, verify=False):
    saved = base.read(output) if verify else None
    if saved is not None:
        core.authenticate(saved, base.ROOT)
        core.require(saved['status'] == 'CONSTANT_SPLIT_GRID_OUTWARD_PRODUCER' and
                     saved['policy'] == grid.POLICY, 'Wrong grid producer/policy')
        m = saved['instance']['message_exponent']
        lo, hi = saved['interval']
        witnesses = [(v['tilt'], v['anchor_occupation'], v['anchor_all_one_rows']) for v in saved['shards']]
        grid.bounds(saved['upper_powers'], lo, hi)
    else:
        core.require(not output.exists(), 'Use a fresh output path')
        if witnesses is None:
            core.require(m == 16 and hi == 512, 'Supply an explicit witness bank for other sizes')
            witnesses = DEFAULT_WITNESSES
    spec = core.instance('t128_s19', m)
    core.require(1 <= lo <= hi <= spec['rows'] and witnesses and
        all(type(tilt) is int and -120 <= tilt <= 30 and type(q) is int and type(h) is int and
            1 <= q <= hi and 0 <= h < q for tilt, q, h in witnesses), 'Invalid grid/witness bank')
    if saved is not None:
        core.require(saved['instance'] == spec, 'Wrong instance')
    t, s, ac, kernel = core.inputs.load(spec['configuration'])
    caps = core.inputs.caps_module.caps()
    ctx.prec = 512 if verify else 256
    best = [[None]*(q+1) for q in range(lo, hi+1)]
    shards = []
    started = time.monotonic()
    for index, (tilt, q, h) in enumerate(witnesses):
        region = density.engine.regions(t, s, ac, kernel, (-(arb(tilt)/10).exp()).exp(), hi, spec['rows'])
        ps = saved['shards'][index]['p'] if saved is not None else density.choose(region[h:q+1], q-h, caps)[:-1]
        values = grid.evaluate(region, ps, spec, lo, hi, tilt)
        grid.merge(best, values)
        shards.append(dict(tilt=tilt, anchor_occupation=q, anchor_all_one_rows=h, p=ps))
        print(ctx.prec, 'witness', index+1, '/', len(witnesses), (tilt, q, h),
              '50-bit occupancies', passing(best, lo, hi), 'seconds', round(time.monotonic()-started, 2), flush=True)
    if saved is not None:
        core.require(all(all(new <= old for new, old in zip(a, b)) for a, b in zip(best, saved['upper_powers'])),
                     '512-bit replay exceeded retained case bound')
        core.require(passing(saved['upper_powers'], lo, hi) == saved['per_occupancy_50_bit_intervals'],
                     'Wrong saved coverage')
        result = dict(status='CONSTANT_SPLIT_GRID_512_BIT_REPLAY_PASSED', producer_sha256=base.sha(output),
            instance=spec, interval=[lo, hi], policy=grid.POLICY,
            cases_checked=sum(q+1 for q in range(lo, hi+1)), source_sha256=density.dependencies.sources())
        base.write_new(output.with_name(output.stem+'_replay.json'), result)
    else:
        result = dict(status='CONSTANT_SPLIT_GRID_OUTWARD_PRODUCER', instance=spec,
            interval=[lo, hi], policy=grid.POLICY, shards=shards, upper_powers=best,
            per_occupancy_50_bit_intervals=passing(best, lo, hi), source_sha256=density.dependencies.sources(),
            full_distance_proved=False)
        base.write_new(output, result)
    print(result['status'], flush=True)
    return result


def combine(output, sparse_report, producers, bits=40):
    core.require(not output.exists(), 'Use a fresh ledger path')
    verifier.verify_report(sparse_report)
    report = base.read(sparse_report)
    spec = report['instance']
    best, dependencies, accepted = {}, {}, []
    def record(path):
        dependencies[path.resolve().relative_to(base.ROOT).as_posix()] = base.sha(path)
    record(sparse_report)
    for item in report['accepted_certificates']:
        cert, replay = sparse_report.parent/item['certificate'], sparse_report.parent/item['replay']
        values = verifier.load_verified(cert, replay, spec)
        core.merge_rows(best, [v for v in values if v['occupation'] in item['occupancies']], spec['rows'])
        record(cert)
        record(replay)
    for path in producers:
        saved = base.read(path)
        replay_path = path.with_name(path.stem+'_replay.json')
        replay = base.read(replay_path)
        core.authenticate(saved, base.ROOT)
        core.authenticate(replay, base.ROOT)
        lo, hi = saved['interval']
        core.require(saved['status'] == 'CONSTANT_SPLIT_GRID_OUTWARD_PRODUCER' and
            replay['status'] == 'CONSTANT_SPLIT_GRID_512_BIT_REPLAY_PASSED' and
            saved['instance'] == spec == replay['instance'] and replay['interval'] == [lo, hi] and
            saved['policy'] == grid.POLICY == replay['policy'] and
            replay['producer_sha256'] == base.sha(path) and
            replay['cases_checked'] == sum(q+1 for q in range(lo, hi+1)), 'Wrong grid/replay')
        values = grid.bounds(saved['upper_powers'], lo, hi)
        core.require(passing(saved['upper_powers'], lo, hi) == saved['per_occupancy_50_bit_intervals'], 'Wrong coverage')
        selected = [dict(occupation=q, upper=base.encode(value)) for q, value in
            zip(range(lo, hi+1), values) if value <= F(2)**-50]
        core.merge_rows(best, selected, spec['rows'])
        accepted.append(dict(certificate=path.resolve().relative_to(base.ROOT).as_posix(),
                             occupancies=[r['occupation'] for r in selected]))
        record(path)
        record(replay_path)
    record(Path(__file__))
    record(Path(grid.__file__))
    result = dict(status='EXACT_SPLIT_GRID_COVERAGE_LEDGER', instance=spec, coverage=core.summarize(best, spec['rows'], bits),
        accepted_grids=accepted, per_occupancy_bounds=[dict(occupation=q, upper=base.encode(value)) for q, value in sorted(best.items())],
        source_sha256=dependencies, numerical_replays_rerun_by_this_ledger=False)
    base.write_new(output, result)
    print(result['coverage'], flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--m', type=int, default=16)
    parser.add_argument('--first', type=int, default=150)
    parser.add_argument('--last', type=int, default=512)
    parser.add_argument('--witness', action='append', help='tilt:anchor-Q:anchor-h; negative values use --witness=-6:150:0')
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--sparse-report', type=Path)
    parser.add_argument('--grid', type=Path, action='append')
    parser.add_argument('--target-bits', type=int, default=40)
    args = parser.parse_args()
    if args.sparse_report:
        core.require(args.grid and not args.verify, 'Ledger requires replayed grids')
        combine(args.output, args.sparse_report, args.grid, args.target_bits)
    else:
        witnesses = [tuple(map(int, v.split(':'))) for v in args.witness] if args.witness else None
        run(args.output, args.m, args.first, args.last, witnesses, args.verify)

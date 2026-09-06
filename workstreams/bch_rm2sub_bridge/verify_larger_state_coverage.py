"""Exact coverage ledger; excludes nonpassing rows and never fills gaps."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import bridge as base


def verify(tag, range_tags):
    assert tag.isidentifier() and all(name.isidentifier() for name in range_tags)
    q1_path = base.HERE/'generated/larger_t64_s20_q1_outward.json'
    q1_replay = base.HERE/'generated/larger_t64_s20_q1_replay.json'
    paths = [q1_path, q1_replay]
    q1 = base.read(q1_path)
    assert q1['configuration'] == 't64_s20'
    replay = base.read(q1_replay)
    assert replay['status'] == 'INDEPENDENT_512_BIT_Q1_REPLAY_PASSED'
    assert replay['producer_sha256'] == base.sha(q1_path) and replay['coefficients_checked'] == 92
    best = {1:base.decode(q1['Q1_upper'])}
    for source, digest in q1['local_sha256'].items():
        assert base.sha(base.HERE/source) == digest
    for name in range_tags:
        path = base.HERE/'generated'/f'larger_range_{name}_outward.json'
        replay_path = path.with_name(f'larger_range_{name}_replay.json')
        paths += [path, replay_path]
        saved = base.read(path)
        replay = base.read(replay_path)
        assert saved['configuration'] == 't64_s20'
        assert replay['status'] == 'LARGER_STATE_512_BIT_RANGE_REPLAY_PASSED'
        assert replay['producer_sha256'] == base.sha(path)
        lower, upper = saved['occupancy_range']
        assert len(saved['rows']) == upper-lower+1 == replay['rows_checked']
        for source, digest in saved['local_sha256'].items():
            assert base.sha(base.HERE/source) == digest
        for q, row in zip(range(lower, upper+1), saved['rows']):
            assert row['occupation'] == q
            bound = base.decode(row['upper'])
            assert bound > 0
            if bound < F(1,1<<60):
                best[q] = min(best.get(q,bound), bound)
    intervals = []
    for q in sorted(best):
        if intervals and q == intervals[-1][1]+1:
            intervals[-1][1] = q
        else:
            intervals.append([q,q])
    total = sum(best.values(), F(0))
    full = set(best) == set(range(1,8193))
    assert total < F(1,1<<50)
    output = base.HERE/'generated'/f'larger_coverage_{tag}.json'
    base.write_new(output, dict(status='EXACT_LARGER_STATE_COVERAGE_LEDGER', configuration='t64_s20',
        covered_intervals=intervals, covered_occupancies=len(best), covered_union_upper=base.encode(total),
        covered_union_below_2_to_minus_50=True, all_occupations_certified=full,
        full_target_proved=full and total<F(1,1<<40),
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__)]+paths}))
    print('Exact covered intervals', intervals, 'union <2^-50; full target', full, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tag', required=True)
    parser.add_argument('--ranges', nargs='+', required=True)
    args = parser.parse_args()
    verify(args.tag, args.ranges)

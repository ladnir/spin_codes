"""Exact combined coverage ledger; only hash-authenticated replays are accepted."""
from fractions import Fraction as F
from pathlib import Path

import certificate_search_core as core
import verify_certificate_search as verifier

base = core.base


def run():
    report_path = base.HERE/'generated/closure_t128_s19_m16_v1/report_0002.json'
    report = base.read(report_path)
    verifier.verify_report(report_path)
    spec = report['instance']
    best = {}
    dependencies = {report_path.relative_to(base.ROOT).as_posix(): base.sha(report_path),
                    Path(__file__).resolve().relative_to(base.ROOT).as_posix(): base.sha(Path(__file__))}
    for item in report['accepted_certificates']:
        path, replay = report_path.parent/item['certificate'], report_path.parent/item['replay']
        values = verifier.load_verified(path, replay, spec)
        core.merge_rows(best, [v for v in values if v['occupation'] in item['occupancies']], spec['rows'])
        dependencies[path.relative_to(base.ROOT).as_posix()] = base.sha(path)
        dependencies[replay.relative_to(base.ROOT).as_posix()] = base.sha(replay)
    accepted = []
    for path in sorted((base.HERE/'generated').glob('t128_s19_m16_q*_constant_v*.json')):
        saved = base.read(path)
        if saved.get('status') != 'COMPLETE_OCCUPANCY_CERTIFIED_50_BITS':
            continue
        replay_path = path.with_name(path.stem+'_replay.json')
        if not replay_path.exists():
            continue
        receipt = base.read(replay_path)
        core.authenticate(saved, base.ROOT)
        core.authenticate(receipt, base.ROOT)
        q = saved['occupation']
        core.require(saved['instance'] == spec and receipt['instance'] == spec and receipt['occupation'] == q,
                     'Cross-instance constant certificate')
        core.require(receipt['status'] == 'CONSTANT_SPLIT_512_BIT_REPLAY_PASSED' and
                     receipt['producer_sha256'] == base.sha(path), 'Missing constant-split replay')
        core.require(len(saved['upper_powers']) == q+1, 'Missing all-one-row case')
        bound = sum((F(2)**p for p in saved['upper_powers']), F(0))
        core.require(bound == base.decode(saved['union_upper']) and bound <= F(1, 1 << 50), 'Wrong case union')
        core.merge_rows(best, [dict(occupation=q, upper=base.encode(bound))], spec['rows'])
        accepted.append(dict(occupation=q, certificate=path.relative_to(base.ROOT).as_posix()))
        dependencies[path.relative_to(base.ROOT).as_posix()] = base.sha(path)
        dependencies[replay_path.relative_to(base.ROOT).as_posix()] = base.sha(replay_path)
    coverage = core.summarize(best, spec['rows'], 40)
    output = base.HERE/'generated/t128_s19_m16_combined_constant_coverage_v1.json'
    result = dict(status='EXACT_COMBINED_COVERAGE_LEDGER', instance=spec, coverage=coverage,
        constant_certificates=accepted, source_sha256=dependencies,
        numerical_replays_rerun_by_this_ledger=False)
    if output.exists():
        core.require(base.read(output) == result, 'Use a new ledger version when coverage changes')
    else:
        base.write_new(output, result)
    print(coverage['status'], coverage['covered_intervals'], 'gaps', coverage['uncovered_intervals'])
    print('covered contribution margin', coverage['covered_margin_bits'], 'full coverage', coverage['complete_coverage'])


if __name__ == '__main__':
    run()

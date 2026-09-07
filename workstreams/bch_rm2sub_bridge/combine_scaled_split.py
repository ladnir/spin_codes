"""Exact full/partial union of authenticated sparse and selected-split receipts."""
import argparse
from pathlib import Path

import search_scaled_split as split
import verify_certificate_search as verifier

core, base = split.core, split.base


def run(output, sparse_report, grids, bits=40):
    core.require(not output.exists() and grids, 'Existing output or no dense grids')
    verifier.verify_report(sparse_report)
    report = base.read(sparse_report)
    spec, best, dependencies = report['instance'], {}, {}
    def record(path):
        dependencies[path.resolve().relative_to(base.ROOT).as_posix()] = base.sha(path)
    record(sparse_report)
    for item in report['accepted_certificates']:
        cert, replay = sparse_report.parent/item['certificate'], sparse_report.parent/item['replay']
        values = verifier.load_verified(cert, replay, spec)
        core.merge_rows(best, [r for r in values if r['occupation'] in item['occupancies']], spec['rows'])
        record(cert)
        record(replay)
    for path in grids:
        saved, bounds, _ = split.load(path, require_replay=True)
        core.require(saved['instance'] == spec, 'Cross-instance dense grid')
        lo, hi = saved['interval']
        values = split.passing(bounds, lo, hi)
        core.merge_rows(best, [dict(occupation=q, upper=base.encode(v)) for q, v in values.items()], spec['rows'])
        record(path)
        record(path.with_name(path.stem+'_replay.json'))
    dependencies.update(split.density.dependencies.sources())
    result = dict(status='EXACT_SCALED_SPLIT_COVERAGE_LEDGER', instance=spec, coverage=core.summarize(best, spec['rows'], bits),
        per_occupancy_bounds=[dict(occupation=q, upper=base.encode(v)) for q, v in sorted(best.items())],
        source_sha256=dependencies, numerical_replays_rerun_by_this_ledger=False)
    base.write_new(output, result)
    coverage = result['coverage']
    print(coverage['status'], 'coverage', coverage['covered_intervals'], 'gaps', coverage['uncovered_intervals'],
        'covered margin', coverage['covered_margin_bits'], 'full', coverage['complete_coverage'], flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--sparse-report', type=Path, required=True)
    p.add_argument('--grid', type=Path, action='append', required=True)
    p.add_argument('--target-bits', type=int, default=40)
    a = p.parse_args()
    run(a.output, a.sparse_report, a.grid, a.target_bits)

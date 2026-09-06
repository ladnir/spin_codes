"""Audit finite-grid coverage, union provenance, reports, and snapshots.

Completion here means a complete diagnostic grid, not useful positive bounds
at every tuple and not an outward certificate. Default execution fails if
any native tuple lacks full evaluated occupation coverage.
"""
import argparse
import csv
import json
import math
import re
import sqlite3

import numpy as np

import read_grid_receipts
import run_complete_q1_grid as grid
import spectrum_events
import landscape_snapshot


def digest(path):
    return landscape_snapshot.sha(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--allow-partial',action='store_true')
    args=parser.parse_args();here=grid.HERE
    planned,observations,_=read_grid_receipts.snapshot()
    native={grid.key(r):r for r in planned if r['native']}
    if len(planned)!=3120 or len(native)!=3108 or set(observations)!=set(native):raise ValueError('finite grid scope mismatch')
    catalog=json.loads((here/'catalog.json').read_text());registered=set();dependencies={}
    for spec in catalog['activation_sources']:
        receipt=json.loads((here/spec['manifest']).read_text());path=here/spec['path']
        if digest(path)!=receipt['csv_sha256']:raise ValueError('CSV receipt mismatch')
        registered.add(spec['manifest'])
        for name,expected in receipt['source_sha256'].items():
            if name in dependencies and dependencies[name]!=expected:raise ValueError('incompatible receipt versions at one source path')
            dependencies[name]=expected
    for name,expected in dependencies.items():
        if digest(grid.pilot.ROOT/name)!=expected:raise ValueError(f'changed dependency: {name}')
    unregistered=[p.relative_to(here).as_posix() for p in here.glob('activation_occupation_grid*/t*/manifest.json')
                  if p.relative_to(here).as_posix() not in registered]
    events=spectrum_events.load_registered_events()
    with (here/'complete_grid_unions.csv').open(newline='') as handle:
        unions={grid.key(r):r for r in csv.DictReader(handle)}
    if set(unions)!=set(native):raise ValueError('union CSV omits a native tuple')
    details=json.loads((here/'complete_grid_unions.json').read_text())
    full={tuple(r['parameter_key']):r for r in details['full_union_sources']}
    if len(full)!=len(details['full_union_sources']):raise ValueError('duplicate full union')
    db=sqlite3.connect(f'{(here/"spin_landscape.sqlite3").as_uri()}?mode=ro',uri=True);db.row_factory=sqlite3.Row
    try:
        if db.execute('PRAGMA quick_check').fetchone()[0]!='ok':raise ValueError('database check failed')
        if db.execute("SELECT value FROM schema_metadata WHERE key='catalog_sha256'").fetchone()[0]!=digest(here/'catalog.json'):
            raise ValueError('database catalog provenance mismatch')
        if db.execute('SELECT COUNT(*) FROM certified_results').fetchone()[0]:raise ValueError('unexpected current outward-certificate claim')
        checked_sources=0
        for key,record in full.items():
            if key not in native:raise ValueError('full union outside finite grid')
            candidate=native[key];length=(1<<candidate['message_exponent'])//candidate['dimension']
            family='random' if candidate['outer_model']=='random-ensemble-average' else 'rm' if candidate['series'].startswith('RM(') else 'bch'
            if len(set(record['result_ids']))!=len(record['result_ids']):raise ValueError('duplicate union component')
            rows=[];event=record['setup_event_id']
            if event and (event not in events or events[event]['setup_failure_bits']!=record['setup_failure_bits']):
                raise ValueError('full union setup budget mismatch')
            for result_id in record['result_ids']:
                row=db.execute('SELECT * FROM landscape WHERE result_id=?',(result_id,)).fetchone()
                if row is None or row['transfer_review_status']!='activation_aware' or not row['comparison_eligible']:
                    raise ValueError('full union selected an ineligible observation')
                if row['outer_family']!=family:raise ValueError('full union mixed outer families')
                if (row['block_bits'],row['dimension'],row['step_bits'],row['state_bits'],row['message_exponent'],row['map_tag'])!=(
                        candidate['block_bits'],candidate['dimension'],candidate['step_bits'],candidate['state_bits'],candidate['message_exponent'],observations[key]['map_tag']):
                    raise ValueError('full union mixed parameter geometries')
                conditional=row['setup_event_id']
                if conditional and (not event or not spectrum_events.implies(events[event],events[conditional])):
                    raise ValueError('unchecked conditional-event splice')
                rows.append(row);checked_sources+=1
            through=0
            for row in sorted(rows,key=lambda r:r['occupation_min']):
                if row['occupation_min']>through+1:raise ValueError('gap in a claimed full union')
                if not 1<=row['occupation_min']<=row['occupation_max']<=length:raise ValueError('invalid occupation component')
                through=max(through,row['occupation_max'])
            if through!=length:raise ValueError('full union stops before L')
            value=float(np.logaddexp.reduce([-r['margin_bits']*math.log(2) for r in rows]))
            if event:value=float(np.logaddexp(value,-record['setup_failure_bits']*math.log(2)))
            if not math.isclose(value,record['log_bound'],rel_tol=1e-12,abs_tol=1e-8):raise ValueError('union sum mismatch')
            if not math.isclose(-value/math.log(2),float(unions[key]['full_union_margin_bits']),rel_tol=1e-12,abs_tol=1e-8):
                raise ValueError('union CSV margin mismatch')
            expected='positive_diagnostic_union' if value<0 else 'evaluated_bound_too_loose'
            if unions[key]['full_status']!=expected:raise ValueError('full evidence label mismatch')
        row_count=db.execute('SELECT COUNT(*) FROM landscape').fetchone()[0]
    finally:db.close()
    for key,row in unions.items():
        if (key in full)!=(row['full_status']!='incomplete_occupation_coverage'):raise ValueError('full-coverage status mismatch')
    frontiers=json.loads((here/'tested_parameter_frontiers.json').read_text())
    if (frontiers['native_parameter_count']!=3108
            or frontiers['union_snapshot_sha256']!=digest(here/'complete_grid_unions.csv')
            or frontiers['csv_sha256']!=digest(here/'tested_parameter_frontiers.csv')):
        raise ValueError('stale cost frontiers')
    projections=json.loads((here/'exact_family_projection_summary.json').read_text())
    if projections['primary_exact_blocks']!=dict(bch=[8,32,64,128],rm=[8,32,128,512]):raise ValueError('projection anchors changed')
    for field,path in [('projection_csv_sha256','exact_family_q1_projections.csv'),
                       ('holdout_csv_sha256','exact_family_projection_holdouts.csv'),
                       ('training_csv_sha256','exact_family_projection_training.csv'),
                       ('producer_sha256','extrapolate_exact_families.py')]:
        if projections[field]!=digest(here/path):raise ValueError('projection provenance mismatch')
    snapshot=json.loads((here/'landscape_snapshot.json').read_text())
    if snapshot['catalog_sha256']!=digest(here/'catalog.json') or snapshot['row_count']!=row_count:raise ValueError('stale compressed snapshot')
    for record in snapshot['files']:
        if digest(here/record['file'])!=record['plain_sha256'] or digest(here/record['compressed_file'])!=record['compressed_sha256']:
            raise ValueError('snapshot content mismatch')
    tests=(here/'grid_validation.log').read_text();match=re.search(r'Ran (\d+) tests',tests)
    if not match or not tests.rstrip().endswith('OK'):raise ValueError('successful test log missing')
    complete=len(full)==len(native) and not unregistered
    audit=dict(status='COMPLETE_DIAGNOSTIC_GRID' if complete else 'INCOMPLETE_DIAGNOSTIC_GRID',
        native_parameters=len(native),full_evaluated=len(full),full_positive=sum(r['full_status']=='positive_diagnostic_union' for r in unions.values()),
        full_certificates=0,registered_dependencies_checked=len(dependencies),union_components_checked=checked_sources,
        unregistered_completed_batches=unregistered,database_rows=row_count,tests_passed=int(match.group(1)),
        test_log_sha256=digest(here/'grid_validation.log'),
        note='A complete grid includes explicitly weak bounds. Only positive full diagnostics can guide distance selection; none is an outward certificate.')
    (here/'complete_landscape_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))
    if not complete and not args.allow_partial:raise SystemExit('Full finite-grid coverage is incomplete.')


if __name__=='__main__':main()

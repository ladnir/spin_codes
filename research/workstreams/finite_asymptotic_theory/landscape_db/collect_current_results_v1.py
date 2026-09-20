"""Collect current numerical tables without running or changing any producer.

Run the BCH grid audit first. --check verifies that the committed numerical
appendix still matches the local receipts. Generated research inputs stay local.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE/'CURRENT_RESULTS_TABLES.md'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(name):
    return json.loads((HERE/name).read_text(encoding='utf-8'))


def authenticate(data):
    for name, digest in data['source_sha256'].items():
        if sha(ROOT/Path(name)) != digest:
            raise ValueError('Changed source: '+name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    audit = read_json('bch_evidence_grid_audit_v1.json'); authenticate(audit)
    if audit['status'] != 'VERIFIED_COMPLETE_BCH_EVIDENCE_GRID' or audit['geometries'] != 130:
        raise ValueError('Expected the completed 130-point BCH audit')
    receipt = read_json('bch_engineering_evidence_v8.json')
    surfaces = read_json('engineering_surfaces.json'); authenticate(surfaces)
    for name, data in [('bch_engineering_evidence_v8.csv', receipt), ('engineering_surfaces.csv', surfaces)]:
        if sha(HERE/name) != data['csv_sha256']:
            raise ValueError('Changed numerical CSV: '+name)
    with (HERE/'bch_engineering_evidence_v8.csv').open(newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    keys = ('block_bits', 'step_bits', 'state_bits', 'message_exponent')
    indexed = {tuple(int(r[k]) for k in keys): r for r in rows}
    if len(rows) != 130 or len(indexed) != 130:
        raise ValueError('Duplicate or missing BCH geometry')
    if {key[0] for key in indexed} != {64,128}:
        raise ValueError('The refined BCH tables are restricted to blocks 64 and 128')
    used = set()

    def cell(b, t, s, e):
        key = (b, t, s, e); row = indexed[key]; used.add(key)
        if row['evidence'] == 'FIRST_MOMENT_OBSTRUCTION':
            return f"O ({float(row['first_moment_lower_bits']):.6g})"
        margin = float(row['full_margin_bits'])
        if margin <= 0:
            return f'W ({margin:.6g})'
        return f"{margin:.6f} ({float(row['full_aggregation_penalty_bits']):.6g})"

    lines = ['# Current numerical tables', '',
             'Collected from existing receipts; no new numerical search is performed.',
             'Read [the engineering overview](CURRENT_ENGINEERING_RESULTS.md) for interpretation.', '',
             '## Complete BCH-64/128 grid', '',
             'Positive entries are full margin in bits, followed by the full-bound loss relative to Q1 in parentheses.',
             'W gives a weak full upper bound and its negative margin. O gives a positive lower exponent a,',
             'meaning the bad-word first moment is at least 2^a; it does not assert actual code failure.',
             'Every original geometry appears in these tables.', '',
             f"Audit: {audit['useful_full_bounds']} useful full bounds, {audit['weak_full_bounds']} weak full bounds, "
             f"{audit['first_moment_obstructions']} first-moment obstructions; no missing geometries.", '',
             '### Message size at T=64, S=20', '',
             '| log2 K | BCH-64: margin (loss) | BCH-128: margin (loss) |',
             '| ---: | ---: | ---: |']
    for e in range(12, 27):
        lines.append(f'| {e} | {cell(64,64,20,e)} | {cell(128,64,20,e)} |')
    lines += ['', '### Epoch/state sweep at K=2^20', '',
              '| T | S | BCH-64: margin (loss) | BCH-128: margin (loss) |',
              '| ---: | ---: | ---: | ---: |']
    for t, minimum in ((64,7),(128,8),(256,9)):
        for s in range(minimum,21):
            lines.append(f'| {t} | {s} | {cell(64,t,s,20)} | {cell(128,t,s,20)} |')
    lines += ['', '### Additional message/state slices at T=64', '',
              'The exponent-20 slice is included above.', '',
              '| BCH block | S | log2 K=16 | log2 K=18 | log2 K=22 | log2 K=24 |',
              '| ---: | ---: | ---: | ---: | ---: | ---: |']
    for b in (64,128):
        for s in (10,12,16):
            lines.append(f'| {b} | {s} | '+' | '.join(cell(b,64,s,e) for e in (16,18,22,24))+' |')
    if used != set(indexed):
        raise ValueError('The tables do not cover the original BCH grid')
    with (HERE/'engineering_surfaces.csv').open(newline='', encoding='utf-8') as f:
        family = [r for r in csv.DictReader(f)
                  if r['family'] != 'bch' or int(r['block_bits']) <= 128]
    lines += ['', '## Matched family Q1 references at T=64, S=20', '',
              'Every entry in this section is Q1-only, including the random setup charge where shown.',
              'A random mean is an ensemble reference. A caps60 row is conditional on one shared spectrum event;',
              'its parentheses give the margin after adding that event failure once. None is a full-tail claim.', '',
              '| Family | B | log2 K=16 | log2 K=18 | log2 K=20 |',
              '| --- | ---: | ---: | ---: | ---: |']
    for name in ('bch','rm','random_mean','random_caps60'):
        for b in sorted({int(r['block_bits']) for r in family if r['family']==name}):
            values=[]
            for e in (16,18,20):
                matches=[r for r in family if (r['family'],int(r['block_bits']),int(r['step_bits']),int(r['state_bits']),int(r['message_exponent']))==(name,b,64,20,e)]
                if len(matches)!=1:
                    raise ValueError('Missing or duplicate family reference')
                row=matches[0]; value=f"{float(row['margin_bits']):.6f}"
                if row['setup_charged_margin_bits']:
                    value+=f" ({float(row['setup_charged_margin_bits']):.6f})"
                values.append(value)
            lines.append(f'| {name} | {b} | '+' | '.join(values)+' |')
    lines += ['', '## Positive full diagnostics in the earlier broad ledger', '',
              'This ledger is a separate snapshot; it does not contain the newer BCH grid or the separate RM certificate.',
              'Rounded random margins of 60 do not prove a strict greater-than-60-bit bound.', '',
              '| Series | log2 K | T | S | Full margin |', '| --- | ---: | ---: | ---: | ---: |']
    with (HERE/'complete_grid_unions.csv').open(newline='',encoding='utf-8') as f:
        unions=list(csv.DictReader(f))
    positive=[r for r in unions if r['full_union_margin_bits'] and float(r['full_union_margin_bits'])>0]
    if any('bch' in r['series'].lower() and int(r['block_bits'])>128 for r in positive):
        raise ValueError('Review the broad ledger: larger BCH results are outside this collection')
    if len(positive)!=read_json('complete_grid_unions.json')['full_positive_diagnostics']:
        raise ValueError('Broad ledger count mismatch')
    for row in positive:
        lines.append('| '+' | '.join(row[k] for k in ('series','message_exponent','step_bits','state_bits','full_union_margin_bits'))+' |')
    lines += ['', '## Input fingerprints', '',
              'These identify the local snapshots used for this document, not a new proof replay.', '',
              '| Local input | SHA-256 |', '| --- | --- |']
    for name in ('bch_evidence_grid_audit_v1.json','bch_engineering_evidence_v8.csv',
                 'engineering_surfaces.json','engineering_surfaces.csv','complete_grid_unions.json','complete_grid_unions.csv'):
        lines.append(f'| {name} | {sha(HERE/name)} |')
    text='\n'.join(lines)+'\n'
    if args.check:
        if OUTPUT.read_text(encoding='utf-8') != text:
            raise ValueError('Numerical appendix differs from the current snapshots')
        print('Verified numerical appendix: all 130 BCH geometries and matched family references.')
    else:
        OUTPUT.write_text(text,encoding='utf-8',newline='\n')
        print(OUTPUT)


if __name__ == '__main__':
    main()

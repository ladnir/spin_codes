"""Authenticated sequential BCH-64/128 higher-occupation grid.

The first pass evaluates Q2..4 at every existing engineering geometry.
It never interprets sparse dominance as dominance over the full tail.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import occupation_refresh_v1 as transfer
import read_grid_receipts as receipts
import run_occupation_grid as source

HERE = Path(__file__).resolve().parent
ROOT = source.grid.pilot.ROOT
LN2 = math.log(2)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_inputs():
    manifest = json.loads((HERE/'engineering_surfaces.json').read_text())
    dependencies = dict(manifest['source_sha256'])
    for name, digest in dependencies.items():
        if sha(ROOT/name) != digest:
            raise ValueError(f'changed engineering dependency: {name}')
    if sha(HERE/'engineering_surfaces.csv') != manifest['csv_sha256']:
        raise ValueError('engineering CSV changed')
    with (HERE/'engineering_surfaces.csv').open(newline='') as handle:
        rows = [r for r in csv.DictReader(handle) if r['family'] == 'bch' and int(r['block_bits']) in (64, 128)]
    _, observations, _ = receipts.snapshot()
    counts, maps = {}, {}
    for row in observations.values():
        if 'BCH' not in row['series'] or int(row['block_bits']) not in (64, 128):
            continue
        block, t, s = (int(row[k]) for k in ('block_bits', 'step_bits', 'state_bits'))
        if block not in counts:
            counts[block] = source.exact_counts(row)
        if (t, s) not in maps:
            path = ROOT/row['map_source']
            maps[t, s] = json.loads(path.read_text())
            if sha(path) != dependencies[row['map_source']]:
                raise ValueError('map not bound by engineering receipt')
    for path in (Path(__file__), Path(transfer.__file__), Path(receipts.__file__), Path(source.__file__),
                 HERE/'activation_occupation.py', HERE/'activation_q1_refresh.py',
                 HERE/'balanced_occupation.py', HERE/'composition_occupation.py',
                 HERE/'composition_boxes.py', HERE/'typed_dense_boxes.py',
                 HERE/'engineering_surfaces.csv', HERE/'engineering_surfaces.json'):
        dependencies[path.relative_to(ROOT).as_posix()] = sha(path)
    return rows, counts, maps, dependencies


def evaluate(row, counts, config):
    block, t, s, exponent = (int(row[k]) for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
    length = (1 << exponent)//(block//2)
    cutoff = block*length//10
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    epochs = transfer.Epochs(t, s, ac, config['kernel_counts'])
    tilts = np.unique(np.r_[np.arange(-20, 33, 2)/4-math.log(length), np.arange(-20, 5)/4])
    shifts = np.arange(-2, 11)/2
    models = {q: transfer.SparseModel(counts, block, q, 4, 4) for q in (2, 3, 4)}
    best = {q: np.full(len(m.indices), np.inf) for q, m in models.items()}
    selected = {q: np.zeros((len(m.indices), 2)) for q, m in models.items()}
    for tilt in tilts:
        lam = math.exp(tilt)
        regions = transfer.region_logs(epochs.at(lam, 4), t, length, 4)
        for q, model in models.items():
            for shift in shifts:
                values = model.components(regions, cutoff, lam, shift)
                improve = values < best[q]
                best[q][improve] = values[improve]
                selected[q][improve] = (tilt, shift)
    detail = {}
    for q, model in models.items():
        terms = best[q]+model.log_mult
        dominant = int(np.argmax(terms))
        total = model.aggregate(best[q], length)
        detail[str(q)] = dict(log_upper=total, bands=model.bands,
                             compositions=model.indices.tolist(), component_log_upper=best[q].tolist(),
                             witnesses=selected[q].tolist(), dominant_composition=dominant,
                             dominant_shift=float(selected[q][dominant, 1]),
                             dominant_log_tilt=float(selected[q][dominant, 0]),
                             dominant_at_shift_edge=bool(selected[q][dominant, 1] in (shifts[0], shifts[-1])),
                             dominant_at_tilt_edge=bool(selected[q][dominant, 0] in (tilts[0], tilts[-1])))
    q1 = -float(row['margin_bits'])*LN2
    tail = float(np.logaddexp.reduce([d['log_upper'] for d in detail.values()]))
    delta = (float(np.logaddexp(q1, tail))-q1)/LN2
    summary = dict(block_bits=block, step_bits=t, state_bits=s, message_exponent=exponent,
                   q1_margin_bits=float(row['margin_bits']),
                   q2_margin_bits=-detail['2']['log_upper']/LN2,
                   q3_margin_bits=-detail['3']['log_upper']/LN2,
                   q4_margin_bits=-detail['4']['log_upper']/LN2,
                   q2_q4_margin_bits=-tail/LN2, q2_q4_to_q1_log2_ratio=(tail-q1)/LN2,
                   q1_q4_penalty_bits=delta, q1_q4_margin_bits=-float(np.logaddexp(q1, tail))/LN2,
                   full_dominance_status='UNRESOLVED_Q5_AND_HIGHER')
    return dict(summary=summary, occupations=detail)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int)
    args = parser.parse_args()
    rows, counts, maps, dependencies = load_inputs()
    directory = HERE/'bch_dominance_v1'
    directory.mkdir(exist_ok=True)
    fingerprint = hashlib.sha256(json.dumps(dependencies, sort_keys=True).encode()).hexdigest()
    # Establish high-state endpoints, then knees, then the rest of the surface.
    def order(r):
        t, s, e = (int(r[k]) for k in ('step_bits', 'state_bits', 'message_exponent'))
        priority = 0 if s == 20 and e in (12, 16, 20, 24, 26) else 1 if s in (12, 14, 16) else 2
        return priority, int(r['block_bits']), e, t, s
    rows.sort(key=order)
    if args.limit:
        rows = rows[:args.limit]
    lock = HERE/'occupation_grid.lock'
    results = []
    with lock.open('x') as handle:
        handle.write('Sequential BCH dominance producer is active.\n')
    try:
        for ordinal, row in enumerate(rows, 1):
            b, t, s, e = (int(row[k]) for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
            path = directory/f'b{b}_t{t}_s{s}_e{e}.json'
            if path.exists():
                result = json.loads(path.read_text())
                if result['input_fingerprint'] != fingerprint:
                    raise ValueError(f'stale checkpoint {path}')
            else:
                result = evaluate(row, counts[b], maps[t, s])
                result['input_fingerprint'] = fingerprint
                path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
            results.append(result['summary'])
            r = result['summary']
            print(f'{ordinal}/{len(rows)} B{b} t{t} s{s} e{e}: Q1={r["q1_margin_bits"]:.3f}; '
                  f'Q2..4={r["q2_q4_margin_bits"]:.3f}; penalty={r["q1_q4_penalty_bits"]:.6g}', flush=True)
    finally:
        lock.unlink()
    output = HERE/'bch_dominance_sparse.csv'
    with output.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader(); writer.writerows(results)
    receipt = dict(schema='bch-dominance-sparse-v1', status='BINARY64_DIAGNOSTIC',
                   rows=len(results), input_fingerprint=fingerprint, source_sha256=dependencies,
                   csv_sha256=sha(output), checkpoints={p.name: sha(p) for p in directory.glob('*.json')},
                   limitations=['Q2..4 only; no full-tail dominance inference.',
                                'The ratio compares upper-bound contributions, not true event probabilities.',
                                'No outward rounding or interpolation across untested geometry.'])
    (HERE/'bch_dominance_sparse.json').write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    main()

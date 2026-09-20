"""Recompute the largest pilot size on a superset of the original tilt grid."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import activation_q1 as q1
import run_activation_pilot as pilot


def counts_for(row):
    b, k = int(row['block_bits']), int(row['dimension'])
    if row['outer_model'] == 'random-ensemble-average':
        offset = math.log((1 << k) - 1) - math.log((1 << b) - 1)
        return {w: offset + math.log(math.comb(b, w)) for w in range(1, b + 1)}
    key, c = next((key, c) for key, c in pilot.CONSTITUENTS.items() if c.name + ' exact' == row['series'])
    if c.spectrum_path.exists():
        counts = pilot.load_spectrum(c)
    else:
        payload = json.loads((pilot.SMALL / 'spectra' / f'{key}_weight_counts.json').read_text())
        counts = {int(w): int(n) for w, n in payload['weight_counts'].items()}
    return {w: math.log(n) for w, n in counts.items() if w and n}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--message-exponent', type=int, default=24)
    args = parser.parse_args()
    source = args.source.resolve()
    manifest_path = source.parent / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if pilot.sha(source) != manifest['csv_sha256']:
        raise ValueError('source CSV changed')
    for path, digest in manifest['source_sha256'].items():
        if pilot.sha(pilot.ROOT / path) != digest:
            raise ValueError(f'source dependency changed: {path}')
    with source.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = [r for r in reader if int(r['message_exponent']) == args.message_exponent]
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    output_csv = out / 'q1.csv'
    if output_csv.exists():
        raise ValueError('refusing to overwrite screen receipt')
    tilts = np.arange(-160, 1, dtype=float) / 10
    if not set(manifest['log_surprisal_grid']).issubset(set(tilts)):
        raise ValueError('new grid must contain all old witnesses')
    cached_map = None
    cache = {}
    improvements = []
    with output_csv.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            if row['map_source'] != cached_map:
                cached_map = row['map_source']
                cache.clear()
                config = json.loads((pilot.ROOT / cached_map).read_text())
                print(f'refine t{config["step_bits"]} s{config["state_bits"]}', flush=True)
            t, s = config['step_bits'], config['state_bits']
            b, k = int(row['block_bits']), int(row['dimension'])
            L = int(row['outer_rows'])
            if (b, k) not in cache:
                spectrum = {w: n for w, n in enumerate(config['a_counts']) if w and n}
                lam = np.exp(tilts)
                moments = q1.coefficient_logs(*q1.region_logs(*q1.epoch_logs(t, s, spectrum, lam), L // t), b)
                values = np.minimum(0., moments + int(row['bad_weight']) * lam[:, None])
                witness = np.argmin(values, axis=0)
                cache[b, k] = values[witness, np.arange(b + 1)], witness
            best, witness = cache[b, k]
            counts = counts_for(row)
            weights = sorted(counts)
            terms = np.array([math.log(L) + counts[w] + best[w] for w in weights])
            margin = -float(np.logaddexp.reduce(terms)) / math.log(2)
            dominant = weights[int(np.argmax(terms))]
            improvement = margin - float(row['margin_bits'])
            if improvement < -1e-9 or not math.isfinite(margin):
                raise ArithmeticError('larger witness grid made the bound worse')
            improvements.append(improvement)
            row.update(margin_bits=margin, dominant_weight=dominant,
                       dominant_log_surprisal=float(tilts[witness[dominant]]),
                       dominant_witness_at_grid_edge=int(witness[dominant] in (0, len(tilts) - 1)),
                       notes='Activation-aware Q1; expanded grid [-16,0]; supersedes the e24 coarse grid for comparisons; no full-distance claim.')
            writer.writerow(row)
            handle.flush()
    sources = dict(manifest['source_sha256'])
    for path in (Path(__file__), source, manifest_path):
        sources[path.relative_to(pilot.ROOT).as_posix()] = pilot.sha(path)
    manifest.update(row_count=len(rows), log_surprisal_grid=list(tilts),
                    arguments={'message_exponent': args.message_exponent},
                    source_sha256=sources, csv_sha256=pilot.sha(output_csv),
                    maximum_margin_improvement_bits=max(improvements),
                    remaining_dominant_grid_edges=sum(int(r['dominant_witness_at_grid_edge']) for r in rows))
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'refined {len(rows)} rows; max improvement {max(improvements):.6f} bits', flush=True)


if __name__ == '__main__':
    main()

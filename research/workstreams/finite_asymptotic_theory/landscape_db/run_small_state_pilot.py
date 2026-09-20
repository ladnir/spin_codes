"""Extend the frozen pilot chains by prefix restriction; no new map sampling."""
import csv
import json
import math
from pathlib import Path

import numpy as np

import activation_q1 as q1
import run_activation_pilot as pilot
from refine_activation_grid import counts_for

HERE = Path(__file__).resolve().parent
OUT = HERE / 'activation_pilot_small_state_v1'


def main():
    parent = HERE / 'activation_pilot_v1'
    manifest = json.loads((parent / 'manifest.json').read_text())
    for path, digest in manifest['source_sha256'].items():
        if pilot.sha(pilot.ROOT / path) != digest:
            raise ValueError(f'changed parent dependency: {path}')
    if pilot.sha(parent / 'q1.csv') != manifest['csv_sha256']:
        raise ValueError('changed parent CSV')
    with (parent / 'q1.csv').open(newline='') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        templates = {}
        for row in reader:
            templates.setdefault(row['series'], row)
    if (OUT / 'q1.csv').exists():
        raise ValueError('write-once output already exists')
    (OUT / 'maps').mkdir(parents=True, exist_ok=True)
    sources = dict(manifest['source_sha256'])
    for path in (Path(__file__), HERE / 'refine_activation_grid.py', parent / 'q1.csv', parent / 'manifest.json'):
        sources[path.relative_to(pilot.ROOT).as_posix()] = pilot.sha(path)
    tilts = np.arange(-160, 1, dtype=float) / 10
    completed = 0
    with (OUT / 'q1.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for t in (64, 128, 256):
            parent_path = parent / 'maps' / f't{t}_s20.json'
            full = json.loads(parent_path.read_text())
            full_generators = [int(g, 16) for g in full['generator_words_hex']]
            for s in (12, 13, 15):
                generators = full_generators[:s]
                columns = [sum(((g >> j) & 1) << i for i, g in enumerate(generators)) for j in range(t)]
                if pilot.rank(generators) != s or len(set(columns)) != t or 0 in columns:
                    raise ArithmeticError('rank/column audit failed')
                if any((a & b).bit_count() & 1 for a in generators for b in generators):
                    raise ArithmeticError('BA audit failed')
                spectrum = pilot.enumerate_spectrum(generators, s, t)
                kernel = pilot.macwilliams_kernel(spectrum, s)
                if kernel[4] != pilot.kernel_weight_four_count(columns):
                    raise ArithmeticError('kernel audit failed')
                config = dict(step_bits=t, state_bits=s, chain_seed=full['chain_seed'],
                              generator_words_hex=[hex(g) for g in generators], a_counts=spectrum, kernel_counts=kernel,
                              checks=dict(full_rank=True, BA_zero=True, distinct_nonzero_columns=True,
                                          kernel_weight_four_independent=True, prefix_of_parent=True),
                              parent_source=parent_path.relative_to(pilot.ROOT).as_posix(), parent_sha256=pilot.sha(parent_path))
                map_path = OUT / 'maps' / f't{t}_s{s}.json'
                map_path.write_text(json.dumps(config, indent=2) + '\n')
                sources[map_path.relative_to(pilot.ROOT).as_posix()] = pilot.sha(map_path)
                nz = {w: n for w, n in enumerate(spectrum) if w and n}
                for e in (16, 18, 20):
                    cache = {}
                    for template in templates.values():
                        row = dict(template)
                        b, k = int(row['block_bits']), int(row['dimension'])
                        L = (1 << e) // k
                        if (1 << e) % k or L % t:
                            raise ValueError('non-native length')
                        if (b, k) not in cache:
                            lam = np.exp(tilts)
                            moments = q1.coefficient_logs(*q1.region_logs(*q1.epoch_logs(t, s, nz, lam), L // t), b)
                            values = np.minimum(0., moments + ((b * L) // 10) * lam[:, None])
                            witness = np.argmin(values, axis=0)
                            cache[b, k] = values[witness, np.arange(b + 1)], witness
                        best, witness = cache[b, k]
                        counts = counts_for(row)
                        weights = sorted(counts)
                        terms = np.array([math.log(L) + counts[w] + best[w] for w in weights])
                        dominant = weights[int(np.argmax(terms))]
                        margin = -float(np.logaddexp.reduce(terms)) / math.log(2)
                        if not math.isfinite(margin):
                            raise ArithmeticError('nonfinite margin')
                        row.update(message_exponent=e, message_bits=1 << e, output_bits=b * L, outer_rows=L,
                                   bad_weight=(b * L) // 10, step_bits=t, state_bits=s, epochs_per_region=L // t,
                                   margin_bits=margin, dominant_weight=dominant,
                                   dominant_log_surprisal=float(tilts[witness[dominant]]),
                                   dominant_witness_at_grid_edge=int(witness[dominant] in (0, len(tilts) - 1)),
                                   a_minimum_distance=min(nz), kernel_minimum_distance=next(w for w,n in enumerate(kernel) if w and n),
                                   kernel_weight_four=kernel[4], map_tag='nested-' + pilot.sha(map_path)[:16],
                                   map_source=map_path.relative_to(pilot.ROOT).as_posix(),
                                   notes='Activation-aware Q1; exact prefix of frozen s20 pilot map; no full-distance claim.')
                        writer.writerow(row)
                        completed += 1
                    handle.flush()
                    print(f't{t} s{s} e{e}: {completed} rows', flush=True)
    manifest.update(row_count=completed, arguments={'states': [12,13,15], 'message_exponents':[16,18,20]},
                    log_surprisal_grid=list(tilts), source_sha256=sources, csv_sha256=pilot.sha(OUT / 'q1.csv'))
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()

"""Fill the finite grid sequentially, reusing authenticated pilot observations.

Each map is a restart boundary. Only a CSV with a verified final manifest is
complete; a pending CSV can be recomputed after an interrupted run.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import activation_q1 as q1
import run_activation_pilot as pilot
from refine_activation_grid import counts_for

HERE = Path(__file__).resolve().parent
OUT = HERE / 'activation_grid_v1'
EXPONENTS = (16, 18, 20, 22, 24)
STEPS = (64, 128, 256)
RANDOM_BLOCKS = (8, 16, 32, 64, 128, 256, 512, 1024)
BASES = ('activation_pilot_v1', 'activation_pilot_scaling_v1',
         'activation_pilot_e24_refined_v1', 'activation_pilot_small_state_v1')


def verify(directory):
    receipt = json.loads((directory / 'manifest.json').read_text())
    if pilot.sha(directory / 'q1.csv') != receipt['csv_sha256']:
        raise ValueError(f'changed CSV: {directory}')
    for path, digest in receipt['source_sha256'].items():
        if pilot.sha(pilot.ROOT / path) != digest:
            raise ValueError(f'changed dependency: {path}')
    return receipt


def key(row):
    return (row['series'], int(row['step_bits']), int(row['state_bits']),
            int(row['message_exponent']))


def inventory():
    existing, maps, sources, templates = {}, {}, {}, {}
    fields = None
    for name in BASES:
        directory = HERE / name
        receipt = verify(directory)
        sources.update(receipt['source_sha256'])
        for path in (directory / 'manifest.json', directory / 'q1.csv'):
            sources[path.relative_to(pilot.ROOT).as_posix()] = pilot.sha(path)
        with (directory / 'q1.csv').open(newline='') as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            for row in reader:
                if name == 'activation_pilot_scaling_v1' and int(row['message_exponent']) == 24:
                    continue
                if key(row) in existing:
                    raise ValueError(f'duplicate preferred pilot tuple: {key(row)}')
                existing[key(row)] = row
                maps.setdefault((int(row['step_bits']), int(row['state_bits'])),
                                pilot.ROOT / row['map_source'])
                templates.setdefault(row['series'], row)
    models = [dict(series=c.name + ' exact', outer_model='fixed',
                   block_bits=c.block_bits, dimension=c.dimension)
              for c in pilot.CONSTITUENTS.values()]
    models += [dict(series=f'random full-rank [{b},{b//2}] reused',
                    outer_model='random-ensemble-average', block_bits=b, dimension=b//2)
               for b in RANDOM_BLOCKS]
    for path in (Path(__file__), HERE / 'refine_activation_grid.py'):
        sources[path.relative_to(pilot.ROOT).as_posix()] = pilot.sha(path)
    return existing, maps, sources, models, fields


def tuples(models):
    for t in STEPS:
        for s in range(t.bit_length(), 21):
            for e in EXPONENTS:
                for model in models:
                    row = dict(model, step_bits=t, state_bits=s, message_exponent=e)
                    L, remainder = divmod(1 << e, model['dimension'])
                    row['native'] = remainder == 0 and L % t == 0
                    yield row


def map_for(t, s, maps):
    if (t, s) in maps:
        return maps[t, s]
    parent_path = HERE / 'activation_pilot_v1' / 'maps' / f't{t}_s20.json'
    parent = json.loads(parent_path.read_text())
    generators = [int(g, 16) for g in parent['generator_words_hex'][:s]]
    columns = [sum(((g >> j) & 1) << i for i, g in enumerate(generators)) for j in range(t)]
    if pilot.rank(generators) != s or len(set(columns)) != t or 0 in columns:
        raise ArithmeticError('map rank/columns failed')
    if any((a & b).bit_count() & 1 for a in generators for b in generators):
        raise ArithmeticError('BA failed')
    spectrum = pilot.enumerate_spectrum(generators, s, t)
    kernel = pilot.macwilliams_kernel(spectrum, s)
    if kernel[4] != pilot.kernel_weight_four_count(columns):
        raise ArithmeticError('independent kernel weight-four audit failed')
    payload = dict(step_bits=t, state_bits=s, chain_seed=parent['chain_seed'],
                   generator_words_hex=[hex(g) for g in generators],
                   a_counts=spectrum, kernel_counts=kernel,
                   checks=dict(full_rank=True, BA_zero=True, distinct_nonzero_columns=True,
                               kernel_weight_four_independent=True, prefix_of_parent=True),
                   parent_source=parent_path.relative_to(pilot.ROOT).as_posix(),
                   parent_sha256=pilot.sha(parent_path))
    path = OUT / 'maps' / f't{t}_s{s}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(payload, indent=2) + '\n'
    if path.exists() and path.read_text() != encoded:
        raise ValueError(f'changed map: {path}')
    if not path.exists():
        path.write_text(encoded)
    return path


def run_map(t, s, rows, maps, sources, fields):
    directory = OUT / f't{t}_s{s}'
    if (directory / 'manifest.json').exists():
        receipt = verify(directory)
        with (directory / 'q1.csv').open(newline='') as handle:
            actual = [key(r) for r in csv.DictReader(handle)]
        if sorted(actual) != sorted(key(r) for r in rows):
            raise ValueError('completed batch does not match the grid')
        return receipt['row_count']
    if (directory / 'q1.csv').exists():
        raise ValueError(f'CSV without final receipt requires inspection: {directory}')
    path = map_for(t, s, maps)
    config = json.loads(path.read_text())
    nz = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    kernel = config['kernel_counts']
    tilts = np.arange(-160, 1, dtype=float) / 10
    lam = np.exp(tilts)
    epoch = q1.epoch_logs(t, s, nz, lam)
    directory.mkdir(parents=True, exist_ok=True)
    pending = directory / 'q1.pending.csv'
    count = 0
    # Cache one geometry at a time. Spectra of equal length share the transfer.
    rows = sorted(rows, key=lambda r: (r['message_exponent'], r['block_bits'], r['series']))
    previous = None
    with pending.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            e, b, d = row['message_exponent'], row['block_bits'], row['dimension']
            L = (1 << e) // d
            if previous != (e, b, d):
                moments = q1.coefficient_logs(*q1.region_logs(*epoch, L // t), b)
                values = np.minimum(0., moments + ((b * L) // 10) * lam[:, None])
                witness = np.argmin(values, axis=0)
                best = values[witness, np.arange(b + 1)]
                previous = (e, b, d)
            counts = counts_for(row)
            weights = sorted(counts)
            terms = np.array([math.log(L) + counts[w] + best[w] for w in weights])
            dominant = weights[int(np.argmax(terms))]
            margin = -float(np.logaddexp.reduce(terms)) / math.log(2)
            if not math.isfinite(margin):
                raise ArithmeticError('nonfinite screen')
            result = {k: v for k, v in row.items() if k != 'native'}
            result.update(message_bits=1 << e, output_bits=b * L, outer_rows=L,
                          bad_weight=(b * L) // 10, epochs_per_region=L // t,
                          margin_bits=margin, dominant_weight=dominant,
                          dominant_log_surprisal=float(tilts[witness[dominant]]),
                          dominant_witness_at_grid_edge=int(witness[dominant] in (0, len(tilts)-1)),
                          a_minimum_distance=min(nz),
                          kernel_minimum_distance=next(w for w,n in enumerate(kernel) if w and n),
                          kernel_weight_four=kernel[4], map_tag='nested-' + pilot.sha(path)[:16],
                          map_source=path.relative_to(pilot.ROOT).as_posix(),
                          notes='Complete finite grid Q1; exact frozen-chain prefix; nearest binary64; higher occupations separate.')
            writer.writerow(result)
            count += 1
            handle.flush()
    dependencies = dict(sources)
    dependencies[path.relative_to(pilot.ROOT).as_posix()] = pilot.sha(path)
    receipt = dict(schema='activation-aware-q1-grid-v1', status='BINARY64_DIAGNOSTIC',
                   occupation=1, transfer_review_status='activation_aware', row_count=count,
                   distance_target='1/10', bad_event='weight <= floor(N/10)',
                   log_surprisal_grid=list(tilts), source_sha256=dependencies,
                   csv_sha256=pilot.sha(pending), numpy_version=np.__version__,
                   limitations=['Q1 only; higher occupations are separate obligations.',
                                'Random rows are ensemble expectations for one reused full-rank code.',
                                'One unselected nested chain per t; no outward certificate.'])
    # Publish the receipt last. Incomplete files are never included in the index.
    (directory / 'manifest.pending.json').write_text(json.dumps(receipt, indent=2) + '\n')
    pending.replace(directory / 'q1.csv')
    (directory / 'manifest.pending.json').replace(directory / 'manifest.json')
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan-only', action='store_true')
    parser.add_argument('--limit-maps', type=int)
    args = parser.parse_args()
    existing, maps, sources, models, fields = inventory()
    planned = list(tuples(models))
    native = [r for r in planned if r['native']]
    if not set(existing).issubset({key(r) for r in native}):
        raise ValueError('preferred pilot lies outside the finite grid')
    OUT.mkdir(exist_ok=True)
    summary = dict(grid_version=1, message_exponents=list(EXPONENTS), steps=list(STEPS),
                   state_rule='log2(t)+1 through 20 inclusive', random_blocks=list(RANDOM_BLOCKS),
                   models=models, candidate_count=len(planned), native_count=len(native),
                   reused_q1_count=len(existing), new_q1_count=len(native)-len(existing),
                   inapplicable=[dict(r, reason='outer row count is not divisible by epoch size')
                                 for r in planned if not r['native']])
    encoded = json.dumps(summary, indent=2) + '\n'
    plan = OUT / 'grid_plan.json'
    if plan.exists() and plan.read_text() != encoded:
        raise ValueError('grid plan changed; use a new grid version')
    if not plan.exists():
        plan.write_text(encoded)
    print(json.dumps({k:summary[k] for k in ('candidate_count','native_count','reused_q1_count','new_q1_count')}), flush=True)
    if args.plan_only:
        return
    # Exclusive creation prevents two instances of this numerical producer.
    lock = OUT / 'run.lock'
    with lock.open('x') as handle:
        handle.write('Sequential grid runner is active. Remove only after verifying the process has stopped.\n')
    try:
        completed = 0
        for t in STEPS:
            for s in range(t.bit_length(), 21):
                rows = [r for r in native if (r['step_bits'],r['state_bits']) == (t,s) and key(r) not in existing]
                if not rows:
                    continue
                was_complete = (OUT / f't{t}_s{s}' / 'manifest.json').exists()
                count = run_map(t, s, rows, maps, sources, fields)
                print(f't{t} s{s}: {count} new-grid Q1 rows verified', flush=True)
                if not was_complete:
                    completed += 1
                if args.limit_maps is not None and completed >= args.limit_maps:
                    return
    finally:
        lock.unlink()


if __name__ == '__main__':
    main()

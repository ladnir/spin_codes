"""Sequential, activation-aware Q1 pilot on fixed nested RM2Sub maps."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

import activation_q1

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent.parent
SMALL = HERE.parent / 'small_k_replay'
sys.path.insert(0, str(SMALL))
sys.path.insert(0, str(HERE.parent))
from evaluate_exact_spectra_q1_phase import CONSTITUENTS, load_spectrum
from generate_rm2sub_calibration_constituent import (
    coordinate_columns, enumerate_spectrum, generator_words, kernel_weight_four_count,
    macwilliams_kernel, rank, sample_independent_masks)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_maps(directory, states, steps, seed):
    import random
    maps = []
    for t in steps:
        m = t.bit_length() - 1
        if 1 << m != t or m < 5 or min(states) < m + 1:
            raise ValueError('invalid RM2Sub epoch or state dimensions')
        monomials = [(i, j) for i in range(m) for j in range(i + 1, m)]
        if max(states) > m + 1 + len(monomials):
            raise ValueError('state exceeds RM(2,m) capacity')
        # Include the matched-persistence p=26 member at t=128.
        selected_states = sorted(set(states) | ({19} if t == 128 else set()))
        quadratics = sample_independent_masks(random.Random(seed + t), len(monomials), max(selected_states) - m - 1)
        columns = coordinate_columns(m, quadratics, monomials)
        full_generators = generator_words(columns, max(selected_states))
        for s in selected_states:
            generators = full_generators[:s]
            cols = [c & ((1 << s) - 1) for c in columns]
            if rank(generators) != s or len(set(cols)) != t or 0 in cols:
                raise ArithmeticError('rank or distinct-nonzero-column check failed')
            if any((a & b).bit_count() & 1 for a in generators for b in generators):
                raise ArithmeticError('BA is not zero')
            spectrum = enumerate_spectrum(generators, s, t)
            kernel = macwilliams_kernel(spectrum, s)
            if kernel[4] != kernel_weight_four_count(cols):
                raise ArithmeticError('independent kernel weight-four count failed')
            payload = dict(step_bits=t, state_bits=s, chain_seed=seed + t,
                           generator_words_hex=[hex(g) for g in generators],
                           a_counts=spectrum, kernel_counts=kernel,
                           checks=dict(full_rank=True, BA_zero=True, distinct_nonzero_columns=True,
                                       kernel_weight_four_independent=True))
            path = directory / f't{t}_s{s}.json'
            encoded = json.dumps(payload, indent=2) + '\n'
            if path.exists() and path.read_text() != encoded:
                raise ValueError(f'refusing to overwrite different map: {path}')
            path.write_text(encoded, encoding='utf-8')
            maps.append((path, payload))
    return maps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--message-exponents', nargs='+', type=int, default=[16, 20])
    parser.add_argument('--states', nargs='+', type=int, default=[14, 16, 18, 20])
    parser.add_argument('--steps', nargs='+', type=int, default=[64, 128, 256])
    parser.add_argument('--constituents', nargs='+', choices=list(CONSTITUENTS), default=list(CONSTITUENTS))
    parser.add_argument('--random-blocks', nargs='*', type=int, default=[64, 128, 512])
    parser.add_argument('--seed', type=int, default=3390173185)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    map_dir = out / 'maps'
    map_dir.mkdir(exist_ok=True)
    maps = make_maps(map_dir, args.states, args.steps, args.seed)
    families = []
    dependencies = [Path(__file__), Path(activation_q1.__file__),
                    SMALL / 'generate_rm2sub_calibration_constituent.py',
                    SMALL / 'evaluate_exact_spectra_q1_phase.py']
    for key in args.constituents:
        c = CONSTITUENTS[key]
        source = c.spectrum_path
        if source.exists():
            counts = load_spectrum(c)
        else:
            # The curated checkout retains normalized exact spectra instead of
            # the older scripts-directory originals for BCH128 and RM49.
            source = SMALL / 'spectra' / f'{key}_weight_counts.json'
            payload = json.loads(source.read_text(encoding='utf-8'))
            if (payload['length'], payload['dimension']) != (c.block_bits, c.dimension):
                raise ValueError('normalized spectrum parameter mismatch')
            counts = {int(w): int(n) for w, n in payload['weight_counts'].items()}
            if counts.get(0) != 1 or sum(counts.values()) != 1 << c.dimension:
                raise ValueError('normalized spectrum mass mismatch')
            if min(w for w, n in counts.items() if w and n) != c.minimum_distance:
                raise ValueError('normalized spectrum minimum distance mismatch')
        families.append((c.name + ' exact', key, c.block_bits, c.dimension,
                         {w: math.log(n) for w, n in counts.items() if w and n}, 'fixed'))
        dependencies.append(source)
    for b in args.random_blocks:
        if b < 2 or b % 2:
            raise ValueError('random outer must have positive even length')
        k = b // 2
        offset = math.log((1 << k) - 1) - math.log((1 << b) - 1)
        families.append((f'random full-rank [{b},{k}] reused', f'random{b}', b, k,
                         {w: offset + math.log(math.comb(b, w)) for w in range(1, b + 1)}, 'random-ensemble-average'))
    tilts = np.arange(-120, 1, dtype=float) / 10
    csv_path = out / 'q1.csv'
    if csv_path.exists():
        raise ValueError('choose a new output directory; screen receipts are write-once')
    rows = []
    fields = ['series', 'outer_model', 'block_bits', 'dimension', 'message_exponent',
              'message_bits', 'output_bits', 'outer_rows', 'bad_weight', 'step_bits',
              'state_bits', 'epochs_per_region', 'margin_bits', 'dominant_weight',
              'dominant_log_surprisal', 'dominant_witness_at_grid_edge',
              'a_minimum_distance', 'kernel_minimum_distance', 'kernel_weight_four',
              'map_tag', 'map_source', 'notes']
    with csv_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for map_path, config in maps:
            dependencies.append(map_path)
            t, s = config['step_bits'], config['state_bits']
            spectrum = {w: n for w, n in enumerate(config['a_counts']) if w and n}
            kd = next(w for w, n in enumerate(config['kernel_counts']) if w and n)
            for exponent in args.message_exponents:
                # Cache the transfer by geometry: BCH/RM/random at the same size
                # share the complete coefficient calculation. Only spectra differ.
                cache = {}
                for label, key, b, k, counts, model in families:
                    geometry = (b, k)
                    if geometry not in cache:
                        L = (1 << exponent) // k
                        if (1 << exponent) % k or L % t:
                            raise ValueError(f'non-native pilot tuple: {key}, e{exponent}, t{t}')
                        lam = np.exp(tilts)
                        moments = activation_q1.coefficient_logs(*activation_q1.region_logs(
                            *activation_q1.epoch_logs(t, s, spectrum, lam), L // t), b)
                        values = np.minimum(0., moments + ((b * L) // 10) * lam[:, None])
                        witness = np.argmin(values, axis=0)
                        best = values[witness, np.arange(b + 1)]
                        cache[geometry] = (L, best, witness)
                    L, best, witness = cache[geometry]
                    weights = sorted(counts)
                    terms = np.array([math.log(L) + counts[w] + best[w] for w in weights])
                    dominant = weights[int(np.argmax(terms))]
                    margin = -float(np.logaddexp.reduce(terms)) / math.log(2)
                    if not math.isfinite(margin):
                        raise ArithmeticError('nonfinite screen')
                    row = dict(series=label, outer_model=model, block_bits=b, dimension=k,
                               message_exponent=exponent, message_bits=1 << exponent,
                               output_bits=b * L, outer_rows=L, bad_weight=(b * L) // 10,
                               step_bits=t, state_bits=s, epochs_per_region=L // t,
                               margin_bits=margin, dominant_weight=dominant,
                               dominant_log_surprisal=float(tilts[witness[dominant]]),
                               dominant_witness_at_grid_edge=int(witness[dominant] in (0, len(tilts) - 1)),
                               a_minimum_distance=min(spectrum), kernel_minimum_distance=kd,
                               kernel_weight_four=config['kernel_counts'][4],
                               map_tag='nested-' + sha(map_path)[:16],
                               map_source=map_path.relative_to(ROOT).as_posix(),
                               notes='Activation-aware Q1 only; unselected nested chain; nearest binary64; no full-distance claim.')
                    writer.writerow(row)
                    rows.append(row)
                handle.flush()
                print(f't{t} s{s} e{exponent}: {len(rows)} rows complete', flush=True)
    manifest = dict(schema='activation-aware-q1-pilot-v1', status='BINARY64_DIAGNOSTIC',
                    transfer_review_status='activation_aware', occupation=1, row_count=len(rows),
                    distance_target='1/10', bad_event='weight <= floor(N/10)',
                    log_surprisal_grid=list(tilts), arguments={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                    source_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(set(dependencies))},
                    csv_sha256=sha(csv_path), numpy_version=np.__version__,
                    limitations=['Q1 only; all higher occupations remain unscreened in this pilot.',
                                 'No outward arithmetic or new distance certificate.',
                                 'One unselected chain per t; cross-t comparisons include map variation.',
                                 'Random entries average over one full-rank outer reused everywhere.',
                                 'No performance measurements; number of epochs is only an operation-count input.'])
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(f'completed {len(rows)} rows: {csv_path}', flush=True)


if __name__ == '__main__':
    main()

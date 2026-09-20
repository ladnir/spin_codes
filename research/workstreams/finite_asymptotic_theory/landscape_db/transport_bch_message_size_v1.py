"""Reuse fixed-Q sparse partitions at a different K, then search the dense tail.

B, T and S stay fixed. Sparse type partitions depend on Q rather than L;
the target re-evaluates all coefficients, row counts and bounds. Witness
surprisals scale by L_source/L_target. This is a witness proposal, not an
assumption about transfer scaling. Dense integer partitions are rebuilt.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np

import study_bch_dominance_v1 as study
import transport_bch_full_cover_v1 as transport

HERE = Path(__file__).resolve().parent


def target_sparse(source, exponent, source_length, counts, block, t, s, config):
    length = (1 << exponent)//(block//2)
    data = copy.deepcopy(source)
    maximum = min(length, source['occupation_max'])
    data['occupation_max'] = maximum
    data['occupations'] = [r for r in data['occupations'] if r['occupation'] <= maximum]
    shift = math.log(source_length/length)
    for row in data['occupations']:
        row['witness_log_tilt'] += shift
        if row['composition_cover'] is not None:
            for box in row['composition_cover']['boxes']:
                box['witness']['log_tilt'] += shift
    data = transport.sparse_cover(data, counts, block, t, s, length, config)
    data['arguments'].update(exponent=exponent, maximum=maximum)
    total = float(np.logaddexp.reduce([r['log_upper'] for r in data['occupations']]))
    data.update(log_upper=total, margin_bits=-total/math.log(2))
    data.pop('input_fingerprint', None)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--exponent', type=int, action='append', required=True)
    parser.add_argument('--maximum-refinements', type=int, default=300)
    args = parser.parse_args()
    if args.maximum_refinements < 0:
        parser.error('maximum-refinements must be nonnegative')
    reference_path = args.reference.resolve()
    reference = json.loads(reference_path.read_text())
    if reference['status'] != 'VERIFIED_BINARY64_FULL_REFERENCE':
        raise ValueError('A verified source reference is required')
    transport.authenticate(reference)
    b, t, s, source_exponent = (reference[k] for k in
                    ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
    source_length = (1 << source_exponent)//(b//2)
    intervals = transport.select_covers(reference)
    sparse = [(p, d) for p, d in intervals if d['status'] == 'BINARY64_COMPLETE_SPARSE_INTERVAL']
    if len(sparse) != 1 or sparse[0][1]['occupation_min'] != 5:
        raise ValueError('Expected one source sparse interval starting at Q5')
    source_path, source = sparse[0]
    _, spectra, maps, dependencies = study.load_inputs()
    dependencies.update(reference['source_sha256'])
    for path in (Path(__file__), Path(transport.__file__), Path(transport.fast.__file__), reference_path, source_path):
        dependencies[str(path.resolve())] = study.sha(path)
    for exponent in args.exponent:
        if exponent <= 0 or exponent == source_exponent:
            parser.error('Expected a positive target exponent distinct from the source')
        length = (1 << exponent)//(b//2)
        if length < 5 or length % t:
            parser.error('Expected L>=5 and complete epochs')
        if not (HERE/f'bch_dominance_v1/b{b}_t{t}_s{s}_e{exponent}.json').exists():
            parser.error('Target needs an original Q2..4 checkpoint')
    for exponent in dict.fromkeys(args.exponent):
        length = (1 << exponent)//(b//2)
        tag = f'b{b}_t{t}_s{s}_e{exponent}'
        full_path = HERE/f'bch_full_reference_{tag}.json'
        if full_path.exists():
            existing = json.loads(full_path.read_text())
            if (existing['status'] != 'VERIFIED_BINARY64_FULL_REFERENCE' or
                    tuple(existing[k] for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent')) != (b, t, s, exponent)):
                raise ValueError('Existing reference is incompatible')
            transport.authenticate(existing)
            print(f'{tag}: existing verified reference retained', flush=True)
            continue
        maximum = min(length, source['occupation_max']); minimum = maximum+1
        directory = HERE/f'bch_sparse_k_transport_v1_{tag}_q5_{maximum}'
        directory.mkdir(exist_ok=True)
        path = directory/'cover.json'
        origin = dict(reference=str(reference_path), reference_sha256=study.sha(reference_path),
                      source_interval=str(source_path), source_interval_sha256=study.sha(source_path),
                      source_exponent=source_exponent, target_exponent=exponent,
                      witness_log_tilt_shift=math.log(source_length/length))
        if path.exists():
            current = json.loads(path.read_text())
            if current.get('message_size_transport') != origin or current['source_sha256'] != dependencies:
                raise ValueError('Existing sparse transport has different inputs')
            transport.authenticate(current)
        else:
            print(f'{tag}: sparse transport', flush=True)
            current = target_sparse(source, exponent, source_length, spectra[b], b, t, s, maps[t, s])
            current.update(source_sha256=dependencies, message_size_transport=origin,
                           limitations=['Fixed-Q partitions; every target coefficient and bound is recomputed.',
                                        'Only the displayed occupation interval is covered.', 'No outward arithmetic.'])
            temporary = path.with_suffix('.json.tmp')
            temporary.write_text(json.dumps(current, indent=2)+'\n'); temporary.replace(path)
        common = ['--block', str(b), '--step', str(t), '--state', str(s), '--exponent', str(exponent)]

        def run(script, extra, stage, with_geometry=True):
            print(f'{tag}: {stage}', flush=True)
            with (HERE/f'bch_k_transport_v1_{tag}_{stage}.log').open('w') as output:
                subprocess.run([sys.executable, '-u', str(HERE/script), *(common if with_geometry else []), *extra],
                               cwd=HERE, stdout=output, stderr=subprocess.STDOUT, check=True)

        dense_path = HERE/f'bch_dense_v11_{tag}_q{minimum}/cover.json'
        if minimum <= length:
            seed = HERE/f'bch_dense_seed_bulk_{tag}_q{minimum}.json'
            if not seed.exists():
                run('seed_bch_dense_v3.py', ['--minimum', str(minimum), '--nodes', '127' if b == 128 else '63'], 'seed')
            if not dense_path.exists():
                run('close_bch_dense_v11.py', ['--minimum', str(minimum), '--target-bits', '55' if b == 128 else '30',
                                             '--maximum-refinements', str(args.maximum_refinements)], 'dense')
        refined = HERE/f'bch_dominance_v2/{tag}.json'
        if not refined.exists():
            run('refine_bch_sparse_v2.py', ['--checkpoint', f'{tag}.json'], 'refine', with_geometry=False)
        extra = ['--dense-minimum', str(minimum), '--sparse-cover', str(path), '--sparse-checkpoint', str(refined)]
        if minimum <= length:
            extra += ['--dense-cover', str(dense_path)]
        run('verify_bch_full_reference_v5.py', extra, 'verify')
        final = json.loads(full_path.read_text())
        print(f'{tag}: full margin {final["full_margin_bits"]:.9f}; Q1 loss {final["margin_penalty_bits"]:.9g}', flush=True)


if __name__ == '__main__':
    main()

"""Search target-map dense witnesses while retaining verified sparse intervals.

Useful when transporting a larger-state witness yields a weak dense bound.
The complete verifier selects and replays the candidate; a weaker candidate
is archived and the stronger preceding reference is restored.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import complete_bch_reference_batch_v3 as batch
import study_bch_dominance_v1 as study
import transport_bch_full_cover_v1 as transport

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry', type=batch.geometry, action='append', required=True)
    parser.add_argument('--maximum-refinements', type=int, default=300)
    args = parser.parse_args()
    if args.maximum_refinements < 0:
        parser.error('maximum-refinements must be nonnegative')
    requests = []
    for b, t, s, e in dict.fromkeys(args.geometry):
        tag = f'b{b}_t{t}_s{s}_e{e}'
        path = HERE/f'bch_full_reference_{tag}.json'
        reference = json.loads(path.read_text())
        if reference['status'] != 'VERIFIED_BINARY64_FULL_REFERENCE':
            raise ValueError('Expected an existing verified full reference')
        transport.authenticate(reference)
        intervals = transport.select_covers(reference)
        dense = [(p, d) for p, d in intervals if d['status'] == 'BINARY64_COMPLETE_DENSE_INTERVAL']
        if len(dense) != 1:
            raise ValueError('Expected one existing complete dense interval')
        requests.append((b, t, s, e, tag, path, reference, intervals, dense[0][1]['occupation_min']))
    for b, t, s, e, tag, path, reference, intervals, minimum in requests:
        common = ['--block', str(b), '--step', str(t), '--state', str(s), '--exponent', str(e)]

        def run(script, extra, stage):
            print(f'{tag}: {stage}', flush=True)
            with (HERE/f'bch_dense_target_v1_{tag}_{stage}.log').open('w') as output:
                subprocess.run([sys.executable, '-u', str(HERE/script), *common, *extra],
                               cwd=HERE, stdout=output, stderr=subprocess.STDOUT, check=True)

        seed = HERE/f'bch_dense_seed_bulk_{tag}_q{minimum}.json'
        if not seed.exists():
            run('seed_bch_dense_v3.py', ['--minimum', str(minimum), '--nodes', '127' if b == 128 else '63'], 'seed')
        dense_path = HERE/f'bch_dense_v11_{tag}_q{minimum}/cover.json'
        if not dense_path.exists():
            run('close_bch_dense_v11.py', ['--minimum', str(minimum), '--target-bits', '55' if b == 128 else '30',
                                         '--maximum-refinements', str(args.maximum_refinements)], 'dense')
        refined = HERE/f'bch_dominance_v2/{tag}.json'
        if not refined.exists():
            raise ValueError('Existing transported reference needs its refined Q2..4 checkpoint')
        extra = ['--dense-minimum', str(minimum), '--dense-cover', str(dense_path),
                 '--sparse-checkpoint', str(refined)]
        for interval, data in intervals:
            if data['status'] == 'BINARY64_COMPLETE_SPARSE_INTERVAL':
                extra += ['--sparse-cover', str(interval)]
        history = HERE/'bch_full_reference_history'; history.mkdir(exist_ok=True)
        original_bytes = path.read_bytes()
        archive = history/f'{path.stem}_{study.sha(path)}.json'
        if not archive.exists():
            archive.write_bytes(original_bytes)
        run('verify_bch_full_reference_v5.py', extra, 'verify')
        current = json.loads(path.read_text())
        if current['full_margin_bits'] < reference['full_margin_bits']:
            candidate = history/f'{path.stem}_{study.sha(path)}.json'
            if not candidate.exists():
                candidate.write_bytes(path.read_bytes())
            path.write_bytes(original_bytes)
            print(f'{tag}: retained stronger preceding reference', flush=True)
        else:
            print(f'{tag}: full margin {current["full_margin_bits"]:.9f}; '
                  f'Q1 loss {current["margin_penalty_bits"]:.9g} bits', flush=True)


if __name__ == '__main__':
    main()

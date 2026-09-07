"""Produce and replay selected BCH reference geometries, strictly sequentially.

Existing receipts are reused and authenticated by the full verifier. A weak
complete upper bound is retained; reaching the dense target is not required
for coverage. This driver does not turn a binary64 result into a certificate.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def geometry(value):
    try:
        parts = tuple(map(int, value.split(':')))
    except ValueError as error:
        raise argparse.ArgumentTypeError('Expected B:T:S:log2K') from error
    if len(parts) != 4 or parts[0] not in (64, 128) or min(parts[1:]) <= 0:
        raise argparse.ArgumentTypeError('Expected B:T:S:log2K with B=64 or 128')
    return parts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry', type=geometry, action='append', required=True)
    parser.add_argument('--maximum-refinements', type=int, default=300)
    parser.add_argument('--dense-minimum', type=int, default=257)
    args = parser.parse_args()
    if args.maximum_refinements < 0:
        parser.error('maximum-refinements must be nonnegative')
    if args.dense_minimum < 257:
        parser.error('This driver requires a dense cutoff >=257')
    # Validate the whole request before starting any numerical producer.
    for b, t, s, e in args.geometry:
        length = (1 << e)//(b//2)
        if length < 5 or length % t:
            parser.error('This version requires L>=5 and complete epochs')
        if not (HERE/f'bch_dominance_v1/b{b}_t{t}_s{s}_e{e}.json').exists():
            parser.error('The requested geometry needs an original Q2..4 checkpoint first')
    for b, t, s, e in dict.fromkeys(args.geometry):
        length = (1 << e)//(b//2)
        tag = f'b{b}_t{t}_s{s}_e{e}'
        common = ['--block', str(b), '--step', str(t), '--state', str(s), '--exponent', str(e)]

        def run(script, extra, stage, with_geometry=True):
            log = HERE/f'bch_batch_v3_{tag}_{stage}.log'
            print(f'{tag}: {stage}', flush=True)
            with log.open('w') as output:
                subprocess.run([sys.executable, '-u', str(HERE/script), *(common if with_geometry else []), *extra],
                               cwd=HERE, stdout=output, stderr=subprocess.STDOUT, check=True)

        minimum = args.dense_minimum
        seed = HERE/f'bch_dense_seed_bulk_{tag}_q{minimum}.json'
        if length >= minimum and not seed.exists():
            run('seed_bch_dense_v3.py', ['--minimum', str(minimum), '--nodes', '127' if b == 128 else '63'], 'seed')
        dense = HERE/f'bch_dense_v11_{tag}_q{minimum}/cover.json'
        if length >= minimum and not dense.exists():
            run('close_bch_dense_v11.py', ['--minimum', str(minimum), '--target-bits', '55' if b == 128 else '30',
                                         '--maximum-refinements', str(args.maximum_refinements)], 'dense')
        maximum = min(minimum-1, length)
        original = json.loads((HERE/f'bch_dominance_v1/{tag}.json').read_text())
        target = min(55. if b == 128 else 35., original['summary']['q1_margin_bits']+20.)
        sparse = HERE/f'bch_sparse_tail_v5_{tag}_q5_{maximum}/cover.json'
        if not sparse.exists():
            run('close_bch_sparse_tail_v5.py', ['--minimum', '5', '--maximum', str(maximum),
                                              '--target-bits', str(target)], 'sparse')
        refined = HERE/f'bch_dominance_v2/{tag}.json'
        if not refined.exists():
            run('refine_bch_sparse_v2.py', ['--checkpoint', f'{tag}.json'], 'refine', with_geometry=False)
        extra = ['--dense-minimum', str(minimum), '--sparse-cover', str(sparse),
                 '--sparse-checkpoint', str(refined)]
        run('verify_bch_full_reference_v4.py', extra, 'verify')
        result = json.loads((HERE/f'bch_full_reference_{tag}.json').read_text())
        if result['status'] != 'VERIFIED_BINARY64_FULL_REFERENCE':
            raise ValueError('Missing verified full-reference status')
        print(f"{tag}: full margin {result['full_margin_bits']:.9f}; "
              f"Q1 loss {result['margin_penalty_bits']:.9g} bits", flush=True)


if __name__ == '__main__':
    main()

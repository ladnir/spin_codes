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
    args = parser.parse_args()
    if args.maximum_refinements < 0:
        parser.error('maximum-refinements must be nonnegative')
    # Validate the whole request before starting any numerical producer.
    for b, t, s, e in args.geometry:
        length = (1 << e)//(b//2)
        if length < 257 or length % t:
            parser.error('This version requires L>=257 and complete epochs')
    for b, t, s, e in dict.fromkeys(args.geometry):
        tag = f'b{b}_t{t}_s{s}_e{e}'
        common = ['--block', str(b), '--step', str(t), '--state', str(s), '--exponent', str(e)]

        def run(script, extra, stage):
            log = HERE/f'bch_batch_v1_{tag}_{stage}.log'
            print(f'{tag}: {stage}', flush=True)
            with log.open('w') as output:
                subprocess.run([sys.executable, '-u', str(HERE/script), *common, *extra],
                               cwd=HERE, stdout=output, stderr=subprocess.STDOUT, check=True)

        seed = HERE/f'bch_dense_seed_bulk_{tag}_q257.json'
        if not seed.exists():
            run('seed_bch_dense_v3.py', ['--minimum', '257', '--nodes', '127' if b == 128 else '63'], 'seed')
        dense = HERE/f'bch_dense_v11_{tag}_q257/cover.json'
        if not dense.exists():
            run('close_bch_dense_v11.py', ['--minimum', '257', '--target-bits', '55' if b == 128 else '30',
                                         '--maximum-refinements', str(args.maximum_refinements)], 'dense')
        sparse = HERE/f'bch_sparse_tail_v2_{tag}_q5_256/cover.json'
        if not sparse.exists():
            run('close_bch_sparse_tail_v2.py', ['--minimum', '5', '--maximum', '256'], 'sparse')
        run('verify_bch_full_reference_v3.py', [], 'verify')
        result = json.loads((HERE/f'bch_full_reference_{tag}.json').read_text())
        if result['status'] != 'VERIFIED_BINARY64_FULL_REFERENCE':
            raise ValueError('Missing verified full-reference status')
        print(f"{tag}: full margin {result['full_margin_bits']:.9f}; "
              f"Q1 loss {result['margin_penalty_bits']:.9g} bits", flush=True)


if __name__ == '__main__':
    main()

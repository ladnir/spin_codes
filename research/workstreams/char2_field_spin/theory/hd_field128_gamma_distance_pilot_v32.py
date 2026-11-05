"""Remote-only bounded FLOAT distance pilot; never a whole-code certificate.

Project only exact positive witness literals from a pinned source. Rebuild all
coefficients and complete original-weight midpoint-profile sums at the new
distance. Frozen 40% producers and their receipts are not modified or reused.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import resource
import shutil
import signal

ROOT = Path('/tmp/spin-field128-proof-52GBmlZD/bundle/theory')
CONSTRUCTION = 'toeplitz-fullfield-full-rs64-30-gamma16-feedback-alternating-swapped-pairs-k20-rate-half-v1'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--distance', default='9/20')
    parser.add_argument('--mode', choices=('screen', 'search'), default='screen')
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    if Path(__file__).resolve().parent != ROOT or not 2 <= args.r < 8192:
        raise ValueError('Remote research tree and original weight2..8191 required')
    distance = F(args.distance)
    if not F(2, 5) < distance < F(1, 2):
        raise ValueError('Expected a new target strictly between40% and50%')
    if shutil.disk_usage(ROOT).free < 20 << 30:
        raise RuntimeError('Remote free-space reserve')
    resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (160, 165))
    resource.setrlimit(resource.RLIMIT_FSIZE, (1 << 20, 1 << 20))
    signal.alarm(180)
    import hd_field128_gamma_midpoint_dense_neighbors as screen
    import hd_field128_gamma_midpoint_mark_search as lower
    import hd_field128_gamma_midpoint_dense_mark_search as upper
    source = (ROOT/args.source).resolve()
    output = (ROOT/args.output_dir).resolve()
    if not source.is_relative_to(ROOT) or not output.is_relative_to(ROOT) or output == ROOT:
        raise ValueError('Paths outside bounded research tree')
    if source.stat().st_size > 1 << 20:
        raise ValueError('Source larger than1MiB')
    raw = source.read_bytes()
    if sha256(raw).hexdigest() != args.sha256:
        raise ValueError('Source SHA256 mismatch')
    saved = json.loads(raw)
    if saved.get('construction') != CONSTRUCTION or saved.get('state_dimension') != 30:
        raise ValueError('Source is not the same gamma30 construction')
    witness = saved['witness']
    markers, marker_key = lower.enrich.positive_markers(witness)
    z, alpha = F(witness['z']), F(witness['alpha'])
    if not 0 < z < 1 or not 0 < alpha <= 1 or F(float(z)) != z or F(float(alpha)) != alpha:
        raise ValueError('Marks must be positive exact binary64 values')
    grid = witness.get('grid', 128)
    if grid not in (64, 128, 256):
        raise ValueError('Unsupported alias grid')
    output.mkdir(exist_ok=False)
    written = 0

    def write(name, value):
        nonlocal written
        encoded = (json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n').encode()
        written += len(encoded)
        if written > 1 << 20:
            raise RuntimeError('Total JSON output budget1MiB exceeded')
        path = output/name
        with path.open('xb') as stream:
            stream.write(encoded)
        return dict(path=str(path.relative_to(ROOT)), sha256=sha256(encoded).hexdigest(), bytes=len(encoded))

    seed = dict(construction=CONSTRUCTION, state_dimension=30, active_groups=args.r,
        target_distance=str(distance), source_sha256={},
        literal_origin=dict(path=str(source.relative_to(ROOT)), sha256=args.sha256,
            source_marker_key=marker_key, source_active_groups=saved.get('active_groups'),
            source_distance=saved.get('target_distance')),
        witness=dict(z=str(z), alpha=str(alpha), grid=grid,
            markers=[list(map(str, pair)) for pair in markers], state_dimension=30),
        scope='Exact literal proposal only. No inherited numerical bound or certificate.')
    write('seed.json', seed)
    source_seed = output/'seed.json'
    if args.mode == 'screen':
        result = screen.run(source_seed, [args.r], max_seconds=90.)
    else:
        helper = lower if args.r < 4096 else upper
        result = helper.run_search(source_seed, args.r, max_evaluations=24, max_seconds=120.)
    if source.read_bytes() != raw:
        raise RuntimeError('Pinned literal source changed')
    result['pilot_driver_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    result['pilot_limits'] = dict(address_space_bytes=4 << 30, wall_seconds=180,
        cpu_seconds=160, json_bytes=1 << 20, local_generated_data_bytes=0)
    receipt = write('result.json', result)
    if args.mode == 'screen':
        compact = [dict(r=row['active_groups'], log2_upper_screen=row['log2_upper_screen'],
            largest_profile=row['largest_profile']) for row in result['completed_weights']]
    else:
        best = result.get('best')
        compact = None if best is None else dict(r=args.r, log2_upper_screen=best['log2_upper_screen'],
            largest_profile=best['largest_profile'], completed_trials=result['completed_trials'])
    print(json.dumps(dict(receipt=receipt, result=compact, directed=False,
        whole_code_certificate=False), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()

"""Sparse-gap orchestration; frozen arithmetic, explicit original weights."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
import shutil

import hd_field128_gamma_all_pair_certificate as cert
import hd_field128_gamma_midpoint_weight_bands as bands
import hd_field128_gamma_profile_certificate as output

HERE = Path(__file__).resolve()
RESERVE = 750 << 20


def assigned(first, last):
    return (type(first) is int and type(last) is int and first <= last
            and any(lo <= first <= last <= hi
                    for lo, hi in ((65, 100), (140, 240), (273, 383))))


def write_new(path, value):
    path = Path(path).resolve()
    if shutil.disk_usage(path.parent).free < RESERVE:
        raise RuntimeError('Resource stop: fewer than750MiB free; files preserved')
    return output.write_new(path, value)


def positive_proposal(source, expected_hash, weight):
    if not assigned(weight, weight):
        raise ValueError('Weight outside root sparse gaps')
    path = Path(source).resolve()
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != expected_hash:
        raise ValueError('Source byte pin mismatch')
    data = json.loads(raw)
    if (data['active_groups'] != weight or data['state_dimension'] != 30
            or F(data['target_distance']) != F(2, 5)
            or data['construction'] != cert.local.construction(30)):
        raise ValueError('Construction or original weight mismatch')
    w = data['witness']
    z, alpha = F(w['z']), F(w['alpha'])
    if not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError('Invalid z/alpha proposal')
    markers = [[F(2. ** float(F(x))) for x in row]
               for row in w['marker_log_odds_exact']]
    if (not 1 <= len(markers) <= 128 or any(len(row) != 2 or any(
            x <= 0 or not isfinite(float(x)) for x in row) for row in markers)):
        raise ValueError('Invalid positive marker proposal')
    pins = dict(data['source_sha256'])
    pins[HERE.name] = sha256(HERE.read_bytes()).hexdigest()
    for name, digest in pins.items():
        if sha256((HERE.parent / name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Source changed: ' + name)
    if path.read_bytes() != raw:
        raise RuntimeError('Calibration source changed')
    return dict(schema='gamma-v22-root-positive-literal-proposal-v1',
        active_groups=weight, state_dimension=30, target_distance='2/5',
        construction=data['construction'], source_sha256=pins,
        source=dict(path=str(path), sha256=expected_hash),
        witness=dict(z=str(z), alpha=str(alpha), grid=128,
                     markers=[list(map(str, row)) for row in markers], state_dimension=30),
        source_role='Exact positive proposals only; rebuild all numerical caps',
        conversion='Freeze binary64 exponentiation as exact rational literals, not bounds',
        directed=False, whole_code_certificate=False)


def run_band(bank, first, last, folder, witness_source, *, max_nodes=512):
    if type(bank) is not cert.AliasBank:
        raise TypeError('A fresh actual AliasBank is required')
    if (not assigned(first, last) or last - first >= 64
            or type(max_nodes) is not int or max_nodes not in (512, 2048)):
        raise ValueError('Expected an assigned band of at most64weights and bounded nodes')
    driver_bytes = HERE.read_bytes()
    witness_path = Path(witness_source['path']).resolve()
    witness_bytes = witness_path.read_bytes()
    if sha256(witness_bytes).hexdigest() != witness_source['sha256']:
        raise ValueError('Witness byte pin mismatch')
    destination = Path(folder).resolve()
    if shutil.disk_usage(destination.parent).free < RESERVE:
        raise RuntimeError('Resource stop: fewer than750MiB free')

    def check():
        if HERE.read_bytes() != driver_bytes or witness_path.read_bytes() != witness_bytes:
            raise RuntimeError('Driver or witness changed during cover')
        bank.check_pins()

    check()
    destination.mkdir(exist_ok=False)
    try:
        write_new(destination / 'driver.json', dict(
            schema='gamma-v22-root-bounded-band-driver-v1',
            driver=dict(path=str(HERE), sha256=sha256(driver_bytes).hexdigest()),
            witness_source=witness_source, requested_weight_range=[first, last],
            max_nodes=max_nodes, max_prepare_seconds=180, max_cover_seconds=60,
            cell_failure_bits=72, aggregate_failure_bits=54, disk_reserve_bytes=RESERVE,
            numerical_provider='hd_field128_gamma_midpoint_weight_bands.py',
            fresh_bank_rebuilt=False, old_numerical_bounds_reused=False))
        evaluator = bands.Evaluator(bank, first, last, max_seconds=180.,
            witness_source=witness_source,
            progress=lambda row: print(json.dumps(row), flush=True))
        result = evaluator.run_cover(max_nodes=max_nodes, max_seconds=60.,
            cell_failure_bits=72, aggregate_bits=54,
            checkpoint=lambda row: write_new(destination / f"node_{row['sequence']:06d}.json", row),
            context_checkpoint=lambda row: write_new(destination / 'context.json', row))
        check()
        write_new(destination / 'summary.json', result)
        print(json.dumps({key: result[key] for key in ('status', 'requested_weight_range',
            'evaluated_nodes', 'log2_upper_rational', 'stop_reason')}), flush=True)
        return result
    finally:
        check()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--weight', type=int, required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = positive_proposal(args.source, args.sha256, args.weight)
    print(json.dumps(write_new(args.output, result)), flush=True)


if __name__ == '__main__':
    main()

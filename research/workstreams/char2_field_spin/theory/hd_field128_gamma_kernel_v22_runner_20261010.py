"""Guarded orchestration only for the assigned lower gaps; frozen arithmetic."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_all_pair_certificate as cert
import hd_field128_gamma_midpoint_weight_bands as bands
import hd_field128_gamma_profile_certificate as output

HERE = Path(__file__).resolve()
RESERVE = 750 << 20


def assigned(first, last):
    return (type(first) is int and type(last) is int and first <= last
            and (385 <= first <= last <= 510 or 521 <= first <= last <= 757))


def write_new(path, value):
    """Check the reserve immediately before every exclusive artifact write."""
    path = Path(path).resolve()
    if shutil.disk_usage(path.parent).free < RESERVE:
        raise RuntimeError("resource stop: fewer than 750 MiB free; existing files preserved")
    return output.write_new(path, value)


def positive_proposal(source, weight):
    """Freeze finite binary64 exponentiation as literal proposals, not bounds."""
    from math import isfinite
    if not assigned(weight, weight):
        raise ValueError("weight is outside the assigned gaps")
    path = Path(source).resolve()
    raw = path.read_bytes()
    data = json.loads(raw)
    if (data['active_groups'] != weight or data['state_dimension'] != 30
            or F(data['target_distance']) != F(2, 5)
            or data['construction'] != cert.local.construction(30)):
        raise ValueError("construction or original weight mismatch")
    w = data['witness']
    markers = [[F(2. ** float(F(x))) for x in row]
               for row in w['marker_log_odds_exact']]
    if not markers or len(markers) > 128 or any(
            len(row) != 2 or any(x <= 0 or not isfinite(float(x)) for x in row)
            for row in markers):
        raise ValueError("invalid positive marker proposal")
    pins = dict(data['source_sha256'])
    pins[HERE.name] = sha256(HERE.read_bytes()).hexdigest()
    for name, digest in pins.items():
        if sha256((HERE.parent / name).read_bytes()).hexdigest() != digest:
            raise RuntimeError("source changed: " + name)
    if path.read_bytes() != raw:
        raise RuntimeError("calibration source changed")
    return dict(schema='gamma-v22-positive-literal-proposal-v1',
                active_groups=weight, state_dimension=30, target_distance='2/5',
                construction=data['construction'], source_sha256=pins,
                source=dict(path=str(path), sha256=sha256(raw).hexdigest()),
                witness=dict(z=w['z'], alpha=w['alpha'], grid=128,
                             markers=[list(map(str, row)) for row in markers],
                             state_dimension=30),
                source_role='Exact positive proposals only; no numerical bound or coefficient reused',
                conversion='Freeze binary64 2**log as exact positive rational literals; rebuild all caps',
                directed=False, whole_code_certificate=False)


def run_band(bank, first, last, folder, witness_source, *, max_nodes=512):
    if type(bank) is not cert.AliasBank:
        raise TypeError("a fresh actual all-pair AliasBank is required")
    if (not assigned(first, last) or last - first >= 64
            or type(max_nodes) is not int or max_nodes not in (512, 2048)):
        raise ValueError("expected an assigned band of at most 64 weights and bounded nodes")
    driver_bytes = HERE.read_bytes()
    witness_path = Path(witness_source['path']).resolve()
    witness_bytes = witness_path.read_bytes()
    if sha256(witness_bytes).hexdigest() != witness_source['sha256']:
        raise ValueError("witness byte pin mismatch")
    destination = Path(folder).resolve()
    if shutil.disk_usage(destination.parent).free < RESERVE:
        raise RuntimeError("resource stop: fewer than 750 MiB free")

    def check():
        if HERE.read_bytes() != driver_bytes or witness_path.read_bytes() != witness_bytes:
            raise RuntimeError("driver or witness changed during cover")
        bank.check_pins()

    check()
    destination.mkdir(exist_ok=False)
    try:
        write_new(destination / 'driver.json', dict(
            schema='gamma-v22-bounded-band-driver-v1',
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

"""Version 24 portable orchestration with bounded output and frozen arithmetic.

Only the exact seed literals are consumed. Its historical provenance is left
for the local saved-receipt auditor; no numerical bound is imported from it.
All newly emitted provenance links are relative to this theory directory.
This version accepts only 140..383, with at most 64 consecutive weights.
Every artifact write requires 20 GiB free and shares a hard 64 MiB JSON budget.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import re
import shutil
import sys

import flint
import numpy
import scipy

import hd_field128_gamma_all_pair_certificate as cert
import hd_field128_gamma_midpoint_weight_bands as bands
import hd_field128_gamma_profile_certificate as output


HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
RESERVE = 20 << 30
MAX_OUTPUT_BYTES = 64 << 20
FIELD_ORDER = 1 << 128
BUILD_SECONDS, PREPARE_SECONDS, COVER_SECONDS = 600, 180, 60
MAX_TABLE_BYTES = 64 << 20


def relative_name(value, *, suffix=None):
    """One portable basename, not an absolute path or a traversal."""
    if (not isinstance(value, str) or not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', value)
            or value.endswith('.') or (suffix is not None and not value.endswith(suffix))):
        raise ValueError('Expected a portable filename relative to the theory directory')
    return value


def digest_string(value):
    if not isinstance(value, str) or re.fullmatch(r'[0-9a-f]{64}', value) is None:
        raise ValueError('Expected a lowercase SHA-256 hex digest')
    return value


def exact_literal(value):
    """Require canonical rational strings, preserving literal bank metadata."""
    if not isinstance(value, str):
        raise ValueError('Exact witness literals must be rational strings')
    try:
        number = F(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('Invalid exact rational witness literal') from error
    if str(number) != value:
        raise ValueError('Exact rational witness strings must be canonical')
    return number


def validate_seed(seed):
    """Validate parameters without opening or interpreting ancestral links."""
    if not isinstance(seed, dict):
        raise ValueError('The portable seed must be a JSON object')
    if (seed.get('schema') != 'gamma-portable-exact-witness-v1'
            or type(seed.get('field_order')) is not int or seed['field_order'] != FIELD_ORDER
            or type(seed.get('state_dimension')) is not int or seed['state_dimension'] != 30
            or seed.get('target_distance') != '2/5'
            or seed.get('construction') != cert.local.construction(30)
            or type(seed.get('active_groups')) is not int or not 1 <= seed['active_groups'] <= 8192
            or seed.get('directed') is not False or seed.get('whole_code_certificate') is not False):
        raise ValueError('Wrong portable seed construction, parameters, or proposal-only scope')
    witness = seed.get('witness')
    if (not isinstance(witness, dict) or type(witness.get('grid')) is not int
            or witness['grid'] != 128 or type(witness.get('state_dimension')) is not int
            or witness['state_dimension'] != 30):
        raise ValueError('Expected an exact grid-128 state-30 witness')
    z, alpha = exact_literal(witness.get('z')), exact_literal(witness.get('alpha'))
    if not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError('Expected 0 < z < 1 and 0 < alpha <= 1')
    rows = witness.get('markers')
    if not isinstance(rows, list) or not 1 <= len(rows) <= 128:
        raise ValueError('Expected 1..128 exact positive marker pairs')
    markers = []
    for row in rows:
        if not isinstance(row, list) or len(row) != 2:
            raise ValueError('A marker must contain exactly two rational literals')
        pair = tuple(map(exact_literal, row))
        if min(pair) <= 0:
            raise ValueError('Markers must be strictly positive')
        markers.append(pair)
    if len(set(markers)) != len(markers):
        raise ValueError('Duplicate markers would change the emitted literal witness')
    return z, alpha, tuple(markers)


def load_seed(name, expected_sha256):
    name = relative_name(name, suffix='.json')
    digest_string(expected_sha256)
    raw = (HERE/name).read_bytes()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Portable seed byte pin mismatch')
    seed = json.loads(raw)
    parameters = validate_seed(seed)
    return seed, parameters, dict(path=name, sha256=expected_sha256)


def validate_band(first, last, max_nodes):
    if (type(first) is not int or type(last) is not int or not 140 <= first <= last <= 383
            or last-first+1 > 64 or type(max_nodes) is not int or max_nodes not in (512, 2048)):
        raise ValueError('Expected at most 64 assigned weights in 140..383 and 512 or 2048 nodes')


def source_pins():
    pins = bands.source_pins()
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_pins(pins, witness_source):
    for name, expected in pins.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Pinned portable runner source changed: '+name)
    if sha256((HERE/witness_source['path']).read_bytes()).hexdigest() != witness_source['sha256']:
        raise RuntimeError('Pinned portable seed changed')


def reserve_guard(directory, planned_bytes=0):
    if type(planned_bytes) is not int or planned_bytes < 0:
        raise ValueError('Expected a nonnegative integer planned write size')
    if shutil.disk_usage(directory).free-planned_bytes < RESERVE:
        raise RuntimeError('Resource stop: planned output would leave fewer than 20 GiB free; existing artifacts preserved')


class OutputBudget:
    """One new directory's artifact bytes; failed writes reserve their full size.

    The byte count uses exactly the frozen writer's UTF-8 JSON serialization.
    It covers driver, bank, context, node and summary JSON, not external logs.
    Exhaustion raises before opening a file; existing partial evidence remains.
    """

    def __init__(self, directory):
        self.directory = Path(directory).resolve()
        self.used_bytes = 0

    def write_new(self, path, value):
        path = Path(path)
        if path.parent.resolve() != self.directory:
            raise ValueError('An artifact must be directly inside its run directory')
        if path.exists():
            raise FileExistsError(path)
        size = len((json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8'))
        if size > MAX_OUTPUT_BYTES-self.used_bytes:
            raise RuntimeError('Resource stop: 64 MiB total artifact budget; existing partial evidence preserved')
        reserve_guard(path.parent, size)
        # Reserve first: an interrupted/failed write cannot release its charge.
        self.used_bytes += size
        # Ignore the writer's absolute return pin: emit relative links only.
        output.write_new(path, value)


def runtime_metadata():
    """Execution environment observations, not authenticated source pins."""
    return dict(python_version=sys.version, platform=platform.platform(),
        python_flint_version=flint.__version__, flint_version=flint.__FLINT_VERSION__,
        numpy_version=numpy.__version__, scipy_version=scipy.__version__,
        PYTHONPATH=os.environ.get('PYTHONPATH'),
        OPENBLAS_NUM_THREADS=os.environ.get('OPENBLAS_NUM_THREADS'),
        scope='Execution metadata only; not source authentication or a portability proof')


def driver_record(pins, witness_source, first, last, max_nodes, seed_weight):
    return dict(schema='gamma-portable-bounded-band-driver-v24',
        driver=dict(path=NAME, sha256=pins[NAME]), source_sha256=pins,
        witness_source=witness_source, source_calibration_weight=seed_weight,
        requested_weight_range=[first, last], max_nodes=max_nodes,
        max_build_seconds=BUILD_SECONDS, max_prepare_seconds=PREPARE_SECONDS,
        max_cover_seconds=COVER_SECONDS, max_table_bytes=MAX_TABLE_BYTES,
        precision=256, coefficient_bits=80, grid=128,
        field_order=FIELD_ORDER, state_dimension=30, target_distance='2/5',
        construction=cert.local.construction(30), outer_measure='invertible-toeplitz',
        cell_failure_bits=72, aggregate_failure_bits=54, disk_reserve_bytes=RESERVE,
        max_output_bytes=MAX_OUTPUT_BYTES,
        output_budget_scope='Sum of exact UTF-8 exclusive JSON writes in this run; external console logs excluded',
        numerical_provider='hd_field128_gamma_midpoint_weight_bands.py',
        fresh_bank_required=True, fresh_bank_rebuilt=False, old_numerical_bounds_reused=False,
        source_role='Exact seed literals only; fresh original-weight coefficient and band replay',
        directed=False, whole_code_certificate=False)


def run(seed_name, expected_sha256, first, last, output_dir, *, max_nodes=512):
    validate_band(first, last, max_nodes)
    destination = HERE/relative_name(output_dir)
    seed, (z, alpha, markers), witness_source = load_seed(seed_name, expected_sha256)
    pins = source_pins()
    check_pins(pins, witness_source)
    bank = None
    try:
        reserve_guard(HERE, MAX_OUTPUT_BYTES)
        destination.mkdir(exist_ok=False)
        artifacts = OutputBudget(destination)
        artifacts.write_new(destination/'driver.json', driver_record(
            pins, witness_source, first, last, max_nodes, seed['active_groups']))
        progress = lambda row: print(json.dumps(row, allow_nan=False), flush=True)
        bank = cert.AliasBank(z, alpha, markers, state=30, q=FIELD_ORDER,
            distance=F(2, 5), grid=128, precision=256, bits=80,
            outer_measure='invertible-toeplitz', max_seconds=BUILD_SECONDS,
            max_table_bytes=MAX_TABLE_BYTES, progress=progress)
        check_pins(pins, witness_source)
        bank.check_pins()
        artifacts.write_new(destination/'bank.json', dict(bank.metadata(), witness_source=witness_source,
            runner_source_sha256={NAME: pins[NAME]}, fresh_bank_rebuilt=True,
            runtime=runtime_metadata(), directed=False, whole_code_certificate=False))
        evaluator = bands.Evaluator(bank, first, last, max_seconds=PREPARE_SECONDS,
            witness_source=witness_source, progress=progress)
        result = evaluator.run_cover(max_nodes=max_nodes, max_seconds=COVER_SECONDS,
            cell_failure_bits=72, aggregate_bits=54,
            checkpoint=lambda row: artifacts.write_new(destination/f"node_{row['sequence']:06d}.json", row),
            context_checkpoint=lambda row: artifacts.write_new(destination/'context.json', row))
        check_pins(pins, witness_source)
        bank.check_pins()
        artifacts.write_new(destination/'summary.json', result)
        return result
    finally:
        check_pins(pins, witness_source)
        if bank is not None:
            bank.check_pins()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', required=True)
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--last', type=int, required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--max-nodes', type=int, choices=(512, 2048), default=512)
    args = parser.parse_args()
    result = run(args.seed, args.expected_sha256, args.first, args.last,
                 args.output_dir, max_nodes=args.max_nodes)
    print(json.dumps({key: result[key] for key in ('status', 'requested_weight_range',
        'evaluated_nodes', 'log2_upper_rational', 'stop_reason')}, allow_nan=False), flush=True)
    return 0 if result['complete_requested_band'] else 2


if __name__ == '__main__':
    raise SystemExit(main())

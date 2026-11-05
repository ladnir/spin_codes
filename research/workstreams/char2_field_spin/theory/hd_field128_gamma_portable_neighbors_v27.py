"""Bounded portable neighbor FLOAT screens; all numerical helpers stay frozen."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from contextlib import contextmanager
from copy import deepcopy
from hashlib import sha256
import importlib
import json
from math import isfinite
from pathlib import Path
import re
import shutil
import signal
import sys

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
RESERVE, MAX_OUTPUT_BYTES, MEMORY_BYTES = 20 << 30, 1 << 20, 4 << 30
WALL_SECONDS, SCREEN_SECONDS = 120, 90.
HELPERS = (
    ('hd_field128_gamma_portable_band_runner_v24', 'e69c97a06032bbe06aade7d9ed707f9353b633c464e70d740d8647ddea58ed1b'),
    ('hd_field128_gamma30_d40_r913_portable_refit_v24', '73c9dc7dc364b5ab465e799e1a496a09a588df8e94c51521f8007654100c6792'),
    ('hd_field128_gamma_midpoint_neighbors', '108aa7a3e94753d62f2b5a485531673937f2ed4f0d554e4be86dabb42d1f909d'),
)


def validate_weights(weights):
    weights = tuple(weights)
    if (not 1 <= len(weights) <= 8 or len(set(weights)) != len(weights)
            or any(type(r) is not int or not 2 <= r < 4096 for r in weights)):
        raise ValueError('Expected 1..8 distinct original weights in 2..4095')
    return weights


def reserve_guard(directory, planned):
    if type(planned) is not int or planned < 0:
        raise ValueError('Expected nonnegative planned output bytes')
    if shutil.disk_usage(directory).free-planned < RESERVE:
        raise OSError('Resource stop: planned output breaches the 20 GiB reserve')


@contextmanager
def remote_limits():
    """Install process-wide limits before any numerical module imports."""
    if not sys.platform.startswith('linux'):
        raise RuntimeError('Actual screens require Linux; pure tests mock this context')
    import resource
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    ceiling = min([MEMORY_BYTES]+[n for n in (soft, hard) if n != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (ceiling, ceiling))
    previous = signal.getsignal(signal.SIGALRM)
    # Not TimeoutError: frozen helpers may catch cooperative timeouts.
    signal.signal(signal.SIGALRM, lambda unused_signum, unused_frame: os._exit(124))
    signal.setitimer(signal.ITIMER_REAL, WALL_SECONDS)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def helpers():
    for name, digest in HELPERS:
        if sha256((HERE/(name+'.py')).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Frozen helper changed: '+name)
    return tuple(importlib.import_module(name) for name, unused in HELPERS)


def source_pins(loaded):
    portable, normalizer, neighbors = loaded
    pins = {}
    for closure in (portable.source_pins(), normalizer.adaptive.source_pins(),
                    neighbors.enrich.source_pins()):
        for name, digest in closure.items():
            if name in pins and pins[name] != digest:
                raise RuntimeError('Conflicting source pin: '+name)
            pins[name] = digest
    for name, digest in HELPERS:
        pins[name+'.py'] = digest
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_pins(pins, source):
    for name, digest in pins.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Pinned helper changed: '+name)
    if sha256((HERE/source['path']).read_bytes()).hexdigest() != source['sha256']:
        raise RuntimeError('Portable seed changed')


def calibration(seed, source):
    """Project only validated literals, not any ancestral numerical fields."""
    return dict(schema='gamma-portable-neighbor-literal-calibration-v1',
        construction=seed['construction'], state_dimension=30, target_distance='2/5',
        active_groups=seed['active_groups'], source=deepcopy(source), source_sha256={},
        witness=deepcopy({key: seed['witness'][key]
                         for key in ('z', 'alpha', 'markers', 'grid', 'state_dimension')}),
        directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, old_numerical_bounds_reused=False,
        source_role='Exact positive literals only; all requested-weight coefficients and bounds rebuilt')


class Outputs:
    def __init__(self, directory, normalizer):
        self.directory, self.normalize, self.used, self.pins = directory, normalizer, 0, []

    def write(self, name, value):
        if not isinstance(name, str) or re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*\.json', name) is None:
            raise ValueError('Expected one JSON basename')
        path = self.directory/name
        if path.exists():
            raise FileExistsError(path)
        value = self.normalize(value, HERE)
        raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')
        if len(raw) > MAX_OUTPUT_BYTES-self.used:
            raise OSError('Resource stop: 1 MiB cumulative JSON limit')
        reserve_guard(self.directory, len(raw))
        self.used += len(raw)
        with path.open('xb') as stream:
            stream.write(raw)
        pin = dict(path=path.relative_to(HERE).as_posix(), sha256=sha256(raw).hexdigest())
        self.pins.append(pin)
        return pin

    def check(self):
        for pin in self.pins:
            if sha256((HERE/pin['path']).read_bytes()).hexdigest() != pin['sha256']:
                raise RuntimeError('A new immutable output changed')


def validate_result(result, seed, weights):
    """Fail closed on truncated or reordered complete-profile claims."""
    if (result.get('requested_weights') != list(weights)
            or result.get('construction') != seed['construction']
            or result.get('state_dimension') != 30 or result.get('target_distance') != '2/5'
            or result.get('directed') is not False or result.get('whole_code_certificate') is not False
            or result.get('old_numerical_bounds_reused') is not False
            or result.get('no_weight_interpolation') is not True):
        raise ArithmeticError('Wrong FLOAT screen scope')
    if any(result['witness'][key] != seed['witness'][key] for key in ('z', 'alpha', 'markers', 'grid')):
        raise ArithmeticError('Exact parent literals changed')
    rows = result['completed_weights']
    if len(rows) > len(weights):
        raise ArithmeticError('Unexpected completed weights')
    for row, weight in zip(rows, weights):
        if (row.get('active_groups') != weight
                or row.get('profile_pairs_evaluated') != (weight//2+1)**2
                or row.get('complete_selected_weight_screen') is not True
                or row.get('all_profiles_retained') is not True
                or row.get('directed') is not False or row.get('whole_code_certificate') is not False
                or not isfinite(row['log2_upper_screen'])):
            raise ArithmeticError('Incomplete original-weight FLOAT objective')
    pending = list(weights[len(rows):])
    if result.get('pending_weights') != pending:
        raise ArithmeticError('Pending weights differ from the complete prefix')
    if result.get('status') not in ('complete', 'timeout') or (result['status'] == 'complete') != (not pending):
        raise ArithmeticError('Completion status contradicts full-profile coverage')
    return not pending


def run(seed_name, expected_sha256, weights, output_dir):
    weights = validate_weights(weights)
    with remote_limits():
        loaded = helpers()
        portable, normalizer, neighbors = loaded
        seed, unused_parameters, source = portable.load_seed(seed_name, expected_sha256)
        destination = HERE/portable.relative_name(output_dir)
        pins = source_pins(loaded)
        check_pins(pins, source)
        reserve_guard(HERE, MAX_OUTPUT_BYTES)
        destination.mkdir(exist_ok=False)
        outputs = Outputs(destination, normalizer.normalize_new_metadata)
        try:
            outputs.write('driver.json', dict(schema='gamma-portable-neighbor-driver-v27',
                source=source, source_sha256=pins, requested_weights=list(weights),
                screen_seconds=SCREEN_SECONDS, hard_wall_seconds=WALL_SECONDS,
                memory_bytes=MEMORY_BYTES, minimum_free_bytes=RESERVE, maximum_output_bytes=MAX_OUTPUT_BYTES,
                output_budget_scope='All new JSON; caller must separately bound console logs',
                directed=False, whole_code_certificate=False))
            proposal = outputs.write('calibration.json', calibration(seed, source))
            result = neighbors.run(HERE/proposal['path'], weights, max_seconds=SCREEN_SECONDS)
            complete = validate_result(result, seed, weights)
            check_pins(pins, source)
            outputs.check()
            result.update(portable_adapter_schema='gamma-portable-neighbors-FLOAT-v27',
                portable_source=source, adapter_source_sha256=pins, calibration=proposal,
                all_requested_weights_completed=complete,
                hard_wall_seconds=WALL_SECONDS, maximum_output_bytes=MAX_OUTPUT_BYTES)
            pin = outputs.write('neighbors.json', result)
            check_pins(pins, source)
            outputs.check()
            return dict(status=result['status'], result=pin,
                complete_weights=[row['active_groups'] for row in result['completed_weights']],
                pending_weights=result['pending_weights'], whole_code_certificate=False)
        finally:
            check_pins(pins, source)
            outputs.check()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', required=True)
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--weights', nargs='+', type=int, required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    result = run(args.seed, args.expected_sha256, args.weights, args.output_dir)
    print(json.dumps(result, allow_nan=False), flush=True)
    return 0 if result['status'] == 'complete' else 2


if __name__ == '__main__':
    raise SystemExit(main())

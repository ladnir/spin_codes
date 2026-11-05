"""Bounded r913 FLOAT retuning; frozen numerical helpers remain unchanged."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from contextlib import contextmanager
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import importlib
import json
from pathlib import Path
import re
import shutil
import signal
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
PARENT = 'hd_field128_gamma30_d40_r913_portable_refit_v24_20261010/proposal.json'
PARENT_SHA256 = '1f28de050db6812f90fb1c5fef320b01593c4d53c742ef080e7cc12115a2ab79'
BASE = 'hd_field128_gamma30_d40_r913_portable_refit_v24'
BASE_SHA256 = '73c9dc7dc364b5ab465e799e1a496a09a588df8e94c51521f8007654100c6792'
SEARCH = 'hd_field128_gamma_midpoint_mark_search'
SEARCH_SHA256 = 'b753db0272abfcbf8a964cb3298170e75f4e5f8821ca36ceba8e3b6ff5c3f69f'
CONSTRUCTION = 'toeplitz-fullfield-full-rs64-30-gamma16-feedback-alternating-swapped-pairs-k20-rate-half-v1'
RESERVE, MAX_OUTPUT_BYTES, MEMORY_BYTES = 20 << 30, 32 << 20, 4 << 30
WALL_SECONDS = 180
WEIGHT, PROFILE_PAIRS = 913, 457**2
MAX_EVALUATIONS, SEARCH_SECONDS = 32, 90.
MAX_FITS, ADAPTIVE_SECONDS, REPLAY_SECONDS = 8, 45., 20.
TARGET_BITS = 80


def encode(value):
    return (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')


def reserve_guard(directory, planned):
    if type(planned) is not int or planned < 0:
        raise ValueError('Expected nonnegative planned bytes')
    if shutil.disk_usage(directory).free-planned < RESERVE:
        raise OSError('Resource stop: planned output breaches 20 GiB reserve')


@contextmanager
def remote_limits():
    """Linux process limits, installed before importing numerical libraries."""
    if not sys.platform.startswith('linux'):
        raise RuntimeError('Actual run requires Linux resource limits; local tests mock this context')
    import resource
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    limits = [MEMORY_BYTES]+[n for n in (soft, hard) if n != resource.RLIM_INFINITY]
    ceiling = min(limits)
    resource.setrlimit(resource.RLIMIT_AS, (ceiling, ceiling))
    old_handler = signal.getsignal(signal.SIGALRM)
    # A hard deadline must escape the frozen helpers' TimeoutError handlers.
    # Abrupt exit retains immutable trials but does not constitute completion.
    signal.signal(signal.SIGALRM, lambda unused_signum, unused_frame: os._exit(124))
    signal.setitimer(signal.ITIMER_REAL, WALL_SECONDS)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old_handler)


def load_parent():
    raw = (HERE/PARENT).read_bytes()
    if sha256(raw).hexdigest() != PARENT_SHA256:
        raise ValueError('Reviewed v24 exact79 parent changed')
    source = json.loads(raw)
    if (source.get('active_groups') != WEIGHT or source.get('state_dimension') != 30
            or source.get('target_distance') != '2/5' or source.get('construction') != CONSTRUCTION):
        raise ValueError('Wrong original-r913 gamma30/d40 parent')
    witness = source['witness']
    markers = witness['markers']
    if (witness.get('grid') != 128 or witness.get('state_dimension') != 30
            or len(markers) != 79 or any(not isinstance(p, list) or len(p) != 2 for p in markers)):
        raise ValueError('Expected exact79 grid128 marker pairs')
    for literal in (witness['z'], witness['alpha'], *(v for pair in markers for v in pair)):
        if not isinstance(literal, str) or str(F(literal)) != literal or F(literal) <= 0:
            raise ValueError('Expected canonical positive rational literals')
    if not F(witness['z']) < 1 or not F(witness['alpha']) <= 1:
        raise ValueError('Invalid z/alpha domain')
    return source, raw


def helpers():
    for name, digest in ((BASE, BASE_SHA256), (SEARCH, SEARCH_SHA256)):
        if sha256((HERE/(name+'.py')).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Frozen helper changed: '+name)
    return importlib.import_module(BASE), importlib.import_module(SEARCH)


def source_pins(source=None, loaded=None):
    source = load_parent()[0] if source is None else source
    base, search = helpers() if loaded is None else loaded
    pins = base.source_pins(source)
    for name, digest in search.source_pins().items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = digest
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_pins(pins, raw):
    if (HERE/PARENT).read_bytes() != raw:
        raise RuntimeError('Exact parent bytes changed')
    for name, digest in pins.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Helper source changed: '+name)


class Outputs:
    def __init__(self, directory, normalizer):
        self.directory, self.normalize = directory, normalizer
        self.used, self.pins = 0, []

    def write(self, name, value):
        if not isinstance(name, str) or re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*\.json', name) is None:
            raise ValueError('Expected one JSON basename')
        path = self.directory/name
        if path.exists():
            raise FileExistsError(path)
        raw = encode(self.normalize(value, HERE))
        if self.used+len(raw) > MAX_OUTPUT_BYTES:
            raise OSError('Resource stop: 32 MiB cumulative JSON output limit')
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
                raise RuntimeError('New immutable output changed')


def full(row):
    if row is not None and (row.get('profile_pairs_evaluated') != PROFILE_PAIRS
            or row.get('active_groups') != WEIGHT or row.get('complete_selected_weight_screen') is not True
            or row.get('all_profiles_retained') is not True):
        raise ArithmeticError('Incomplete original-r913 objective')
    return row


def same(row, expected):
    return row is not None and all(row[k] == expected[k]
        for k in ('regional_logs_sha256', 'log2_upper_screen'))


def run(output_dir='hd_field128_gamma30_d40_r913_portable_retune_v25_20261010'):
    if (not isinstance(output_dir, str) or re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', output_dir) is None
            or output_dir.endswith('.')):
        raise ValueError('Expected one portable directory basename')
    started = monotonic()
    source, raw = load_parent()
    reserve_guard(HERE, MAX_OUTPUT_BYTES)
    with remote_limits():
        base, search = helpers()
        pins = source_pins(source, (base, search))
        check_pins(pins, raw)
        destination = HERE/output_dir
        destination.mkdir(exist_ok=False)
        outputs = Outputs(destination, base.normalize_new_metadata)
        def save(name, value):
            check_pins(pins, raw)
            return outputs.write(name, value)
        def seconds(cap, reserved):
            available = WALL_SECONDS-5-(monotonic()-started)-reserved
            if available <= 0:
                raise TimeoutError('No stage budget remains')
            return min(cap, available)
        common = dict(construction=CONSTRUCTION, state_dimension=30, target_distance='2/5', active_groups=WEIGHT,
            source_sha256=pins, directed=False, selected_weight_certificate=False, whole_code_certificate=False,
            old_numerical_bounds_reused=False, no_weight_interpolation=True)
        report = dict(common, schema='gamma-r913-portable-retune-links-FLOAT-v25', status='no-complete-replay',
            source=dict(path=PARENT, sha256=PARENT_SHA256), complete_score=None, exact_replay_match=None,
            required_normal_process_exit=True)
        trials = {}
        try:
            save('driver.json', dict(common, schema='gamma-r913-portable-retune-driver-v25',
                parent=report['source'], driver=dict(path=NAME, sha256=pins[NAME]),
                max_search_evaluations=MAX_EVALUATIONS, search_seconds=SEARCH_SECONDS,
                max_adaptive_fits=MAX_FITS, adaptive_seconds=ADAPTIVE_SECONDS, replay_seconds=REPLAY_SECONDS,
                max_calls_per_fit=128, target_bits=TARGET_BITS, whole_job_seconds=WALL_SECONDS,
                memory_address_space_bytes=MEMORY_BYTES, disk_reserve_bytes=RESERVE,
                maximum_output_bytes=MAX_OUTPUT_BYTES, expected_complete_profile_pairs=PROFILE_PAIRS,
                runtime=base.runtime_metadata(), hard_stop_semantics='Exit124 or resource exception is not completion'))
            seed = save('seed.json', dict(common, schema='gamma-original-weight-positive-literal-proposal-v1',
                source=report['source'], witness=deepcopy(source['witness']), marker_count=79,
                source_role='Exact parent literals only; no numerical bounds inherited'))
            def checkpoint(row):
                full(row)
                if row['witness']['markers'] != source['witness']['markers']:
                    raise ArithmeticError('Fixed-marker search changed exact marker literals')
                pin = save(f"evaluation_{row['evaluation_index']:03d}.json", row)
                trials[row['evaluation_index']] = (pin, deepcopy(row))
                return pin
            result = search.run_search(HERE/seed['path'], WEIGHT, max_evaluations=MAX_EVALUATIONS,
                max_seconds=seconds(SEARCH_SECONDS, ADAPTIVE_SECONDS+REPLAY_SECONDS), checkpoint=checkpoint)
            report['search'] = save('search_summary.json', result)
            best = full(result['best'])
            if best is not None:
                best_pin, saved_best = trials[best['evaluation_index']]
                if best != saved_best:
                    raise ArithmeticError('Best search row does not match immutable trial')
                report['best_trial'] = best_pin
                fits = 0 if best['log2_upper_screen'] <= -TARGET_BITS else MAX_FITS
                fitted = base.adaptive.run(HERE/best_pin['path'], WEIGHT,
                    max_seconds=seconds(ADAPTIVE_SECONDS, REPLAY_SECONDS), max_fits=fits, max_calls=128,
                    target_bits=TARGET_BITS)
                fitted_pin = save('adaptive.json', fitted)
                report.update(adaptive=fitted_pin, adaptive_status=fitted['status'])
                baseline = full(fitted['history'][0]) if fitted.get('history') else None
                if not same(baseline, best):
                    raise ArithmeticError('Fresh adaptive baseline differs from best search trial')
                complete = full(fitted['last_complete_step'])
                if complete is None:
                    raise ArithmeticError('Missing completely scored witness')
                indices, witness = fitted['winning_marker_indices_at_last_complete'], fitted['witness']
                if (any(witness[k] != best['witness'][k] for k in ('z', 'alpha', 'grid', 'state_dimension'))
                        or witness['markers'][:len(best['witness']['markers'])] != best['witness']['markers']):
                    raise ArithmeticError('Adaptive stage changed fixed literals or existing marker prefix')
                if (not isinstance(indices, list) or not indices or indices != sorted(set(indices))
                        or any(type(i) is not int or not 0 <= i < len(witness['markers']) for i in indices)):
                    raise ArithmeticError('Invalid completed winner selection')
                proposal = save('proposal.json', dict(common, schema='gamma-selected-positive-literal-proposal-v1',
                    source=fitted_pin, selected_source_marker_indices=indices, source_marker_count=len(witness['markers']),
                    marker_count=len(indices), witness=dict(z=witness['z'], alpha=witness['alpha'],
                    markers=[witness['markers'][i] for i in indices], grid=128, state_dimension=30),
                    source_role='Last completely scored literal subset; fresh replay required'))
                replay = base.adaptive.run(HERE/proposal['path'], WEIGHT,
                    max_seconds=seconds(REPLAY_SECONDS, 0), max_fits=0, max_calls=128, target_bits=TARGET_BITS)
                replay_pin = save('replay.json', replay)
                matched = replay.get('pending_step') is None and same(full(replay['last_complete_step']), complete)
                report.update(proposal=proposal, replay=replay_pin, marker_count=len(indices),
                    exact_replay_match=matched, complete_score=complete['log2_upper_screen'])
                if not matched:
                    raise ArithmeticError('Fresh selected-subset replay mismatch')
                report['status'] = 'FLOAT-target-met' if complete['log2_upper_screen'] <= -TARGET_BITS else 'FLOAT-budget-exhausted'
            check_pins(pins, raw)
            outputs.check()
            report['elapsed_seconds'] = monotonic()-started
            report['bytes_written_before_summary'] = outputs.used
            pin = save('summary.json', report)
            print(json.dumps(dict(summary=pin, status=report['status'], score=report['complete_score'])), flush=True)
            return report
        except (TimeoutError, MemoryError, ArithmeticError) as error:
            # No normal summary or completed result is emitted on a failed replay.
            report.update(status='partial-resource-or-validation-stop', error_type=type(error).__name__,
                error=str(error), complete_score=None, exact_replay_match=False)
            save('partial_summary.json', report)
            return report
        finally:
            check_pins(pins, raw)
            outputs.check()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', default='hd_field128_gamma30_d40_r913_portable_retune_v25_20261010')
    answer = run(parser.parse_args().output_dir)
    if answer['status'] == 'partial-resource-or-validation-stop':
        raise SystemExit(2)

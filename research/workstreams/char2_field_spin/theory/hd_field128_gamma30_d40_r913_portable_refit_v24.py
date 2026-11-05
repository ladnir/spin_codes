"""Portable, bounded r913 FLOAT refit; frozen numerical helpers are unchanged.

Only newly created metadata links inside theory become relative paths. The
byte-pinned parent and its original Windows ancestry are never rewritten.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform
import re
import shutil
import sys

import flint
import numpy
import scipy
import hd_field128_gamma_midpoint_adaptive as adaptive

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
PARENT = 'hd_field128_gamma30_d40_r913_gap_mark_search_20261010/proposal.json'
PARENT_SHA256 = 'fd99500cdaca668b4b59cc17f742f1c567ca1bd4376ad519fb04f9557239aa80'
ADAPTIVE_SHA256 = 'f0cbbe986c67a1de8cc6d49edbe7b60d80136c164cfaf0e9397c90c36d3f44dd'
RESERVE, MAX_OUTPUT_BYTES = 20 << 30, 16 << 20
WEIGHT, PROFILE_PAIRS = 913, 457**2
MAX_FITS, MAX_SECONDS, MAX_CALLS, TARGET_BITS = 15, 90., 128, 4096


def encode(value):
    return (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')


def normalize_new_metadata(value, root=None):
    """Copy a new result; relativize only absolute path+sha links under root."""
    root = HERE.resolve() if root is None else Path(root).resolve()
    if isinstance(value, list):
        return [normalize_new_metadata(child, root) for child in value]
    if not isinstance(value, dict):
        return value
    result = {key: normalize_new_metadata(child, root) for key, child in value.items()}
    if isinstance(value.get('path'), str) and isinstance(value.get('sha256'), str):
        path = Path(value['path'])
        if path.is_absolute():
            try:
                result['path'] = path.resolve().relative_to(root).as_posix()
            except ValueError:
                pass  # External ancestry remains unchanged, never opened here.
    return result


def reserve_guard(directory, planned_bytes):
    if type(planned_bytes) is not int or planned_bytes < 0:
        raise ValueError('Expected a nonnegative planned byte count')
    if shutil.disk_usage(directory).free-planned_bytes < RESERVE:
        raise OSError('Resource stop: planned output would breach the 20 GiB reserve')


class OutputBudget:
    """Exclusive new JSON only; failed writes retain their byte-budget charge."""
    def __init__(self, directory):
        self.directory, self.used_bytes, self.pins = Path(directory).resolve(), 0, []

    def write_new(self, name, value):
        if not isinstance(name, str) or re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*\.json', name) is None:
            raise ValueError('Expected one portable JSON basename')
        path = self.directory/name
        if path.exists():
            raise FileExistsError(path)
        raw = encode(normalize_new_metadata(value))
        if len(raw) > MAX_OUTPUT_BYTES-self.used_bytes:
            raise OSError('Resource stop: 16 MiB cumulative JSON output budget')
        reserve_guard(self.directory, len(raw))
        self.used_bytes += len(raw)
        with path.open('xb') as stream:
            stream.write(raw)
        digest = sha256(raw).hexdigest()
        self.pins.append((path, digest))
        return dict(path=path.relative_to(HERE.resolve()).as_posix(), sha256=digest)

    def check(self):
        for path, digest in self.pins:
            if sha256(path.read_bytes()).hexdigest() != digest:
                raise RuntimeError('A new immutable output changed')


def load_parent():
    raw = (HERE/PARENT).read_bytes()
    if sha256(raw).hexdigest() != PARENT_SHA256:
        raise ValueError('Reviewed exact67 parent changed')
    source = json.loads(raw)
    if (source.get('active_groups') != WEIGHT or source.get('state_dimension') != 30
            or source.get('target_distance') != '2/5'
            or source.get('construction') != adaptive.provider.gamma.construction(30)):
        raise ValueError('Wrong original-weight gamma30/d40 source')
    witness = source['witness']
    if witness.get('grid') != 128 or witness.get('state_dimension') != 30 or len(witness['markers']) != 67:
        raise ValueError('Wrong exact67 grid128 parent witness')
    for value in (witness['z'], witness['alpha'], *(x for pair in witness['markers'] for x in pair)):
        if not isinstance(value, str) or str(F(value)) != value or F(value) <= 0:
            raise ValueError('Noncanonical or nonpositive exact literal')
    if (not F(witness['z']) < 1 or not F(witness['alpha']) <= 1
            or any(not isinstance(pair, list) or len(pair) != 2 for pair in witness['markers'])):
        raise ValueError('Invalid exact witness domain')
    return source, raw


def source_pins(source=None):
    if sha256(Path(adaptive.__file__).read_bytes()).hexdigest() != ADAPTIVE_SHA256:
        raise RuntimeError('Frozen adaptive helper changed')
    source = load_parent()[0] if source is None else source
    pins = adaptive.source_pins()
    for name, digest in source['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = digest
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_pins(pins, raw_parent):
    if (HERE/PARENT).read_bytes() != raw_parent:
        raise RuntimeError('Exact parent bytes changed')
    for name, digest in pins.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Helper source changed: '+name)


def runtime_metadata():
    return dict(python_version=sys.version, platform=platform.platform(),
        python_flint_version=flint.__version__, flint_version=flint.__FLINT_VERSION__,
        numpy_version=numpy.__version__, scipy_version=scipy.__version__,
        PYTHONPATH=os.environ.get('PYTHONPATH'), OPENBLAS_NUM_THREADS=os.environ.get('OPENBLAS_NUM_THREADS'),
        scope='Execution observations only; not authentication or a numerical certificate')


def complete_step(result):
    row = result.get('last_complete_step')
    if row is not None and (row.get('profile_pairs_evaluated') != PROFILE_PAIRS
                            or row.get('complete_selected_weight_screen') is not True):
        raise ArithmeticError('Incomplete original913 profile sum')
    return row


def run(output_dir='hd_field128_gamma30_d40_r913_portable_refit_v24_20261010'):
    if (not isinstance(output_dir, str)
            or re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', output_dir) is None
            or output_dir.endswith('.')):
        raise ValueError('Expected one portable output-directory basename')
    source, raw_parent = load_parent()
    pins = source_pins(source)
    check_pins(pins, raw_parent)
    reserve_guard(HERE, MAX_OUTPUT_BYTES)
    destination = HERE/output_dir
    destination.mkdir(exist_ok=False)
    artifacts = OutputBudget(destination)
    def save(name, value):
        check_pins(pins, raw_parent)
        return artifacts.write_new(name, value)
    try:
        common = dict(construction=source['construction'], state_dimension=30, target_distance='2/5',
            active_groups=WEIGHT, source_sha256=pins, directed=False,
            selected_weight_certificate=False, whole_code_certificate=False, no_weight_interpolation=True)
        save('driver.json', dict(common, schema='gamma-r913-portable-refit-driver-v24',
            driver=dict(path=NAME, sha256=pins[NAME]), parent=dict(path=PARENT, sha256=PARENT_SHA256),
            max_fits=MAX_FITS, max_seconds=MAX_SECONDS, max_calls=MAX_CALLS,
            replay_max_seconds=MAX_SECONDS, target_bits=TARGET_BITS,
            expected_complete_profile_pairs=PROFILE_PAIRS,
            disk_reserve_bytes=RESERVE, maximum_output_bytes=MAX_OUTPUT_BYTES,
            budget_scope='All newly written JSON; external console logs excluded', runtime=runtime_metadata()))
        witness = deepcopy({key: source['witness'][key] for key in ('z', 'alpha', 'markers', 'grid', 'state_dimension')})
        seed_pin = save('seed.json', dict(common, schema='gamma-original-weight-positive-literal-proposal-v1',
            source=dict(path=PARENT, sha256=PARENT_SHA256), witness=witness, marker_count=67,
            old_numerical_bounds_reused=False, source_role='Exact parent literals only; all caps rebuilt'))
        result = adaptive.run(destination/'seed.json', WEIGHT, max_seconds=MAX_SECONDS,
                              max_fits=MAX_FITS, max_calls=MAX_CALLS, target_bits=TARGET_BITS)
        result_pin = save('adaptive.json', result)
        complete = complete_step(result)
        report = dict(common, schema='gamma-r913-portable-refit-links-FLOAT-v24',
            seed=seed_pin, adaptive=result_pin, status=result['status'],
            expected_complete_profile_pairs=PROFILE_PAIRS, exact_replay_match=None,
            complete_score=None, old_numerical_bounds_reused=False)
        if complete is not None:
            indices = result['winning_marker_indices_at_last_complete']
            witness, markers = result['witness'], result['witness']['markers']
            if (not isinstance(indices, list) or not indices or indices != sorted(set(indices))
                    or any(type(i) is not int or not 0 <= i < len(markers) for i in indices)):
                raise ArithmeticError('Invalid completed first-winner selection')
            proposal_pin = save('proposal.json', dict(common,
                schema='gamma-selected-positive-literal-proposal-v1', source=result_pin,
                selected_source_marker_indices=indices, source_marker_count=len(markers), marker_count=len(indices),
                witness=dict(z=witness['z'], alpha=witness['alpha'], markers=[markers[i] for i in indices],
                             grid=witness['grid'], state_dimension=30),
                source_role='Exact literal selection from last complete step; fresh replay required',
                old_numerical_bounds_reused=False))
            replay = adaptive.run(destination/'proposal.json', WEIGHT, max_seconds=MAX_SECONDS,
                                  max_fits=0, max_calls=MAX_CALLS, target_bits=TARGET_BITS)
            fresh = complete_step(replay)
            matched = (replay.get('pending_step') is None and fresh is not None
                and fresh['regional_logs_sha256'] == complete['regional_logs_sha256']
                and fresh['log2_upper_screen'] == complete['log2_upper_screen'])
            replay['independent_replay_comparison'] = dict(source=result_pin,
                complete_score_and_regional_hash_match=matched, comparison_is_not_a_directed_certificate=True)
            replay_pin = save('replay.json', replay)
            report.update(proposal=proposal_pin, replay=replay_pin, marker_count=len(indices),
                          complete_score=complete['log2_upper_screen'], exact_replay_match=matched)
            if not matched:
                save('failed_replay_summary.json', report)
                raise ArithmeticError('Fresh subset replay differs from last complete step')
        pin = save('summary.json', report)
        print(json.dumps(dict(summary=pin, status=report['status'], score=report['complete_score'],
                              exact_replay_match=report['exact_replay_match'])), flush=True)
        return report
    finally:
        check_pins(pins, raw_parent)
        artifacts.check()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', default='hd_field128_gamma30_d40_r913_portable_refit_v24_20261010')
    run(parser.parse_args().output_dir)

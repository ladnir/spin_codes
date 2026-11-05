"""Fresh913 endpoint calibration; signed geometric marks are proposals only."""

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from shutil import disk_usage

import hd_field128_gamma_midpoint_adaptive as adaptive
import hd_field128_gamma_midpoint_calibrate as calibrate


PARENTS = {
    1000: ('hd_field128_gamma30_d40_r1000_recalibration_cycles_20261010/local_mark_search/evaluation_038.json',
           '51951e40daaaa95df10c7bf3df5afbab3b37e77793a3588e9bc818b0fda1d97d'),
    1200: ('hd_field128_gamma30_d40_r1200_retained_bank_20261010/proposal.json',
           'e6a55e8891b2032298c5dfd3436e386867b3b17c4d6ffbc20ec3300caa8903a8')}


def write_new(path, value):
    if disk_usage(Path(__file__).resolve().parent).free < 750 * 1024**2:
        raise OSError('Artifact guard: less than750MiB disk space remains')
    return adaptive.provider.search.write_new(path, value)


def run(weight):
    if weight not in (913,):
        raise ValueError('Only the assigned913 endpoint is allowed')
    low, high = 1000, 1200
    fraction = F(weight-low, high-low)
    here, runner = Path(__file__).resolve().parent, Path(__file__).resolve()
    pins, parents, rows = adaptive.source_pins(), {}, []
    def merge(added):
        for name, digest in added.items():
            if name in pins and pins[name] != digest:
                raise RuntimeError('Conflicting source pin: '+name)
            pins[name] = digest
    merge(calibrate.source_pins())
    sources = []
    for calibration_weight in (low, high):
        name, digest = PARENTS[calibration_weight]
        path = here/name
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise RuntimeError('Pinned neighboring positive proposal changed')
        data = json.loads(raw)
        witness = data['witness']
        z, alpha = F(witness['z']), F(witness['alpha'])
        if (not 0 < z < 1 or not 0 < alpha <= 1 or F(float(z)) != z or F(float(alpha)) != alpha
                or not witness['markers']
                or any(len(pair) != 2 or any(F(x) <= 0 for x in pair) for pair in witness['markers'])):
            raise ValueError('Invalid exact positive parent literals')
        parents[path] = raw
        sources.append(dict(path=str(path), sha256=digest, source_calibration_weight=calibration_weight))
        rows.append(witness)
        merge(data.get('source_sha256', {}))
    pins[runner.name] = sha256(runner.read_bytes()).hexdigest()
    def check():
        if any(path.read_bytes() != raw for path, raw in parents.items()):
            raise RuntimeError('Neighboring literal proposal changed')
        for name, digest in pins.items():
            if sha256((here/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Source closure changed: '+name)
    check()
    marks = [2.**((1-float(fraction))*log2(float(F(rows[0][key])))
                  +float(fraction)*log2(float(F(rows[1][key])))) for key in ('z', 'alpha')]
    if not 0 < marks[0] < 1 or not 0 < marks[1] <= 1:
        raise ArithmeticError('Invalid finite geometric mark proposal')
    output = here/f'hd_field128_gamma30_d40_r{weight}_gap_calibration_20261010'
    output.mkdir(exist_ok=False)
    provenance = dict(sources=sources, high_fraction=str(fraction),
        rule='Signed geometric z/alpha proposal only; negative fraction is extrapolation of literals, never bounds',
        actual_binary64_z=str(F(marks[0])), actual_binary64_alpha=str(F(marks[1])))
    try:
        fresh = calibrate.run(weight, *marks, max_seconds=90., source=provenance)
    except (TimeoutError, MemoryError) as error:
        check()
        write_new(output/'calibration_failure.json', dict(
            active_groups=weight, source=provenance, source_sha256=dict(sorted(pins.items())),
            error_type=type(error).__name__, message=str(error), complete_selected_weight_screen=False,
            directed=False, whole_code_certificate=False))
        raise
    fresh['source_sha256'] = dict(sorted(pins.items()))
    check()
    fresh_pin = write_new(output/'fresh_calibration.json', fresh)
    logs = fresh['witness']['marker_log_odds_exact']
    markers = [[str(F(2.**float(F(value)))) for value in pair] for pair in logs]
    if (not 1 <= len(markers) <= 128
            or any(len(pair) != 2 or any(F(x) <= 0 or not isfinite(float(F(x))) for x in pair) for pair in markers)):
        raise ArithmeticError('Invalid exact positive marker conversion')
    seed = dict(schema='gamma-original-weight-positive-literal-proposal-v1',
        construction=adaptive.provider.gamma.construction(30), state_dimension=30,
        target_distance='2/5', active_groups=weight, source=fresh_pin,
        source_sha256=dict(sorted(pins.items())), source_marker_logs_exact=logs,
        conversion='F(2.**float(F(log))) for each saved log literal; all resulting caps rebuilt',
        witness=dict(z=str(F(marks[0])), alpha=str(F(marks[1])), markers=markers, grid=128, state_dimension=30),
        marker_count=len(markers), old_numerical_bounds_reused=False, directed=False,
        whole_code_certificate=False, no_weight_interpolation=True)
    seed_pin = write_new(output/'seed.json', seed)
    result = adaptive.run(output/'seed.json', weight, max_seconds=90., max_fits=17,
                          max_calls=128, target_bits=4096)
    check()
    result_pin = write_new(output/'adaptive.json', result)
    complete = result['last_complete_step']
    expected_profiles = (weight//2+1)**2
    report = dict(schema='gamma-assigned-gap-calibration-links-FLOAT-v1', active_groups=weight,
        construction=seed['construction'], state_dimension=30, target_distance='2/5',
        fresh_calibration=fresh_pin, seed=seed_pin, adaptive=result_pin,
        source_sha256=dict(sorted(pins.items())), status=result['status'],
        directed=False, selected_weight_certificate=False, whole_code_certificate=False,
        no_weight_interpolation=True, expected_complete_profile_pairs=expected_profiles)
    if complete is not None:
        if complete['profile_pairs_evaluated'] != expected_profiles:
            raise ArithmeticError('Missing profiles in completed original-weight sum')
        report['complete_score'] = complete['log2_upper_screen']
        indices = result['winning_marker_indices_at_last_complete']
        witness, positives = result['witness'], result['witness']['markers']
        if (not indices or indices != sorted(set(indices))
                or any(type(i) is not int or not 0 <= i < len(positives) for i in indices)):
            raise ArithmeticError('Invalid complete first-winner selection')
        proposal = dict(schema='gamma-selected-positive-literal-proposal-v1',
            construction=seed['construction'], state_dimension=30, target_distance='2/5', active_groups=weight,
            source=result_pin, source_sha256=dict(sorted(pins.items())),
            selected_source_marker_indices=indices, source_marker_count=len(positives), marker_count=len(indices),
            witness=dict(z=witness['z'], alpha=witness['alpha'], markers=[positives[i] for i in indices],
                         grid=witness['grid'], state_dimension=30),
            source_role='Exact positive literal selection only; no numerical bound imported',
            old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
            whole_code_certificate=False, no_weight_interpolation=True)
        proposal_pin = write_new(output/'proposal.json', proposal)
        replay = adaptive.run(output/'proposal.json', weight, max_seconds=90., max_fits=0,
                              max_calls=128, target_bits=4096)
        row = replay['last_complete_step']
        matched = (replay['pending_step'] is None and row is not None
                   and row['profile_pairs_evaluated'] == expected_profiles
                   and row['regional_logs_sha256'] == complete['regional_logs_sha256']
                   and row['log2_upper_screen'] == complete['log2_upper_screen'])
        replay['independent_replay_comparison'] = dict(source=result_pin,
            expected_regional_logs_sha256=complete['regional_logs_sha256'],
            expected_log2_upper_screen=complete['log2_upper_screen'],
            complete_score_and_regional_hash_match=matched, comparison_is_not_a_directed_certificate=True)
        check()
        replay_pin = write_new(output/'replay.json', replay)
        report.update(proposal=proposal_pin, replay=replay_pin, marker_count=len(indices), exact_replay_match=matched)
        if not matched:
            write_new(output/'failed_replay_summary.json', report)
            raise ArithmeticError('Fresh subset replay differs from completed source')
    check()
    pin = write_new(output/'summary.json', report)
    print(json.dumps(dict(summary=pin, status=report['status'], score=report.get('complete_score'),
                          markers=report.get('marker_count'), exact_replay_match=report.get('exact_replay_match'))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weight', type=int, choices=(913,), required=True)
    run(parser.parse_args().weight)

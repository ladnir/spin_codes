"""Bounded fresh outer-edge calibration from exact r6400 marker literals only."""

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_dense_adaptive as adaptive


def run(weight):
    if weight not in (5984, 6677):
        raise ValueError('Only the two explicitly requested endpoint calibrations are allowed')
    here, runner = Path(__file__).resolve().parent, Path(__file__).resolve()
    source_path = here/'hd_field128_gamma30_d40_r6400_winner_replay_20261010/proposal.json'
    raw = source_path.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != 'd7719fe4a2fafbc3a926bbb7189524094cb8fc45821e308cb2d676d941d7826f':
        raise RuntimeError('Reviewed exact12 parent changed')
    source = json.loads(raw)
    if (source['active_groups'] != 6400 or source['state_dimension'] != 30
            or source['target_distance'] != '2/5'
            or source['construction'] != adaptive.provider.gamma.construction(30)
            or len(source['witness']['markers']) != 12):
        raise ValueError('Unexpected literal parent scope')
    pins = adaptive.source_pins()
    for name, value in source['source_sha256'].items():
        if name in pins and pins[name] != value:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = value
    pins[runner.name] = sha256(runner.read_bytes()).hexdigest()
    def check():
        if source_path.read_bytes() != raw:
            raise RuntimeError('Exact literal parent changed')
        for name, value in pins.items():
            if sha256((here/name).read_bytes()).hexdigest() != value:
                raise RuntimeError('Source closure changed: '+name)
    check()
    output = here/f'hd_field128_gamma30_d40_r{weight}_endpoint_calibration_20261010'
    output.mkdir(exist_ok=False)
    seed = dict(schema='gamma-original-weight-positive-literal-proposal-v1',
        construction=source['construction'], state_dimension=30, target_distance='2/5', active_groups=weight,
        source=dict(path=str(source_path), sha256=digest), source_sha256=dict(sorted(pins.items())),
        source_calibration_active_groups=6400, new_requested_active_groups=weight,
        witness=dict(source['witness']), marker_count=12,
        source_role='Exact positive literals only; no numerical bound or calibration guarantee transferred',
        old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, no_weight_interpolation=True)
    seed_pin = adaptive.provider.search.write_new(output/'seed.json', seed)
    result = adaptive.run(output/'seed.json', weight, max_seconds=90., max_fits=17,
                          max_calls=128, target_bits=4096)
    check()
    result_pin = adaptive.provider.search.write_new(output/'adaptive.json', result)
    complete = result['last_complete_step']
    expected_profiles = (weight//2-max(0, weight-4096)+1)**2
    report = dict(schema='gamma-original-weight-endpoint-calibration-links-FLOAT-v1', active_groups=weight,
        construction=source['construction'], state_dimension=30, target_distance='2/5', seed=seed_pin,
        adaptive=result_pin, source_sha256=dict(sorted(pins.items())), status=result['status'],
        directed=False, selected_weight_certificate=False, whole_code_certificate=False,
        no_weight_interpolation=True, expected_complete_profile_pairs=expected_profiles)
    if complete is not None:
        if complete['profile_pairs_evaluated'] != expected_profiles:
            raise ArithmeticError('Original-r domain was not completely summed')
        report['complete_score'] = complete['log2_upper_screen']
        if complete['log2_upper_screen'] <= -80:
            indices = result['winning_marker_indices_at_last_complete']
            witness, markers = result['witness'], result['witness']['markers']
            if (not indices or indices != sorted(set(indices))
                    or any(type(i) is not int or not 0 <= i < len(markers) for i in indices)):
                raise ArithmeticError('Invalid complete first-winner selection')
            proposal = dict(schema='gamma-dense-selected-positive-literal-proposal-v1',
                construction=source['construction'], state_dimension=30, target_distance='2/5', active_groups=weight,
                source=result_pin, source_sha256=dict(sorted(pins.items())),
                selected_source_marker_indices=indices, source_marker_count=len(markers), marker_count=len(indices),
                witness=dict(z=witness['z'], alpha=witness['alpha'], markers=[markers[i] for i in indices],
                             grid=witness['grid'], state_dimension=30),
                source_role='Exact positive literal selection only; no numerical cap or score inherited',
                selection_role='Fresh complete first-winner selection, pending independent subset replay',
                old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
                whole_code_certificate=False, no_weight_interpolation=True)
            proposal_pin = adaptive.provider.search.write_new(output/'proposal.json', proposal)
            replay = adaptive.run(output/'proposal.json', weight, max_seconds=90., max_fits=0,
                                  max_calls=128, target_bits=4096)
            fresh = replay['last_complete_step']
            matched = (replay['pending_step'] is None and fresh is not None
                       and fresh['profile_pairs_evaluated'] == expected_profiles
                       and fresh['regional_logs_sha256'] == complete['regional_logs_sha256']
                       and fresh['log2_upper_screen'] == complete['log2_upper_screen'])
            replay['independent_replay_comparison'] = dict(source=result_pin,
                expected_regional_logs_sha256=complete['regional_logs_sha256'],
                expected_log2_upper_screen=complete['log2_upper_screen'],
                complete_score_and_regional_hash_match=matched, comparison_is_not_a_directed_certificate=True)
            check()
            replay_pin = adaptive.provider.search.write_new(output/'replay.json', replay)
            report.update(proposal=proposal_pin, replay=replay_pin, marker_count=len(indices), exact_replay_match=matched)
            if not matched:
                adaptive.provider.search.write_new(output/'failed_replay_summary.json', report)
                raise ArithmeticError('Fresh subset replay differs from completed source')
    check()
    pin = adaptive.provider.search.write_new(output/'summary.json', report)
    print(json.dumps(dict(summary=pin, status=report['status'], score=report.get('complete_score'),
                          markers=report.get('marker_count'), exact_replay_match=report.get('exact_replay_match'))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weight', type=int, choices=(5984, 6677), required=True)
    run(parser.parse_args().weight)

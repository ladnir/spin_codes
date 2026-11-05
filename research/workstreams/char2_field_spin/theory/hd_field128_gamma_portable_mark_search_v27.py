"""Portable fixed-marker FLOAT search and fresh one-weight replay, never a certificate.

Only exact seed literals are retargeted. Frozen search optimizes the complete
original-weight midpoint objective; the final replay rebuilds its coefficients.
No trial files: the bounded search summary retains all completed-trial metadata.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import importlib
import json
from math import isfinite
from pathlib import Path

import hd_field128_gamma_portable_neighbors_v27 as portable

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
PORTABLE_SHA256 = '169ee5c7719873ce4528ec9e9a5eea1e87e8c866ae32e713abde0d3678339a94'
SEARCH_NAME = 'hd_field128_gamma_midpoint_mark_search'
SEARCH_SHA256 = 'b753db0272abfcbf8a964cb3298170e75f4e5f8821ca36ceba8e3b6ff5c3f69f'
MAX_BANKS, SEARCH_SECONDS, REPLAY_SECONDS = 32, 75., 15.


def validate_target(weight):
    if type(weight) is not int or not 2 <= weight <= 4095:
        raise ValueError('Expected one original target weight in 2..4095')
    return weight


def helpers():
    for name, digest in ((portable.NAME, PORTABLE_SHA256), (SEARCH_NAME+'.py', SEARCH_SHA256)):
        if sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Frozen helper changed: '+name)
    loaded = portable.helpers()
    return loaded, importlib.import_module(SEARCH_NAME)


def source_pins(loaded, search):
    pins = portable.source_pins(loaded)
    for name, digest in search.source_pins().items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Conflicting source pin: '+name)
        pins[name] = digest
    pins[SEARCH_NAME+'.py'] = SEARCH_SHA256
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def calibration(seed, source, weight):
    validate_target(weight)
    result = portable.calibration(seed, source)
    result.update(schema='gamma-portable-retargeted-mark-calibration-v27',
        active_groups=weight, source_calibration_active_groups=seed['active_groups'],
        retargeted_literal_proposal=True, no_weight_interpolation=True)
    # The frozen search uses a bounded binary64 chart; do not round its initial
    # marks silently. Exact positive marker literals themselves stay untouched.
    z, alpha = (F(result['witness'][key]) for key in ('z', 'alpha'))
    if (not isfinite(float(z)) or not isfinite(float(alpha))
            or F(float(z)) != z or F(float(alpha)) != alpha
            or not 2.**-1022 <= float(z) <= 2.**-1e-8
            or not 2.**-32 <= float(alpha) <= 1):
        raise ValueError('Initial z/alpha must be exact binary64 values in the frozen search chart')
    return result


def validate_search(result, proposal, weight):
    """Check the completed best trial and its fixed-marker search lineage."""
    expected = dict(schema='gamma-midpoint-fixed-marker-mark-search-FLOAT-v1',
        active_groups=weight, construction=proposal['construction'], state_dimension=30,
        target_distance='2/5', directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, old_numerical_bounds_reused=False,
        common_potential_used=False, max_evaluations=MAX_BANKS, max_seconds=SEARCH_SECONDS)
    if any(result.get(k) != v for k, v in expected.items()):
        raise ArithmeticError('Wrong fixed-marker search scope')
    witness = proposal['witness']
    if (result.get('initial_witness') != witness
            or result.get('fixed_positive_marker_literals') != witness['markers']
            or result.get('initial_z_preserved_exactly') != witness['z']
            or result.get('initial_alpha_preserved_exactly') != witness['alpha']):
        raise ArithmeticError('Search changed initial literals or fixed marker bank')
    history, best, attempted = result['history'], result['best'], result['attempted_banks']
    if (type(attempted) is not int or not 0 <= attempted <= MAX_BANKS
            or len(history) != result['completed_trials'] or len(history) > attempted):
        raise ArithmeticError('Invalid search bank accounting')
    indices = [row['evaluation_index'] for row in history]
    if (len(set(indices)) != len(indices)
            or any(type(i) is not int or not 0 <= i < attempted for i in indices)
            or any(not isfinite(row['log2_upper_screen']) for row in history)):
        raise ArithmeticError('Invalid completed trial history')
    if not history:
        if best is not None:
            raise ArithmeticError('Best trial without complete history')
        return None
    checks = dict(active_groups=weight, construction=proposal['construction'], state_dimension=30,
        target_distance='2/5', complete_selected_weight_screen=True, all_profiles_retained=True,
        profile_pairs_evaluated=(weight//2+1)**2, directed=False,
        selected_weight_certificate=False, whole_code_certificate=False,
        old_numerical_bounds_reused=False, common_potential_used=False)
    if not isinstance(best, dict) or any(best.get(k) != v for k, v in checks.items()):
        raise ArithmeticError('Best trial is not a complete original-weight objective')
    if (best['log2_upper_screen'] != min(row['log2_upper_screen'] for row in history)
            or best['witness']['markers'] != witness['markers']
            or best['witness']['grid'] != 128 or best['witness']['state_dimension'] != 30):
        raise ArithmeticError('Best selection or fixed witness differs')
    matching = [row for row in history if row['evaluation_index'] == best['evaluation_index']]
    if (len(matching) != 1 or any(matching[0][key] != best['witness'][key] for key in ('z', 'alpha'))
            or matching[0]['log2_upper_screen'] != best['log2_upper_screen']):
        raise ArithmeticError('Best trial does not match history')
    for key in ('z', 'alpha'):
        value = F(best['witness'][key])
        if str(value) != best['witness'][key] or not 0 < value <= 1 or F(float(value)) != value:
            raise ArithmeticError('Invalid exact best mark')
    if F(best['witness']['z']) >= 1:
        raise ArithmeticError('Best z must be less than one')
    return best


def best_proposal(proposal, best, search_pin):
    result = deepcopy(proposal)
    result.update(schema='gamma-portable-best-fixed-marker-proposal-v27', source=deepcopy(search_pin),
        source_sha256={}, witness=deepcopy({key: best['witness'][key]
            for key in ('z', 'alpha', 'markers', 'grid', 'state_dimension')}),
        selected_search_evaluation=best['evaluation_index'],
        source_role='Selected exact literals only; fresh no-fit replay required; no bound inherited')
    return result


def run(seed_name, expected_sha256, weight, output_dir):
    validate_target(weight)
    # Reuse the reviewed hard alarm/address-space guard before heavy imports.
    with portable.remote_limits():
        loaded, search = helpers()
        validator, normalizer, neighbors = loaded
        seed, unused_parameters, source = validator.load_seed(seed_name, expected_sha256)
        proposal = calibration(seed, source, weight)
        destination = HERE/validator.relative_name(output_dir)
        pins = source_pins(loaded, search)
        portable.check_pins(pins, source)
        portable.reserve_guard(HERE, portable.MAX_OUTPUT_BYTES)
        destination.mkdir(exist_ok=False)
        outputs = portable.Outputs(destination, normalizer.normalize_new_metadata)
        try:
            outputs.write('driver.json', dict(schema='gamma-portable-mark-search-driver-v27',
                source=source, source_sha256=pins, target_weight=weight,
                source_calibration_active_groups=seed['active_groups'], max_search_banks=MAX_BANKS,
                search_seconds=SEARCH_SECONDS, fresh_replay_seconds=REPLAY_SECONDS,
                maximum_fresh_replay_banks=1, hard_wall_seconds=portable.WALL_SECONDS,
                memory_bytes=portable.MEMORY_BYTES, minimum_free_bytes=portable.RESERVE,
                maximum_output_bytes=portable.MAX_OUTPUT_BYTES,
                output_budget_scope='All new JSON; caller must separately bound console logs',
                trial_files_written=False, directed=False, whole_code_certificate=False))
            calibration_pin = outputs.write('calibration.json', proposal)
            result = search.run_search(HERE/calibration_pin['path'], weight,
                max_evaluations=MAX_BANKS, max_seconds=SEARCH_SECONDS)
            best = validate_search(result, proposal, weight)
            portable.check_pins(pins, source)
            outputs.check()
            result.update(portable_adapter_schema='gamma-portable-mark-search-FLOAT-v27',
                          adapter_source_sha256=pins, portable_source=source)
            search_pin = outputs.write('search.json', result)
            selected_pin = replay_pin = None
            complete, same_replay = False, False
            status = 'search-incomplete'
            if best is not None:
                selected = best_proposal(proposal, best, search_pin)
                selected_pin = outputs.write('best_proposal.json', selected)
                replay = neighbors.run(HERE/selected_pin['path'], (weight,), max_seconds=REPLAY_SECONDS)
                complete = portable.validate_result(replay, selected, (weight,))
                if complete:
                    row = replay['completed_weights'][0]
                    same_replay = (row['log2_upper_screen'] == best['log2_upper_screen']
                        and row['regional_logs_sha256'] == best['regional_logs_sha256'])
                status = ('complete' if same_replay else 'replay-mismatch') if complete else 'replay-incomplete'
                portable.check_pins(pins, source)
                outputs.check()
                replay.update(portable_adapter_schema='gamma-portable-mark-search-replay-FLOAT-v27',
                    adapter_source_sha256=pins, selected_proposal=selected_pin,
                    exact_score_and_regional_hash_match_search=same_replay)
                replay_pin = outputs.write('replay.json', replay)
            summary = dict(schema='gamma-portable-mark-search-summary-v27', status=status,
                construction=seed['construction'], state_dimension=30, target_distance='2/5',
                active_groups=weight, source=source, source_sha256=pins,
                calibration=calibration_pin, search=search_pin, best_proposal=selected_pin,
                replay=replay_pin, complete_original_weight_replay=complete,
                exact_score_and_regional_hash_match_search=same_replay,
                log2_upper_screen=best['log2_upper_screen'] if same_replay else None,
                profile_pairs_evaluated=(weight//2+1)**2 if complete else None,
                selected_weight_certificate=False, whole_code_certificate=False, directed=False,
                scope='Selected original-weight FLOAT diagnostic only; no directed or neighboring-weight claim')
            portable.check_pins(pins, source)
            outputs.check()
            summary_pin = outputs.write('summary.json', summary)
            return dict(status=status, summary=summary_pin, best_proposal=selected_pin,
                        log2_upper_screen=summary['log2_upper_screen'], whole_code_certificate=False)
        finally:
            portable.check_pins(pins, source)
            outputs.check()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', required=True)
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--target-weight', type=int, required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    result = run(args.seed, args.expected_sha256, args.target_weight, args.output_dir)
    print(json.dumps(result, allow_nan=False), flush=True)
    return 0 if result['status'] == 'complete' else 2


if __name__ == '__main__':
    raise SystemExit(main())

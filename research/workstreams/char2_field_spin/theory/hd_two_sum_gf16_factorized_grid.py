"""Nine fixed z/p controls of the fullr18..32 factorized G16 bound.

No mark/potential fitting. Three regional banks serve nine exact witnesses.
The first failed inherited-witness artifact and helper remain pinned.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize_scalar

from hd_two_sum_gf16_factorized import (
    forced_nonzero_regions, contraction_envelope, dense_arity_logs,
    group_bounds, all_cohort_band, logsum)


DIRECTORY = Path(__file__).resolve().parent
FIRST_NAME = 'hd_two_sum_gf16_factorized_r18_32_frozen_float.json'
FIRST_SHA256 = 'c898f203b9656fc86cb41dc7acc05fa1d5a5b4118a989b15f9933209bfa93175'
HELPER_NAME = 'hd_two_sum_gf16_factorized.py'
HELPER_SHA256 = 'd34f054a6a1a07cec5ae97c75acdf475763a5c5e5ccab85d8be05cdd281e5da2'
MARKS = (F(99, 100), F(49, 50), F(24, 25))
LOG_POTENTIALS = (F(-64), F(-96), F(-128))


def source_scope():
    for name, digest in ((FIRST_NAME, FIRST_SHA256), (HELPER_NAME, HELPER_SHA256)):
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Pinned first control changed: {name}')
    result = json.loads((DIRECTORY/FIRST_NAME).read_bytes())
    if (result['active_group_band'] != [18, 32] or result['directed'] is not False
            or result['all_constituent_row_sets'] is not True
            or result['dense_arities_summed_after_group_word_power'] is not True):
        raise ValueError('Wrong inherited all-cohort scope')
    for key in ('entries', 'affine_trials', 'exact_dense_moment_sha256', 'group_log_bounds',
                'regional_log_contractions', 'failure_log2_screen', 'envelope_minimum_slack',
                'log2_a_exact_float', 'log2_b_exact_float', 'affine_slope_bounds',
                'affine_search_success', 'affine_evaluations', 'winning_slope_hits_boundary',
                'raw_regional_banks', 'bank_seconds', 'elapsed_seconds', 'source_sha256', 'status',
                'mark', 'log2_potential'):
        result.pop(key, None)
    result.update(complete_requested_band=False, first_control_filename=FIRST_NAME,
        first_control_sha256=FIRST_SHA256, factorized_helper_sha256=HELPER_SHA256,
        fixed_marks=[str(z) for z in MARKS], fixed_log2_potentials=[str(p) for p in LOG_POTENTIALS],
        expected_witnesses=9, z_p_fits=0, reused_regional_banks=True)
    return result


def evaluate(regions, weights, mark, log_potential, check_time):
    caps = contraction_envelope(regions, 32, float(log_potential))
    trials = []
    def objective(slope):
        check_time()
        intercept = float(np.max(caps-slope*np.arange(33)))
        groups = group_bounds(slope, weights)
        bound, _ = all_cohort_band(intercept, groups['total'], 18, 32, mark, float(log_potential))
        trials.append(dict(log2_b=float(slope), log2_a=intercept, failure_log2_screen=bound))
        return bound
    limits = (-64., 8.)
    search = minimize_scalar(objective, bounds=limits, method='bounded',
                             options=dict(maxiter=64, xatol=1e-9))
    objective(limits[0])
    objective(limits[1])
    best = min(trials, key=lambda row: row['failure_log2_screen'])
    a, b = best['log2_a'], best['log2_b']
    groups = group_bounds(b, weights)
    bound, rows = all_cohort_band(a, groups['total'], 18, 32, mark, float(log_potential))
    slack = a+b*np.arange(33)-caps
    if np.min(slack) < -1e-10:
        raise ArithmeticError('Complete regional affine envelope check failed')
    return dict(mark=str(mark), log2_potential=str(log_potential), failure_log2_screen=bound,
        entries=rows, group_log_bounds=groups, regional_log_contractions=caps.tolist(),
        log2_a_exact_float=str(F(a)), log2_b_exact_float=str(F(b)),
        envelope_contacts=[int(m) for m in range(33) if slack[m] < 1e-5],
        envelope_minimum_slack=float(np.min(slack)), affine_trials=trials,
        affine_evaluations=len(trials), affine_search_success=bool(search.success),
        winning_slope_hits_boundary=min(abs(b-edge) for edge in limits) < 1e-7)


def run(*, max_seconds=20, progress=None):
    if not 0 < max_seconds <= 20:
        raise ValueError('At most20 seconds for the entire fixed grid')
    result = source_scope()
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Bounded nine-pair diagnostic exhausted its budget')
    result.update(witnesses=[], raw_regional_banks=0, compute_budget_seconds=max_seconds)
    try:
        weights, pins = dense_arity_logs()
        result['exact_dense_moment_sha256'] = pins
        for mark in MARKS:
            check_time()
            regions = forced_nonzero_regions(float(mark), 32)
            result['raw_regional_banks'] += 1
            for logp in LOG_POTENTIALS:
                check_time()
                row = evaluate(regions, weights, mark, logp, check_time)
                result['witnesses'].append(row)
                if progress:
                    progress(dict(completed=len(result['witnesses']), total=9, mark=str(mark),
                        log2_potential=str(logp), failure_log2_screen=row['failure_log2_screen'],
                        envelope_contacts=row['envelope_contacts'], elapsed_seconds=perf_counter()-started))
        check_time()
        best = min(result['witnesses'], key=lambda row: row['failure_log2_screen'])
        selected = []
        for r in range(18, 33):
            bank = min(result['witnesses'], key=lambda row: row['entries'][r-18]['failure_log2_screen'])
            selected.append(dict(bank['entries'][r-18], mark=bank['mark'], log2_potential=bank['log2_potential']))
        result.update(status='completed', complete_requested_band=True, completed_witnesses=9,
            best_global_witness=dict(mark=best['mark'], log2_potential=best['log2_potential']),
            best_global_failure_log2_screen=best['failure_log2_screen'], entries=selected,
            failure_log2_screen=logsum(row['failure_log2_screen'] for row in selected),
            fixed_witness_selection_per_weight=True, witness_bank_penalty=1)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error), completed_witnesses=len(result['witnesses']))
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite fixed diagnostic evidence')
    result = run(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                     if key not in ('witnesses', 'entries', 'exact_dense_moment_sha256')}, allow_nan=False))

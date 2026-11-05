"""Bounded FLOAT all-cohort cut gate, using a pinned regional cap vector.

Keep c0<=1 separate from c_m<=C*B^m for m>0. Auxiliary iid allowed
coordinates replace C^union-weight. Holder across16 columns then sums all
fixed group types without resampling them. Fresh directed replay is needed.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize_scalar

from hd_two_sum_gf16_factorized import dense_arity_logs, logsum, times
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
INPUT_NAME = 'hd_two_sum_gf16_factorized_r18_32_grid9_float.json'
INPUT_SHA256 = '127dedbfab246d2a831913ec074b902e1b682b68fa2823b856d1b1fa1a7e47fb'
MARK, LOGP = F(99, 100), F(-64)
FLOOR_LOG_MARKERS = (0., 1., 4., 8., 16.)


def source_caps():
    raw = (DIRECTORY/INPUT_NAME).read_bytes()
    if sha256(raw).hexdigest() != INPUT_SHA256:
        raise ValueError('Pinned FLOAT regional diagnostic changed')
    data = json.loads(raw)
    if (data['construction'] != CONSTRUCTION or data['active_group_band'] != [18, 32]
            or data['directed'] is not False or data['raw_regional_banks'] != 3
            or data['regional_cap_amplitude_law'] != 'independent-uniform-nonzero'
            or data['regional_alpha'] != 1):
        raise ValueError('Wrong regional proposal provider or scope')
    selected = [row for row in data['witnesses']
                if row['mark'] == str(MARK) and row['log2_potential'] == str(LOGP)]
    if len(selected) != 1:
        raise ValueError('Exactly one frozen regional witness required')
    caps = np.array(selected[0]['regional_log_contractions'], float)
    if caps.shape != (33,) or np.any(~np.isfinite(caps)) or caps[0] > 0:
        raise ValueError('Complete finite caps and c0<=1 required; no empty surcharge may be discarded')
    return caps


def cut_logs(log_c, regions=16):
    """Exact binomial formula, evaluated in FLOAT logs; no zero flooring."""
    if not np.isfinite(log_c) or log_c > 0 or type(regions) is not int or regions < 1:
        raise ValueError('0<C<=1 and positive number of regions required')
    if log_c == 0:
        return np.array([-np.inf]*regions+[0.])
    log_not_c = float(np.log(-np.expm1(log_c*np.log(2.)))/np.log(2.))
    return np.array([log2(comb(regions, k))+k*log_c+times(regions-k, log_not_c)
                     for k in range(regions+1)])


def type_logs(log_b, log_x, log_y, dense_weights):
    """F_tau(k) for2 sparse types and14 dense arities, each once-only."""
    if (not all(np.isfinite(x) for x in (log_b, log_x, log_y)) or min(log_x, log_y) < 0
            or set(dense_weights) != set(range(3, 17))):
        raise ValueError('Finite slope, floor markers>=1 and all14 dense arities required')
    logs = np.empty((16, 17))
    qlog = log2(255)
    for k in range(17):
        full = 16*log_b if k == 16 else -np.inf
        single = np.logaddexp2(0., qlog+log_x+full)
        logs[0, k] = 4/16-qlog/2-9*log_x/16+single
        hole = log2(16 if k == 16 else 1)+15*log_b if k >= 15 else -np.inf
        double = qlog+2*log_x+np.logaddexp2(log2(239)+full, hole)
        pair = logsum((0., 1+qlog+log_x+full, double))
        logs[1, k] = log2(120)/16-qlog-18*log_x/16+pair
        for a in range(3, 17):
            weights = dense_weights[a]
            if len(weights) != 17 or np.any(np.isnan(weights)) or np.any(np.isposinf(weights)):
                raise ValueError('Complete finite-or-zero cardinality weights required')
            nonempty = logsum(log2(comb(k, w))+weights[w]+w*log_b for w in range(1, k+1))
            logs[a-1, k] = -9*log_y/16+np.logaddexp2(weights[0], log_y+nonempty)
    if np.any(~np.isfinite(logs)):
        raise ArithmeticError('Every type must retain a positive empty-column comparison')
    return logs


def bound_at(caps, weights, log_b, floor_log_marker, first=18, last=32):
    if (len(caps) != 33 or np.any(~np.isfinite(caps)) or caps[0] > 0
            or not 18 <= first <= last <= 32):
        raise ValueError('Complete cap vector and bounded band required')
    log_c = float(np.max(caps[1:]-log_b*np.arange(1, 33)))
    if log_c > 0:
        raise ValueError('The nonempty envelope requires C<=1')
    pi = cut_logs(log_c)
    types = type_logs(log_b, floor_log_marker, floor_log_marker, weights)
    group_moment = np.logaddexp2.reduce(16*types, axis=0)
    prefix = -LIMIT*log2(float(MARK))-log2(255)-float(LOGP)
    rows = []
    for r in range(first, last+1):
        terms = pi+(r/16)*group_moment
        dominant = int(np.argmax(terms))
        rows.append(dict(active_groups=r,
            failure_log2_screen=prefix+log2(comb(8192, r))+16*logsum(terms),
            dominant_allowed_size=dominant,
            dominant_type_at_peak=int(np.argmax(types[:, dominant]))))
    slack = log_c+log_b*np.arange(1, 33)-caps[1:]
    if np.min(slack) < -1e-10:
        raise ArithmeticError('Nonempty regional envelope failed its occupancy check')
    return dict(failure_log2_screen=logsum(row['failure_log2_screen'] for row in rows), entries=rows,
        log2_b_exact_float=str(F(float(log_b))), log2_c_exact_float=str(F(log_c)),
        floor_log_marker=str(F(float(floor_log_marker))), envelope_minimum_slack=float(np.min(slack)),
        envelope_contacts=[int(m) for m in range(1, 33) if slack[m-1] < 1e-5],
        type_order=['singleton', 'pair']+[f'dense-arity-{a}' for a in range(3, 17)],
        log2_type_weights=types.tolist(), log2_cut_probabilities=[None if np.isneginf(x) else float(x) for x in pi])


def run(*, max_seconds=20, progress=None):
    if not 0 < max_seconds <= 20:
        raise ValueError('At most20 seconds for the entire bounded cut gate')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Bounded cut gate exhausted its budget')
    caps = source_caps()
    weights, pins = dense_arity_logs()
    # This is the exact feasible endpoint for the computed FLOAT caps.
    lower = float(np.max(caps[1:]/np.arange(1, 33)))
    lower = float(np.nextafter(lower, np.inf))
    upper = 8.
    result = dict(construction=CONSTRUCTION, labels=list(LABELS), field_order=256,
        forbidden_weight=LIMIT, active_group_band=[18, 32], total_columns=16, constituent_dimension=8,
        all_source_supports=True, all_cohort_splits=True, all_constituent_row_sets=True,
        arbitrary_dense_group_count=True, source_counts='C(8192,r)', pair_swaps=True,
        fixed_rowsets_across_columns=True, fixed_types_across_columns=True,
        type_resampling_assumed=False, dense_arity_moments_separate=True,
        auxiliary_cut_identity=True, holder_columns=16, allowed_size_states=17,
        projective_normalizations=1, endpoint_potential_factors=1,
        dense_pergroup_floor=9, singleton_floor=9, pair_aggregate_floor=18,
        floor_method='positive-Chernoff-markers-at-least-one',
        floor_log_marker_schedules=list(FLOOR_LOG_MARKERS),
        mark=str(MARK), log2_potential=str(LOGP), z_p_fits=0,
        input_filename=INPUT_NAME, input_sha256=INPUT_SHA256,
        input_is_saved_float_proposal=True, fresh_directed_replay_required=True,
        raw_regional_banks_rebuilt=0, exact_dense_moment_sha256=pins,
        empty_cap_log2=float(caps[0]), empty_surcharge_discarded=False,
        log2_b_search_bounds=[lower, upper], maximum_slope_evaluations_per_schedule=42,
        fresh_selected_reconstructions_per_schedule=1, maximum_total_envelope_evaluations_per_schedule=43,
        complete_requested_band=False, directed=False, whole_code_certificate=False,
        compute_budget_seconds=max_seconds, schedules=[])
    try:
        if not lower < upper:
            raise ValueError('No authorized feasible slope interval')
        for marker in FLOOR_LOG_MARKERS:
            check_time()
            trials = []
            def objective(slope):
                check_time()
                row = bound_at(caps, weights, slope, marker)
                trials.append({key: row[key] for key in ('failure_log2_screen', 'log2_b_exact_float', 'log2_c_exact_float')})
                return row['failure_log2_screen']
            search = minimize_scalar(objective, bounds=(lower, upper), method='bounded',
                                     options=dict(maxiter=40, xatol=1e-8))
            objective(lower)
            objective(upper)
            best = min(trials, key=lambda row: row['failure_log2_screen'])
            chosen = bound_at(caps, weights, float(F(best['log2_b_exact_float'])), marker)
            chosen.update(slope_evaluations=len(trials), slope_trials=trials,
                slope_search_success=bool(search.success),
                slope_hits_boundary=min(abs(float(F(chosen['log2_b_exact_float']))-x) for x in (lower, upper)) < 1e-7)
            result['schedules'].append(chosen)
            if progress:
                progress(dict(completed=len(result['schedules']), total=len(FLOOR_LOG_MARKERS),
                    floor_log_marker=marker, failure_log2_screen=chosen['failure_log2_screen'],
                    envelope_contacts=chosen['envelope_contacts'], elapsed_seconds=perf_counter()-started))
        check_time()
        selected = []
        for r in range(18, 33):
            bank = min(result['schedules'], key=lambda row: row['entries'][r-18]['failure_log2_screen'])
            selected.append(dict(bank['entries'][r-18], floor_log_marker=bank['floor_log_marker'],
                log2_b_exact_float=bank['log2_b_exact_float'], log2_c_exact_float=bank['log2_c_exact_float']))
        result.update(status='completed', complete_requested_band=True, entries=selected,
            failure_log2_screen=logsum(row['failure_log2_screen'] for row in selected),
            witness_selection_per_weight=True, extra_witness_bank_factor=False)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded cut artifact')
    result = run(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                     if key not in ('schedules', 'entries', 'exact_dense_moment_sha256')}, allow_nan=False))

"""Saved F1 quiet comparison: row sets and word coefficients before alpha.

Fix the source G, but do not union over its arity8 row-set assignments yet.
The saved two-atom curve majorizes each fixed assignment's unconditional
quiet/all-regional-00 first moment. Sum assignments, keep the whole-word
input-weight tail, and only then apply the scalar fractional inequality.
The source count is still outside alpha. The direct first-moment endpoint
alpha=1 uses the SAME unpowered provider and the SAME scalar comparison.

This covers one positive component of an enumerator majorant, not the
complete failure event, a probability lower bound, or a distance theorem.
"""
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_late_alpha_quiet import (ALPHA, ARITY, COLUMNS, DIMENSION,
    GROUPS, MARK, MINIMUM, log_power)
from hd_two_sum_late_alpha_moment import (checked_curve,
    quiet_mixture_coefficients, two_atom_envelope)
from hd_two_sum_late_alpha_source_f1 import load_source_count
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
ARTIFACT = 'hd_two_sum_late_alpha_source_f1_gate_float.json'
ARTIFACT_SHA = '862d9f5dd08c304d70f282563bfb5130e5ce35423e7153d3b68721337838b045'
SOURCE_PINS = {
    'hd_two_sum_late_alpha_source_f1.py': '77d3f7da4f7439766d0bb18b53cebcd0960697094c0e30361b064884b5ee1987',
    'hd_two_sum_late_alpha_quiet.py': '8a753036a9be44759206b5beb3baa044f28d80b265a9304da81bdde6e8f11ae3',
    'hd_two_sum_late_alpha_moment.py': 'fe9a8936e99f82a4aec31f903bc55ad6ebfec7fa8a67fdb30f6e37bc4f2ca8e6',
}


def checked_artifact(data, source_count):
    required = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=1,
        mark=str(MARK), alpha=str(ALPHA), columns=COLUMNS,
        constituent_dimension=DIMENSION, fresh_regional_alpha=1,
        complete_requested_controls=True, all_sources=False,
        all_sources_in_declared_stratum=True, all_fixed_arity_rowset_vectors=True,
        rowsets_resampled=False, directed=False, whole_code_certificate=False,
        source_count=str(source_count), rowset_count=f'{comb(16, ARITY)}^{GROUPS}')
    if any(data.get(key) != value for key, value in required.items()):
        raise ValueError('Exact frozen F1 source and fixed-arity quiet scope required')
    if data.get('source_sha256') != SOURCE_PINS['hd_two_sum_late_alpha_source_f1.py']:
        raise ValueError('Frozen source identity changed')
    curve = checked_curve([-np.inf if x is None else x
                           for x in data['line_envelope']['log2_curves']], GROUPS)
    if data.get('two_atom_envelope', {}).get('secant_endpoints') != [5, 18]:
        raise ValueError('Only the frozen regional moment secant is allowed')
    return curve


def read_source():
    raw = (DIRECTORY/ARTIFACT).read_bytes()
    if sha256(raw).hexdigest() != ARTIFACT_SHA:
        raise ValueError('Frozen F1 artifact changed')
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Frozen helper changed: '+name)
    count = load_source_count()
    data = json.loads(raw)
    return data, checked_artifact(data, count), count


def order_comparison(raw, *, columns, floor, alpha, log_source_count,
                     log_rowset_count, log_global_marker, log_word_normalizer,
                     check_time=lambda: None):
    """Compare the two orders and alpha=1 using identical raw coefficients.

    global_marker is -LIMIT*log2(z)-log2(255). The word normalizer is
    8*a*r*log2(255). Row-set and source counts are distinct arguments so
    that no source Jensen or extra row-set choice is silently introduced.
    """
    raw = np.asarray(raw, dtype=float)
    if (raw.ndim != 1 or not len(raw) or np.any(np.isnan(raw))
            or np.any(np.isposinf(raw)) or type(columns) is not int or columns < 1
            or type(floor) is not int or not 0 <= floor <= columns*(len(raw)-1)
            or not 0 < alpha <= 1 or any(not np.isfinite(value) for value in
                (log_source_count, log_rowset_count, log_global_marker, log_word_normalizer))):
        raise ValueError('Finite positive comparison data and a reachable tail required')
    check_time()
    whole = log_power(raw, columns, check_time)
    raw_tail = float(np.logaddexp2.reduce(whole[floor:]))
    expanded = log_power(alpha*raw-alpha*log_word_normalizer/columns,
                         columns, check_time)
    expanded_tail = float(np.logaddexp2.reduce(expanded[floor:]))
    log_mean = log_rowset_count+log_global_marker-log_word_normalizer+raw_tail
    old = log_source_count+log_rowset_count+alpha*log_global_marker+expanded_tail
    degree_late = log_source_count+log_rowset_count+alpha*(
        log_global_marker-log_word_normalizer+raw_tail)
    return dict(raw_word_tail_log2=raw_tail,
        fixed_source_all_rowsets_mean_log2=log_mean,
        coefficientwise_alpha_rowset_union_log2=old,
        whole_word_alpha_rowset_union_log2=degree_late,
        whole_word_and_rowset_alpha_log2=log_source_count+alpha*log_mean,
        direct_first_moment_alpha1_log2=log_source_count+log_mean,
        degree_order_gain_bits=old-degree_late,
        rowset_order_gain_bits=(1-alpha)*log_rowset_count,
        log2_source_count=log_source_count, log2_rowset_count=log_rowset_count,
        log2_global_marker=log_global_marker, log2_word_normalizer=log_word_normalizer)


def run(*, max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit positive budget at most30 seconds required')
    source, curve, count = read_source()
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Saved whole-word ordering control exhausted its budget')
    result = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=1,
        mark=str(MARK), old_alpha=str(ALPHA), endpoint_alpha='1',
        scope='Fmax1 sources; quiet moments and all regional boundary states zero only',
        source_count=str(count), rowset_count=f'{comb(16, ARITY)}^{GROUPS}',
        all_sources=False, all_sources_in_declared_stratum=True,
        all_fixed_arity_rowset_vectors=True, rowsets_fixed_across_columns=True,
        source_count_outside_alpha=True, rowset_sum_inside_alpha=True,
        whole_word_tail_before_alpha=True, aggregate_input_weight_floor=MINIMUM*ARITY*GROUPS,
        individual_row_floors=False, all_rows_seen_constraint=False,
        raw_regional_banks_built=0, new_moment_candidates=0, parameter_fits=0,
        directed=False, whole_code_certificate=False, full_quiet_event_bounded=False,
        failure_probability_lower_bound=False, artifact=ARTIFACT, artifact_sha256=ARTIFACT_SHA,
        reviewed_source_sha256=SOURCE_PINS, compute_budget_seconds=max_seconds,
        complete_requested_control=False)
    try:
        moments, envelope = two_atom_envelope(curve, 5, 18)
        if any(abs(envelope[name]-source['two_atom_envelope'][name]) > 1e-12
               for name in ('log2_nonzero_atom', 'log2_probability')):
            raise ArithmeticError('Frozen regional envelope did not reproduce')
        raw = quiet_mixture_coefficients(GROUPS, ARITY, moments, check_time=check_time)
        comparison = order_comparison(raw, columns=COLUMNS, floor=MINIMUM*ARITY*GROUPS,
            alpha=float(ALPHA), log_source_count=log2(count),
            log_rowset_count=GROUPS*log2(comb(16, ARITY)),
            log_global_marker=-LIMIT*log2(float(MARK))-log2(255),
            log_word_normalizer=DIMENSION*ARITY*GROUPS*log2(255), check_time=check_time)
        if abs(comparison['coefficientwise_alpha_rowset_union_log2']-
               source['moment_control']['majorant_term_log2']) > 1e-8:
            raise ArithmeticError('Frozen comparison did not reproduce')
        result['two_atom_envelope'] = envelope
        result['comparison'] = comparison
        check_time()
        result['complete_requested_control'] = True
        if progress:
            progress(comparison)
    except TimeoutError as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--max-seconds', default=30, type=float)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen control')
    answer = run(max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(complete=answer['complete_requested_control'], elapsed_seconds=answer['elapsed_seconds'])))

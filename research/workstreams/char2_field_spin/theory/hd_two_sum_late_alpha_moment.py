"""Two saved-curve quiet controls, no regional provider or fitted bank.

The first scope has ONE common row set A across the18 source groups and
uses the saved raw R1_00 curve. The second covers all fixed arity8 row-set
vectors, using one two-atom positive moment envelope of that curve.
Neither includes nonquiet outer moments or live regional boundary states.

Z2(w)<=255**(w-2) for w>0 gives a positive group mixture: eta=255**-2,
(1-eta)+eta*prod_i(1+255*x_i). For h selected groups and common A, averaging
regional activity gives F_h(X)**8. For arbitrary fixed A_g, the regional
incidences n_i obey0<=n_i<=h and sum n_i=8h. A positive moment envelope
makes F_n coefficientwise strongly log-convex, hence prod_i F_n_i<=F_h**8.
The row-set factor differs between the two scopes and stays outside alpha.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log, log1p, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_late_alpha_quiet import (ALPHA, ARITY, COLUMNS, DIMENSION,
    GROUPS, MARK, MINIMUM, column_coefficients, log_power)
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
SOURCE = 'hd_two_sum_late_alpha_quiet_gate_float.json'
SOURCE_SHA = '1f9bcb13d4a5120a130e3b23dfad264ffe25a3fc42f9b5c2dc4e8488732d4521'
QUIET_HELPER_SHA = '8a753036a9be44759206b5beb3baa044f28d80b265a9304da81bdde6e8f11ae3'


def checked_curve(values, maximum):
    curve = np.asarray(values, dtype=float)
    if (type(maximum) is not int or maximum < 1 or curve.shape != (maximum+1,)
            or np.any(np.isnan(curve)) or np.any(np.isposinf(curve))
            or curve[0] != 0.):
        raise ValueError('Complete log curve with R0=1 required')
    return curve


def marked_regional_polynomial(count, curve, q=256):
    """F_n(X)=sum_m C(n,m)(q-1)^m R_m X^m."""
    if type(count) is not int or not 0 <= count < len(curve) or q < 4:
        raise ValueError('A covered incidence count and valid q are required')
    return np.array([log2(comb(count, m))+m*log2(q-1)+curve[m]
                     for m in range(count+1)])


def quiet_mixture_coefficients(groups, arity, curve, *, q=256, check_time=lambda: None):
    """Positive raw column comparison, before denominator and alpha.

    This function alone is COMMON-A only. The caller may claim arbitrary
    row-set scope only when curve consists of moments of one probability
    measure; mere numerical agreement with a convex sequence is not used.
    """
    if (type(groups) is not int or groups < 1 or type(arity) is not int
            or not 2 <= arity <= min(q, 16) or type(q) is not int or q < 4
            or q & (q-1)):
        raise ValueError('Positive groups and admissible arity/field required')
    curve = checked_curve(curve, groups)
    eta = 1/(q-1)**2
    log_eta, log_other = -2*log2(q-1), log1p(-eta)/log(2)
    result = np.full(groups*arity+1, -np.inf)
    for count in range(groups+1):
        check_time()
        log_weight = log2(comb(groups, count))+count*log_eta+(groups-count)*log_other
        factor = log_power(marked_regional_polynomial(count, curve, q), arity, check_time)
        np.logaddexp2(result[:len(factor)], factor+log_weight, out=result[:len(factor)])
    return result


def two_atom_envelope(curve, left=5, right=18):
    """ONE frozen secant, covered by p=max_m R_m/B^m. No alternatives fit.

    The represented probability law is L=B with probability p and L=0
    otherwise. Thus E L^0=1 and E L^m=p B^m for all m>0.
    """
    curve = checked_curve(curve, len(curve)-1)
    if (type(left) is not int or type(right) is not int
            or not 1 <= left < right < len(curve)
            or not np.isfinite(curve[left]) or not np.isfinite(curve[right])):
        raise ValueError('Two fixed finite positive-occupancy secant endpoints required')
    slope = float((curve[right]-curve[left])/(right-left))
    log_probability = float(np.max(curve[1:]-slope*np.arange(1, len(curve))))
    if not np.isfinite(log_probability) or log_probability > 0:
        raise ArithmeticError('The prescribed two-atom envelope is not a probability law')
    moments = log_probability+slope*np.arange(len(curve))
    moments[0] = 0.
    slack = moments-curve
    if np.min(slack) < -1e-9:
        raise ArithmeticError('The moment envelope does not cover every saved occupancy')
    return moments, dict(secant_endpoints=[left, right], log2_nonzero_atom=slope,
        log2_nonzero_atom_as_exact_float=str(F(slope)),
        log2_probability=log_probability, log2_probability_as_exact_float=str(F(log_probability)),
        moment_order_zero=1, contacts=np.flatnonzero(slack < 1e-9).tolist(),
        minimum_float_slack=float(np.min(slack)), probability_at_most_one=True)


def word_tail(raw, *, alpha=float(ALPHA), groups=GROUPS, arity=ARITY,
              columns=COLUMNS, dimension=DIMENSION, minimum=MINIMUM,
              q=256, check_time=lambda: None):
    normalizer = alpha*dimension*arity*groups/columns*log2(q-1)
    powered = log_power(alpha*np.asarray(raw)-normalizer, columns, check_time)
    return float(np.logaddexp2.reduce(powered[minimum*arity*groups:]))


def read_source(directory=DIRECTORY):
    directory = Path(directory)
    raw = (directory/SOURCE).read_bytes()
    if sha256(raw).hexdigest() != SOURCE_SHA:
        raise ValueError('Frozen fresh-R1 control changed')
    if sha256((directory/'hd_two_sum_late_alpha_quiet.py').read_bytes()).hexdigest() != QUIET_HELPER_SHA:
        raise ValueError('Frozen late-alpha helper changed')
    data = json.loads(raw)
    required = dict(active_groups=GROUPS, arity=ARITY, mark=str(MARK), alpha=str(ALPHA),
        columns=COLUMNS, constituent_dimension=DIMENSION, fresh_regional_alpha=1,
        raw_regional_banks=1, complete_requested_control=True, directed=False,
        whole_code_certificate=False, source_sha256=QUIET_HELPER_SHA)
    if any(data.get(key) != value for key, value in required.items()):
        raise ValueError('Only the exact frozen r18/a8 quiet control is allowed')
    curve = [float('-inf') if value is None else value
             for value in data['regional_envelope']['log2_curves']]
    return data, checked_curve(curve, GROUPS)


def run(*, max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit positive budget at most30 seconds required')
    source, curve = read_source()
    start = perf_counter()
    def check_time():
        if perf_counter()-start >= max_seconds:
            raise TimeoutError('Saved-curve moment gate exhausted its budget')
    prefix = (log2(comb(8192, GROUPS))+float(ALPHA)*(-LIMIT*log2(float(MARK))-log2(255)))
    old_late, _ = column_coefficients(GROUPS, ARITY,
        source['regional_envelope']['log2_slope'], float(ALPHA), check_time=check_time)
    old_word = log_power(old_late, COLUMNS, check_time)
    old_tail = float(np.logaddexp2.reduce(old_word[MINIMUM*ARITY*GROUPS:]))
    if abs(old_tail-source['comparison']['late_column_log2']) > 1e-8:
        raise ArithmeticError('Frozen baseline reproduction failed')
    result = dict(source=SOURCE, source_sha256=SOURCE_SHA, frozen_helper_sha256=QUIET_HELPER_SHA,
        active_groups=GROUPS, arity=ARITY, mark=str(MARK), alpha=str(ALPHA),
        scope='zero moments and zero regional boundary states only',
        raw_regional_banks_built=0, fitted_marks=0, requested_controls=2,
        moment_secant_candidates=1, individual_row_floors=False,
        aggregate_input_weight_floor=MINIMUM*ARITY*GROUPS,
        auxiliary_group_mixture='eta=1/255^2; new group labels are never sampled',
        synthetic_weights_one_and_two_admitted=True,
        row_sets_fixed_across_columns=True, source_count_once=True,
        source_and_rowset_counts_outside_alpha=True, projective_normalizations=1,
        directed=False, whole_code_certificate=False, full_quiet_event_bounded=False,
        failure_probability_lower_bound=False, compute_budget_seconds=max_seconds,
        complete_requested_controls=False, controls=[])
    try:
        moments, envelope = two_atom_envelope(curve)
        result['two_atom_envelope'] = envelope
        for label, values, row_exponent in (
            ('common-A-raw-curve', curve, 1),
            ('all-fixed-A-two-atom', moments, GROUPS)):
            check_time()
            raw = quiet_mixture_coefficients(GROUPS, ARITY, values, check_time=check_time)
            tail = word_tail(raw, check_time=check_time)
            class_prefix = prefix+row_exponent*log2(comb(16, ARITY))
            item = dict(scope=label, rowset_count=f'{comb(16, ARITY)}^{row_exponent}',
                raw_curve_used=label == 'common-A-raw-curve',
                arbitrary_fixed_rowset_vectors=label == 'all-fixed-A-two-atom',
                coefficientwise_moment_balancing=label == 'all-fixed-A-two-atom',
                raw_column_coefficient_degree81_log2=float(raw[81]),
                normalized_column_coefficient_degree81_log2=float(raw[81]-72*log2(255)),
                normalized_word_log2=tail, log2_prefix=class_prefix,
                majorant_term_log2=class_prefix+tail,
                old_matched_scope_majorant_log2=class_prefix+old_tail,
                new_minus_old_bits=tail-old_tail)
            result['controls'].append(item)
            if progress:
                progress(dict(item, elapsed_seconds=perf_counter()-start))
        check_time()
        result['complete_requested_controls'] = True
    except TimeoutError as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-start
    result['implementation_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--max-seconds', default=30, type=float)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen control')
    result = run(max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(complete=result['complete_requested_controls'],
                         elapsed_seconds=result['elapsed_seconds'])))

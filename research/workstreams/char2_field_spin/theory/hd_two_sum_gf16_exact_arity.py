"""FLOAT exact-total-arity norms for the fixed-type state-path bound.

For each local state-path class n, form P_n(u)=sum_a F_a(n)^16*u^a.
Restricting [u^A]P_n(u)^r BEFORE its16th root is exactly the type-vector
restriction in the path Minkowski step. No total-arity Chernoff marker
or cohort enumeration is needed. The existing regional and floor
envelopes remain unchanged. This module never builds a regional bank.
"""
from fractions import Fraction as F
from math import comb, log2

import numpy as np

import hd_two_sum_gf16_state_envelope as state
from hd_two_sum_gf16_sparse import CONSTRUCTION
from hd_xor_multigroup import LIMIT


MAX_GROUPS = 32


def checked_requests(type_logs, intervals, columns):
    """Intervals are inclusive and keyed by the exact active-group count."""
    values = np.asarray(type_logs, dtype=float)
    if (values.ndim != 2 or not 1 <= values.shape[1] <= 16
            or values.shape[0] < 1 or np.any(np.isnan(values)) or np.any(np.isposinf(values))
            or type(columns) is not int or not 1 <= columns <= 16
            or not isinstance(intervals, dict) or not intervals):
        raise ValueError('Finite-or-zero path/type logs and explicit nonempty intervals required')
    result = {}
    for groups, interval in intervals.items():
        if (type(groups) is not int or not 1 <= groups <= MAX_GROUPS
                or not isinstance(interval, (tuple, list)) or len(interval) != 2
                or any(type(bound) is not int for bound in interval)):
            raise ValueError('At most32 groups and two exact integer arity endpoints required')
        lower, upper = interval
        if not 0 <= lower <= upper <= values.shape[1]*groups:
            raise ValueError('Arity interval must stay in the natural0..types*groups domain')
        result[groups] = lower, upper
    return values, result


def restricted_type_norms(type_logs, intervals, *, columns=16, check_time=None):
    """Return log2 of each restricted ell-columns norm, batched over paths.

    Multiplication uses only positive additions in logarithmic arithmetic.
    Truncation discards degrees above every requested endpoint. Degrees
    never decrease, so those coefficients cannot affect a requested sum.
    All requested group counts share the same sequence of products.
    """
    values, intervals = checked_requests(type_logs, intervals, columns)
    paths, types = values.shape
    maximum = max(intervals)
    cap = max(upper for _, upper in intervals.values())
    coefficients = np.zeros((paths, 1))
    scaled = columns*values
    result = {}
    for groups in range(1, maximum+1):
        if check_time is not None:
            check_time()
        degree = min(types*groups, cap)
        following = np.full((paths, degree+1), -np.inf)
        for arity in range(1, min(types, degree)+1):
            width = min(coefficients.shape[1], degree+1-arity)
            following[:, arity:arity+width] = np.logaddexp2(
                following[:, arity:arity+width],
                coefficients[:, :width]+scaled[:, arity-1, None])
        coefficients = following
        if groups in intervals:
            lower, upper = intervals[groups]
            result[groups] = state.logsum(coefficients[:, lower:upper+1], axis=1)/columns
    return result


def kernel_from_restricted_norms(classes, path_a, norms):
    """Keep exact path multiplicities and ordered endpoint labels."""
    path_a, norms = np.asarray(path_a), np.asarray(norms)
    if (path_a.shape != (len(classes),) or norms.shape != path_a.shape
            or np.any(~np.isfinite(path_a)) or np.any(np.isnan(norms))
            or np.any(np.isposinf(norms))):
        raise ValueError('Complete path intercepts and finite-or-zero restricted norms required')
    terms = path_a+norms+np.log2([row[3] for row in classes])
    result = np.full((2, 2), -np.inf)
    for i in range(2):
        for j in range(2):
            indices = [index for index, row in enumerate(classes) if row[:2] == (i, j)]
            if indices:
                result[i, j] = state.logsum(terms[indices])
    return result


def exact_low_kernel(prepared, groups, maximum_arity, *, columns=16, check_time=None):
    """Use the low helper's unchanged prepared path/type arrays, without t."""
    norm = restricted_type_norms(prepared['typed'], {groups: (0, maximum_arity)},
        columns=columns, check_time=check_time)[groups]
    return kernel_from_restricted_norms(prepared['classes'], prepared['path_a'], norm)


def evaluate_prepared(prepared, groups, cutoff, mark, alpha, *, check_time=None):
    """One explicit A_tot<cutoff class; no high-side or band coverage claim."""
    mark, alpha = F(mark), F(alpha)
    if (type(groups) is not int or not 1 <= groups <= MAX_GROUPS
            or type(cutoff) is not int or not 1 <= cutoff <= 16*groups+1
            or not 0 < mark < 1 or not 0 < alpha <= 1
            or np.shape(prepared['typed']) != (len(prepared['classes']), 16)):
        raise ValueError('Explicit G16 group/cutoff scope and valid frozen witness required')
    kernel = exact_low_kernel(prepared, groups, cutoff-1, check_time=check_time)
    power = state.matrix_power_log(kernel, 16)
    prefix = float(alpha)*(-LIMIT*log2(float(mark))-log2(255))
    bound = prefix+log2(comb(8192, groups))+float(state.logsum(power[0]))
    return dict(construction=CONSTRUCTION, active_groups=groups, high_arity_cutoff=cutoff,
        low_total_active_constituents_max=cutoff-1, mark=str(mark), alpha=str(alpha),
        source_counts='C(8192,r)', source_supports=comb(8192, groups),
        all_cohorts_within_low_total_arity=True, all_constituent_row_sets_within_low_total_arity=True,
        exact_total_arity_coefficient_restriction=True, total_arity_chernoff_marker_used=False,
        fixed_types_across_columns=True, type_norm_root=16, outer_columns=16,
        regional_banks_rebuilt=0, fixed_regional_and_floor_envelopes=True,
        initial_state=0, terminal_vector=[1, 1], projective_normalizations=1,
        projective_inside_alpha=True, endpoint_potential_factors=0,
        failure_log2_screen=float(bound), directed=False, high_side_included=False,
        continuous_weight_band_covered=False, whole_code_certificate=False,
        requires_fresh_directed_replay=True)

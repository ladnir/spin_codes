"""Positive whole-group nonzero and aggregate RS-column-floor bounds.

Fix every group type and row set across all outer columns. Split each
column's mask assignment by its labelled group activity bits before the
fixed-type Holder step. If f[a,e,p] is the mass for type a, activity e and
chronological local path p, define U[e,p]=(sum_a f[a,e,p]^columns)^(1/columns).
Holder for each fixed activity-bit array leaves a product of these U's.

We then sum only arrays where every group is seen, or only arrays with
at least minimum*groups active group-column incidences. Both are positive
relaxations of each group's RS floor. Their minimum is valid; neither
asserts the two constraints jointly. No group type is actually resampled.
The branch split can increase Holder slack, so improvement is not assumed.
"""
from fractions import Fraction as F
from math import comb, log2

import numpy as np

import hd_two_sum_gf16_state_envelope as state
from hd_xor_multigroup import LIMIT
from hd_two_sum_gf16_sparse import CONSTRUCTION


def scaled_log(value, exponent):
    """Positive integer powers in logs, including the exact zero exponent."""
    return np.zeros_like(value) if exponent == 0 else exponent*value


def prepare_branches(prepared, weights, *, columns=16):
    """Reuse the path/cut elementary coefficients; never rebuild a provider.

    weights[a][w] is a per-labelled-mask log weight. The elementary
    coefficient counts masks of weight w on the path's allowed positions.
    The same type weights must be fixed across all paths and columns.
    """
    classes = prepared.get('classes', ())
    elementary = np.asarray(prepared.get('elementary'), dtype=float)
    path_a = np.asarray(prepared.get('path_a'), dtype=float)
    if (type(columns) is not int or not 1 <= columns <= 16 or not classes
            or elementary.ndim != 2 or elementary.shape[0] != len(classes)
            or path_a.shape != (len(classes),) or np.any(~np.isfinite(path_a))
            or np.any(np.isnan(elementary)) or np.any(np.isposinf(elementary))
            or not isinstance(weights, dict) or set(weights) != set(range(1, len(weights)+1))
            or not 1 <= len(weights) <= 16):
        raise ValueError('Complete path/cut coefficients and consecutive fixed type laws required')
    laws = np.stack([np.asarray(weights[a], dtype=float) for a in range(1, len(weights)+1)])
    if (laws.shape != (len(weights), elementary.shape[1])
            or np.any(np.isnan(laws)) or np.any(np.isposinf(laws))):
        raise ValueError('Finite-or-zero per-mask laws at every cardinality required')
    empty = elementary[:, :1]+laws[:, 0][None, :]
    nonempty = np.stack([state.logsum(elementary[:, 1:]+law[1:], axis=1) for law in laws], axis=1)
    branches = np.stack((empty, nonempty), axis=1)
    u = state.logsum(columns*branches, axis=2)/columns
    return dict(classes=classes, path_a=path_a, log_u=u,
        type_branch_logs=branches, columns=columns, fixed_type_weights_across_columns=True)


def column_moments(branches, groups):
    """M_l is the ordered endpoint matrix for l specified active groups.

    Path and cut multiplicities stay outside Holder. C(groups,l) is NOT
    included: different consumers insert their own exact labelled counts.
    """
    classes = branches['classes']
    path_a, u = np.asarray(branches['path_a']), np.asarray(branches['log_u'])
    if (type(groups) is not int or not 1 <= groups <= 32
            or path_a.shape != (len(classes),) or u.shape != (len(classes), 2)
            or np.any(~np.isfinite(path_a)) or np.any(np.isnan(u)) or np.any(np.isposinf(u))
            or any(row[0] not in (0, 1) or row[1] not in (0, 1)
                   or type(row[3]) is not int or row[3] < 1 for row in classes)):
        raise ValueError('Positive path multiplicities and bounded source group count required')
    base = path_a+np.log2([row[3] for row in classes])
    positions = {(i, j): [k for k, row in enumerate(classes) if row[:2] == (i, j)]
                 for i in range(2) for j in range(2)}
    result = np.full((groups+1, 2, 2), -np.inf)
    for active in range(groups+1):
        terms = base+scaled_log(u[:, 0], groups-active)+scaled_log(u[:, 1], active)
        for (i, j), indices in positions.items():
            if indices:
                result[active, i, j] = state.logsum(terms[indices])
    return result


def checked_moments(moments, columns):
    values = np.asarray(moments, dtype=float)
    if (values.ndim != 3 or values.shape[1:] != (2, 2) or not 2 <= values.shape[0] <= 33
            or np.any(np.isnan(values)) or np.any(np.isposinf(values))
            or type(columns) is not int or not 1 <= columns <= 16):
        raise ValueError('Finite-or-zero moments and at most16 ordered columns required')
    return values, values.shape[0]-1


def seen_all_log(moments, *, columns=16, check_time=None):
    """Track c seen groups exactly, starting from (state0,c0).

    T[c,c+t] = C(r-c,t) sum_j C(c,j) M_(t+j).
    Endpoint matrix multiplication is chronological; final c must equal r.
    """
    moments, groups = checked_moments(moments, columns)
    transition = np.full((2*(groups+1), 2*(groups+1)), -np.inf)
    for seen in range(groups+1):
        if check_time:
            check_time()
        for new in range(groups-seen+1):
            entries = moments[new:new+seen+1]+np.log2([comb(seen, j) for j in range(seen+1)])[:, None, None]
            transition[2*seen:2*seen+2, 2*(seen+new):2*(seen+new)+2] = (
                log2(comb(groups-seen, new))+state.logsum(entries, axis=0))
    current = np.full(2*(groups+1), -np.inf)
    current[0] = 0.
    for _ in range(columns):
        if check_time:
            check_time()
        current = state.logsum(current[:, None]+transition, axis=0)
    return float(state.logsum(current[2*groups:2*groups+2]))


def aggregate_floor_log(moments, *, columns=16, minimum=9, check_time=None):
    """Exact positive coefficient tail of the ordered matrix polynomial.

    K(Y)=sum_l C(r,l) M_l Y^l; retain degrees >=minimum*r in K(Y)^columns.
    This is the aggregate floor only, not an individual floor for each group.
    """
    moments, groups = checked_moments(moments, columns)
    if type(minimum) is not int or not 0 <= minimum <= columns:
        raise ValueError('A nonnegative per-group necessary floor at most columns required')
    coefficients = moments+np.log2([comb(groups, active) for active in range(groups+1)])[:, None, None]
    current = np.full((2, 1), -np.inf)
    current[0, 0] = 0.
    for _ in range(columns):
        if check_time:
            check_time()
        width = current.shape[1]
        following = np.full((2, width+groups), -np.inf)
        for active in range(groups+1):
            for i in range(2):
                for j in range(2):
                    following[j, active:active+width] = np.logaddexp2(
                        following[j, active:active+width], current[i]+coefficients[active, i, j])
        current = following
    return float(state.logsum(current[:, minimum*groups:]))


def evaluate(prepared, weights, groups, mark, alpha, *, source_count=None, check_time=None):
    """One all-cohort FLOAT class; caller may provide an exact source stratum count."""
    mark, alpha = F(mark), F(alpha)
    if (type(groups) is not int or not 1 <= groups <= 32 or not 0 < mark < 1 or not 0 < alpha <= 1
            or prepared.get('regions') != 16 or set(weights) != set(range(1, 17))):
        raise ValueError('Full G16 fixed type laws, valid witnesses and1..32 active groups required')
    count = comb(8192, groups) if source_count is None else source_count
    if type(count) is not int or not 1 <= count <= comb(8192, groups):
        raise ValueError('Positive indexed-source count within the requested source weight required')
    branches = prepare_branches(prepared, weights)
    moments = column_moments(branches, groups)
    seen = seen_all_log(moments, check_time=check_time)
    aggregate = aggregate_floor_log(moments, check_time=check_time)
    prefix = float(alpha)*(-LIMIT*log2(float(mark))-log2(255))+log2(count)
    return dict(construction=CONSTRUCTION, active_groups=groups, mark=str(mark), alpha=str(alpha),
        indexed_source_supports=count, source_count_outside_alpha=True,
        seen_all_failure_log2_screen=prefix+seen,
        aggregate_9r_failure_log2_screen=prefix+aggregate,
        failure_log2_screen=prefix+min(seen, aggregate),
        exact_whole_group_nonzero_bound=True, exact_aggregate_nonempty_incidence_floor=9*groups,
        individual_nine_column_floors_enforced=False, two_constraints_enforced_jointly=False,
        fixed_types_before_activity_pattern_holder=True, fixed_rowsets_across_columns=True,
        activity_branch_holder_exponent=16, activity_branches_resampled=False,
        state_endpoints_ordered=True, initial_state=0, terminal_vector=[1, 1],
        columns=16, seen_states=groups+1, aggregate_maximum_degree=16*groups,
        projective_normalizations=1, endpoint_potential_factors=0,
        regional_banks_rebuilt=0, directed=False, exact_total_arity_restriction=False,
        whole_code_certificate=False, requires_fresh_directed_replay=True)

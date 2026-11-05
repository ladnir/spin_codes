"""Positive off-diagonal support constraints for the fixed-type bound.

For a fixed state path, every01/10 coordinate must be covered by at least
one group mask. After the common symmetric envelope/type replacement,
the numbers covered in these two classes suffice. T_a is the exact
positive coverage transition for one FIXED type a. For p=16,

 ||(T_tau1 ... T_taur)[empty,full]||_{ell_p(tau)}
 <= (U**r)[empty,full],  U[c,d]=||T_a[c,d]||_{ell_p(a)}.

Minkowski over coverage histories proves this inequality. Type vectors
stay fixed across the outer columns; no type is sampled or refreshed.
The original unrestricted path norm is another valid bound, so the
smaller bound is used before adding chronological path contributions.
All arithmetic here is FLOAT. No inclusion-exclusion, regional provider,
fit, or automatic control run is included.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, log2

import numpy as np

import hd_two_sum_gf16_state_emptycut as emptycut
import hd_two_sum_gf16_state_envelope as state
import hd_two_sum_gf16_exact_arity as exact


@lru_cache(81, typed=True)
def coverage_geometry(first, second):
    """Positive binomial map from required-subset weights to transitions.

    Rows of coefficients index feasible coverage transitions. Columns
    index how many01/10 coordinates this group's mask selects. Old and
    newly covered coordinates are distinguished, then counted exactly.
    """
    if (any(type(n) is not int or n < 0 for n in (first, second))
            or first+second > 16):
        raise ValueError('At most16 required coordinates in two classes')
    states = tuple((i, j) for i in range(first+1) for j in range(second+1))
    sources, targets, rows = [], [], []
    for source, (covered1, covered2) in enumerate(states):
        for new1 in range(first-covered1+1):
            for new2 in range(second-covered2+1):
                target = (covered1+new1)*(second+1)+covered2+new2
                factor = comb(first-covered1, new1)*comb(second-covered2, new2)
                row = np.zeros((first+1, second+1))
                for old1 in range(covered1+1):
                    for old2 in range(covered2+1):
                        row[old1+new1, old2+new2] = factor*comb(covered1, old1)*comb(covered2, old2)
                sources.append(source)
                targets.append(target)
                rows.append(row.ravel())
    sources, targets = np.asarray(sources), np.asarray(targets)
    coefficients = np.asarray(rows)
    for value in (sources, targets, coefficients):
        value.setflags(write=False)
    return states, sources, targets, coefficients


def type_transition_logs(weights, required, required_marks, optional_marks):
    """Exact positive cover transitions, up to FLOAT evaluation error.

    weights[a,w] is the mass of ONE labelled mask of weight w for type a.
    It already includes fixed row-set and floor factors. Selecting an
    optional-coordinate subset uses its elementary polynomial once.
    """
    weights = np.asarray(weights, dtype=float)
    required_marks, optional_marks = np.asarray(required_marks, dtype=float), np.asarray(optional_marks, dtype=float)
    if (len(required) != 2 or weights.ndim != 2 or not 1 <= weights.shape[0] <= 16
            or required_marks.shape != (2,) or optional_marks.ndim != 1
            or np.any(~np.isfinite(required_marks)) or np.any(~np.isfinite(optional_marks))
            or np.any(np.isnan(weights)) or np.any(np.isposinf(weights))):
        raise ValueError('Two required classes and fixed finite-or-zero cardinality weights required')
    states, sources, targets, coefficients = coverage_geometry(*required)
    first, second = required
    if first+second+len(optional_marks) >= weights.shape[1] or weights.shape[1] > 17:
        raise ValueError('Cardinality weights must retain the complete original coordinate domain')
    elementary = state.elementary_logs(optional_marks)
    selected = np.empty((weights.shape[0], (first+1)*(second+1)))
    for i in range(first+1):
        for j in range(second+1):
            selected[:, i*(second+1)+j] = (i*required_marks[0]+j*required_marks[1]
                +state.logsum(weights[:, i+j:i+j+len(elementary)]+elementary, axis=1))
    # Binomial coefficients are positive and <=2^16. Scale each type;
    # fall back to logarithmic sums if binary64 would erase a positive
    # input. Thus exact structural zero is never a numerical floor.
    scales = np.max(selected, axis=1)
    active = np.isfinite(scales)
    linear = np.zeros_like(selected)
    linear[active] = np.exp2(selected[active]-scales[active, None])
    erased = np.any(np.isfinite(selected) & (linear == 0))
    values = np.full((weights.shape[0], len(sources)), -np.inf)
    if erased:
        with np.errstate(divide='ignore'):
            log_coefficients = np.log2(coefficients)
        for a in range(weights.shape[0]):
            values[a] = state.logsum(log_coefficients+selected[a], axis=1)
    else:
        positive = linear@coefficients.T
        with np.errstate(divide='ignore'):
            values[active] = np.log2(positive[active])+scales[active, None]
    result = np.full((weights.shape[0], len(states), len(states)), -np.inf)
    result[:, sources, targets] = values
    return result


def cover_norm_from_transitions(transitions, groups, *, columns=16, check_time=None):
    """Minkowski over coverage histories, then one norm per fixed type."""
    transitions = np.asarray(transitions, dtype=float)
    if (transitions.ndim != 3 or transitions.shape[1] != transitions.shape[2]
            or not 1 <= transitions.shape[1] <= 81 or not 1 <= transitions.shape[0] <= 16
            or np.any(np.isnan(transitions)) or np.any(np.isposinf(transitions))
            or type(groups) is not int or not 1 <= groups <= 32
            or type(columns) is not int or not 1 <= columns <= 16):
        raise ValueError('Bounded positive cover matrices and integer norm exponent required')
    kernel = state.logsum(columns*transitions, axis=0)/columns
    current = np.full(len(kernel), -np.inf)
    current[0] = 0.
    for _ in range(groups):
        if check_time is not None:
            check_time()
        current = state.logsum(current[:, None]+kernel, axis=0)
    return float(current[-1])


def prepared_envelope(curves, slopes, weights):
    """Require the actual off-diagonal zeros; do not infer them from a cap."""
    curves = np.asarray(curves, dtype=float)
    prepared = emptycut.prepared_envelope(curves, slopes, weights)
    if not (np.isneginf(curves[1, 0]) and np.isneginf(curves[2, 0])):
        raise ValueError('Both R01(0) and R10(0) must be exact zeros')
    prepared['coverage_type_weights'] = np.stack([np.asarray(weights[a], dtype=float).copy()
        for a in range(1, len(weights)+1)])
    prepared['off_diagonal_empty_support_exact'] = True
    return prepared


def covered_path_norms(prepared, groups, *, columns=16, check_time=None):
    """Min with the old norm PATHWISE, before any path sum or matrix power."""
    if (prepared.get('off_diagonal_empty_support_exact') is not True
            or prepared.get('regions') != 16 or type(groups) is not int or not 1 <= groups <= 32
            or type(columns) is not int or not 1 <= columns <= 16):
        raise ValueError('Explicit complete G16 empty-zero scope and bounded group count required')
    original = groups/columns*state.logsum(columns*prepared['typed'], axis=1)
    repaired = original.copy()
    cache = {}
    matrix_states = []
    for index, counts in enumerate(prepared['allowed_counts']):
        if check_time is not None:
            check_time()
        key = tuple(map(int, counts))
        if key not in cache:
            allowed00, n01, n10, n11 = key
            if n01 == n10 == 0:
                # Unique coverage history; exactly the old type norm.
                cache[key] = float(original[index])
            else:
                optional = np.repeat(prepared['slopes'][[0, 3]], [allowed00, n11])
                transitions = type_transition_logs(prepared['coverage_type_weights'], (n01, n10),
                    prepared['slopes'][[1, 2]], optional)
                matrix_states.append(transitions.shape[1])
                cache[key] = cover_norm_from_transitions(transitions, groups,
                    columns=columns, check_time=check_time)
        repaired[index] = min(original[index], cache[key])
    return repaired, dict(unique_allowed_count_classes=len(cache),
        cover_matrices_built=len(matrix_states), maximum_cover_states=max(matrix_states, default=1),
        improved_path_cut_classes=int(np.sum(repaired < original-1e-10)),
        minimum_path_norm_change=float(np.min(repaired-original)))


def evaluate_prepared(prepared, groups, mark, alpha, *, check_time=None):
    """Unrestricted fixed-type union for ONE active-group weight, FLOAT only."""
    mark, alpha = F(mark), emptycut.fractional.checked_alpha(alpha)
    if not 0 < mark < 1 or np.shape(prepared['typed'])[1:] != (16,):
        raise ValueError('A valid mark and all16 fixed types are required')
    norms, diagnostic = covered_path_norms(prepared, groups, check_time=check_time)
    kernel = exact.kernel_from_restricted_norms(prepared['classes'], prepared['path_a'], norms)
    power = state.matrix_power_log(kernel, 16)
    prefix = float(alpha)*(-emptycut.fractional.LIMIT*log2(float(mark))-log2(255))
    return dict(construction=emptycut.fractional.CONSTRUCTION, active_groups=groups,
        mark=str(mark), alpha=str(alpha), source_counts='C(8192,r)',
        all_cohorts=True, all_constituent_row_sets=True, all_source_supports=True,
        failure_log2_screen=float(prefix+log2(comb(8192, groups))+state.logsum(power[0])),
        off_diagonal_zero_support_enforced=True, positive_coverage_dp=True,
        inclusion_exclusion_used=False, pathwise_min_with_unrestricted_norm=True,
        types_fixed_across_columns=True, type_norm_exponent=16,
        markers_fixed_across_paths_and_columns=True, independently_resampled_types=False,
        empty_00_separated=True, initial_state=0, terminal_vector=[1, 1],
        regional_banks_rebuilt=0, directed=False, requires_fresh_directed_replay=True,
        continuous_weight_band_covered=False, whole_code_certificate=False,
        coverage_diagnostics=diagnostic)

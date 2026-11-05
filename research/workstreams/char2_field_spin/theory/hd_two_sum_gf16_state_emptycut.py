"""FLOAT state-path bound with an exact empty-occupancy branch for R00.

The regional curves already include their fractional exponent. For a
fixed positive slope B, a nonempty intercept C<=1 gives
    R00(m) <= B**m * (C + (1-C)*1[m==0]).
A delta branch forces EVERY group mask to avoid that region. Expanding
only the00 steps preserves each chronological state path. Symmetric
cardinality type laws then depend only on how many00 positions remain
allowed, not which positions they are. Binomial cut multiplicities are
outside the fixed-type Holder norm and its subsequent Minkowski bound.

This helper accepts fixed per-type cardinality weights, including the
new joint pre-H-weight laws. It builds no regional bank, fits no witness,
and changes no pinned provider. Further Minkowski slack means that a
smaller pointwise R00 cap need not improve the final bound automatically.
"""
from fractions import Fraction as F
from math import comb, log2

import numpy as np

import hd_two_sum_gf16_exact_arity as exact
import hd_two_sum_gf16_fractional_state as fractional
import hd_two_sum_gf16_state_envelope as state


def empty_envelope(curves, slopes):
    """Check all occupancies; return C's log and the other affine caps."""
    curves, slopes = np.asarray(curves, dtype=float), np.asarray(slopes, dtype=float)
    intercepts, slack = state.affine_envelope(curves, slopes)
    if curves.shape[1] < 2 or curves[0, 0] != 0.:
        raise ValueError('An exact R00(0)=1 and at least one nonempty occupancy are required')
    log_c = float(np.max(curves[0, 1:]-slopes[0]*np.arange(1, curves.shape[1])))
    if log_c > 0.:
        raise ValueError('This slope needs a nonempty R00 intercept greater than1')
    # log(1-C) uses expm1 near C=1; exact endpoint cases remain exact.
    log_one_minus_c = (0. if np.isneginf(log_c) else
        -np.inf if log_c == 0. else float(np.log(-np.expm1(log_c*np.log(2.)))/np.log(2.)))
    envelope = log_c+slopes[0]*np.arange(curves.shape[1])
    envelope[0] = 0.
    slack[0] = np.inf
    finite = np.isfinite(curves[0])
    slack[0, finite] = envelope[finite]-curves[0, finite]
    slack[0, np.isneginf(envelope) & np.isneginf(curves[0])] = 0.
    if np.min(slack) < -1e-10:
        raise ArithmeticError('Empty-separated envelope failed a complete occupancy check')
    return dict(log_c=log_c, log_one_minus_c=log_one_minus_c,
                intercepts=intercepts, slack=slack, slopes=slopes)


def fractional_cardinality_weights(alpha, dense_weights, log_floors=0., *,
                                   group=16, columns=16, dimension=8, minimum=9):
    """Exactly the existing fractional type laws, as per-mask weight logs.

    Each returned array is indexed by the mask cardinality, not its count.
    No C(group,w) belongs here; the elementary polynomial enumerates masks.
    Sparse type2 keeps its two singleton input categories separate before
    their fractional masses are added, as in the pinned provider.
    """
    alpha = float(fractional.checked_alpha(alpha))
    floors = fractional.floor_array(log_floors, group)
    if (type(group) is not int or not 2 <= group <= 16
            or set(dense_weights) != set(range(3, group+1))
            or type(columns) is not int or not 1 <= minimum <= columns
            or not 1 <= dimension <= columns):
        raise ValueError('Complete fixed cardinality laws and valid constituent parameters required')
    logq = log2(255)
    result = {}
    single = np.full(group+1, -np.inf)
    single[0] = (log2(group)-dimension*alpha*logq-minimum*floors[0])/columns
    single[group] = single[0]+alpha*logq+floors[0]
    result[1] = single
    pair = np.full(group+1, -np.inf)
    pair[0] = (log2(comb(group, 2))-2*dimension*alpha*logq-2*minimum*floors[1])/columns
    pair[group-1] = pair[0]+2*floors[1]+alpha*logq
    pair[group] = pair[0]+np.logaddexp2(1+alpha*logq+floors[1],
                                      2*floors[1]+alpha*(logq+log2(255-group)))
    result[2] = pair
    for arity in range(3, group+1):
        value = np.asarray(dense_weights[arity], dtype=float)
        if value.shape != (group+1,) or np.any(np.isnan(value)) or np.any(np.isposinf(value)):
            raise ValueError('Finite-or-zero dense weights at every mask cardinality required')
        value = value.copy()-minimum*floors[arity-1]/columns
        value[1:] += floors[arity-1]
        result[arity] = value
    return result


def prepared_envelope(curves, slopes, weights, *, regions=16):
    """Prepare all chronological path/cut classes with fixed type weights.

    A caller may pass fractional_cardinality_weights, or the all-arity
    arrays from input_weight.type_weights (possibly one fixed v per type).
    Those choices must remain the same for every path and column. This
    grouping requires the same regional envelope at every position and
    type laws symmetric in the region labels; labelled laws are excluded.
    """
    if (type(regions) is not int or not 1 <= regions <= 16
            or not isinstance(weights, dict) or not 1 <= len(weights) <= 16
            or set(weights) != set(range(1, len(weights)+1))):
        raise ValueError('At most16 regions and consecutive fixed type weights required')
    laws = np.stack([np.asarray(weights[a], dtype=float) for a in range(1, len(weights)+1)])
    if (laws.shape != (len(weights), regions+1) or np.any(np.isnan(laws))
            or np.any(np.isposinf(laws))):
        raise ValueError('Finite-or-zero per-mask cardinality laws required')
    envelope = empty_envelope(curves, slopes)
    log_c, log_d = envelope['log_c'], envelope['log_one_minus_c']
    classes, path_a, allowed_counts = [], [], []
    for initial, final, counts, multiplicity in state.path_classes(regions):
        for allowed in range(counts[0]+1):
            forbidden = counts[0]-allowed
            if (allowed and np.isneginf(log_c)) or (forbidden and np.isneginf(log_d)):
                continue
            intercept = float(np.dot(counts[1:], envelope['intercepts'][1:]))
            if allowed:
                intercept += allowed*log_c
            if forbidden:
                intercept += forbidden*log_d
            classes.append((initial, final, counts, multiplicity*comb(counts[0], allowed)))
            path_a.append(intercept)
            allowed_counts.append((allowed, *counts[1:]))
    allowed_counts = np.asarray(allowed_counts, dtype=np.int32)
    lengths = allowed_counts.sum(axis=1)
    elementary = np.full((len(classes), regions+1), -np.inf)
    # Bucket equal lengths: no fake positive marker replaces a forbidden
    # coordinate, and no change to the pinned elementary_logs is needed.
    for length in range(regions+1):
        indices = np.flatnonzero(lengths == length)
        if not len(indices):
            continue
        markers = np.asarray([np.repeat(envelope['slopes'], allowed_counts[i]) for i in indices])
        elementary[indices, :length+1] = state.elementary_logs(markers)
    typed = np.stack([state.logsum(elementary+law, axis=1) for law in laws], axis=1)
    return dict(classes=tuple(classes), path_a=np.asarray(path_a), typed=typed,
        allowed_counts=allowed_counts, elementary=elementary, regions=regions,
        log_c=log_c, log_one_minus_c=log_d, slopes=envelope['slopes'],
        intercepts=envelope['intercepts'], slack=envelope['slack'],
        original_state_path_classes=len(state.path_classes(regions)),
        state_cut_classes=len(classes), zero_occupancy_00_exact=True,
        fixed_symmetric_type_laws=True, regional_banks_rebuilt=0)


def evaluate_prepared(prepared, groups, mark, alpha, *, cutoff=None, check_time=None):
    """One G16 active-group class, optionally restricted to A_tot<cutoff."""
    mark, alpha = F(mark), fractional.checked_alpha(alpha)
    if (prepared.get('regions') != 16 or np.shape(prepared['typed'])[1:] != (16,)
            or type(groups) is not int or not 1 <= groups <= 32 or not 0 < mark < 1):
        raise ValueError('Complete G16 types, at most32 groups, and a valid weight mark required')
    if cutoff is not None:
        row = exact.evaluate_prepared(prepared, groups, cutoff, mark, alpha, check_time=check_time)
    else:
        if check_time is not None:
            check_time()
        kernel = state.kernel_from_types(prepared['classes'], prepared['path_a'], prepared['typed'], groups)
        power = state.matrix_power_log(kernel, 16)
        prefix = float(alpha)*(-fractional.LIMIT*log2(float(mark))-log2(255))
        row = dict(construction=fractional.CONSTRUCTION, active_groups=groups, mark=str(mark), alpha=str(alpha),
            source_counts='C(8192,r)', all_cohorts=True, all_constituent_row_sets=True,
            failure_log2_screen=float(prefix+log2(comb(8192, groups))+state.logsum(power[0])),
            exact_total_arity_coefficient_restriction=False, initial_state=0, terminal_vector=[1, 1],
            directed=False, whole_code_certificate=False, continuous_weight_band_covered=False,
            requires_fresh_directed_replay=True, regional_banks_rebuilt=0)
    row.update(empty_occupancy_00_separated=True, zero_occupancy_00_exact=True,
        cut_coefficient_powered_again=False, allowed_cut_classes=prepared['state_cut_classes'],
        log2_nonempty_00_intercept=prepared['log_c'],
        log2_empty_00_coefficient=prepared['log_one_minus_c'],
        chronology='regional cap, positive empty-cut expansion, fixed-type Holder/Minkowski',
        type_weights_fixed_across_paths_and_columns=True, labelled_type_laws_supported=False)
    return row

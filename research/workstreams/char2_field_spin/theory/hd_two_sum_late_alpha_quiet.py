"""One bounded FLOAT column-late-alpha quiet-contribution control.

Freeze r18, arity8, z691/700 and alpha1/5. Build one FRESH unpowered
forced-NZ regional bank. Its 00 entries are bounded by B**m with intercept
one, uniformly over source-pair histograms. In the quiet outer branch,
input/output supports agree, so only total activity enters this envelope.

For one fixed row set, Z2(w) is the number of nonzero assignments on each
weight-w mask whose two moments vanish. The late comparison first sums all
labelled masks AND all r groups at each total column weight; it then takes
alpha. The early control powers each mask count before these sums. Both
use exactly the same R1 envelope, total9ar floor, and normalizations.

This is one term of a proposed failure-probability majorant. Other outer
syndromes and state paths are omitted. It is NOT a lower bound on actual
failure probability, not a complete quiet-event bound, and not a distance
certificate. No old fractional regional table is rooted to obtain R1.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_closed_kernel import zero_two_moments
from hd_two_sum_gf16_factorized import forced_nonzero_regions
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS
from hd_xor_multigroup import LIMIT


GROUPS, ARITY, COLUMNS, DIMENSION, MINIMUM = 18, 8, 16, 8, 9
MARK, ALPHA = F(691, 700), F(1, 5)


def _logs(values):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 1 or not len(values) or np.any(np.isnan(values))
            or np.any(np.isposinf(values))):
        raise ValueError('A nonempty finite-or-zero log coefficient array is required')
    return values


def log_convolve(left, right, check_time=lambda: None):
    """Positive ordinary convolution, never FFT or signed cancellation."""
    left, right = _logs(left), _logs(right)
    if len(left) > len(right):
        left, right = right, left
    result = np.full(len(left)+len(right)-1, -np.inf)
    for offset in np.flatnonzero(np.isfinite(left)):
        check_time()
        target = result[offset:offset+len(right)]
        np.logaddexp2(target, left[offset]+right, out=target)
    return result


def log_power(values, exponent, check_time=lambda: None):
    if type(exponent) is not int or exponent < 0:
        raise ValueError('Nonnegative integer exponent required')
    result, base = np.array([0.]), _logs(values)
    while exponent:
        check_time()
        if exponent & 1:
            result = log_convolve(result, base, check_time)
        exponent >>= 1
        if exponent:
            base = log_convolve(base, base, check_time)
    return result


def quiet_mask_logs(arity, log_marker, alpha=1., *, q=256):
    """Mask multiplicity C(a,w) stays OUTSIDE the fractional count power."""
    if (type(arity) is not int or not 2 <= arity <= min(q, 16)
            or not np.isfinite(log_marker) or not 0 < alpha <= 1
            or type(q) is not int or q < 4 or q & (q-1)):
        raise ValueError('Valid arity, field order, marker and exponent required')
    counts = [zero_two_moments(w, q) for w in range(arity+1)]
    if min(counts) < 0 or sum(comb(arity, w)*count for w, count in enumerate(counts)) != q**(arity-2):
        raise ArithmeticError('Two-moment kernel count normalization failed')
    return np.array([log2(comb(arity, w))+alpha*(log2(count)+w*log_marker)
                     if count else -np.inf for w, count in enumerate(counts)])


def column_coefficients(groups, arity, log_marker, alpha, *, q=256,
                        columns=16, dimension=8, check_time=lambda: None):
    """Two normalized column polynomials with the SAME aggregate degrees.

    A whole word has prescribed-support denominator(q-1)**(dimension*a*r).
    Distribute its alpha power uniformly across columns. This is outside
    all mask sums, and it is charged only once by the final column power.
    """
    if (type(groups) is not int or groups < 1 or type(columns) is not int
            or type(dimension) is not int or not 1 <= dimension <= columns
            or not 0 < alpha <= 1):
        raise ValueError('Positive groups/columns and valid dimension/exponent required')
    unpowered = log_power(quiet_mask_logs(arity, log_marker, q=q), groups, check_time)
    late = alpha*unpowered
    early = log_power(quiet_mask_logs(arity, log_marker, alpha, q=q), groups, check_time)
    normalizer = alpha*dimension*arity*groups/columns*log2(q-1)
    return late-normalizer, early-normalizer


def compare_contributions(groups, arity, log_marker, alpha, *, q=256,
                          columns=16, dimension=8, minimum=9,
                          check_time=lambda: None):
    if type(minimum) is not int or not 1 <= minimum <= columns:
        raise ValueError('An aggregate per-row support floor in1..columns is required')
    late, early = column_coefficients(groups, arity, log_marker, alpha, q=q,
        columns=columns, dimension=dimension, check_time=check_time)
    if np.any(late > early+1e-9):
        raise ArithmeticError('Positive subadditivity control failed')
    floor = minimum*arity*groups
    late_word = log_power(late, columns, check_time)
    early_word = log_power(early, columns, check_time)
    return dict(late_column_log2=float(np.logaddexp2.reduce(late_word[floor:])),
                early_mask_log2=float(np.logaddexp2.reduce(early_word[floor:])),
                aggregate_input_weight_floor=floor,
                column_degree=groups*arity, word_degree=groups*arity*columns,
                individual_row_floors=False, all_rows_seen_constraint=False)


def zero_intercept_envelope(regions, maximum):
    """R1_00(0)=1 exactly; B=max_{m>0,f} R1_00(m,f)**(1/m).

    This is a finite envelope construction, not a continuous fit. Exact
    zero entries remain zero in the input. Complete occupancy coverage is
    checked; numerical values are still FLOAT and are not directed bounds.
    """
    required = {(m, f) for m in range(maximum+1) for f in range(m//2+1)}
    if type(maximum) is not int or maximum < 1 or set(regions) != required:
        raise ValueError('Complete regional matching-histogram table required')
    curves = []
    for m in range(maximum+1):
        matrices = np.asarray([regions[m, f] for f in range(m//2+1)])
        if (matrices.shape != (m//2+1, 2, 2) or np.any(np.isnan(matrices))
                or np.any(np.isposinf(matrices))):
            raise ValueError('Finite-or-zero log regional matrices required')
        curves.append(float(np.max(matrices[:, 0, 0])))
    if not np.isfinite(curves[0]) or abs(curves[0]) > 1e-10:
        raise ArithmeticError('Empty-input zero-state transfer must equal one')
    ratios = np.array(curves[1:])/np.arange(1, maximum+1)
    slope = float(np.max(ratios))
    if not np.isfinite(slope):
        raise ArithmeticError('A positive nonempty00 cap is required')
    slack = slope*np.arange(maximum+1)-curves
    if np.min(slack) < -1e-9:
        raise ArithmeticError('Complete occupancy-envelope check failed')
    return dict(log2_intercept=0., log2_slope=slope,
                slope_as_exact_float=str(F(slope)),
                log2_curves=[value if np.isfinite(value) else None for value in curves],
                null_log_curve_means_exact_zero=True,
                contacts=np.flatnonzero(slack < 1e-9).tolist(),
                minimum_float_slack=float(np.min(slack)))


def run(*, max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit positive budget at most30 seconds required')
    start = perf_counter()
    def check_time():
        if perf_counter()-start >= max_seconds:
            raise TimeoutError('Late-alpha quiet control exhausted its budget')
    result = dict(construction=CONSTRUCTION, labels=list(LABELS),
        field_order=256, active_groups=GROUPS, arity=ARITY, mark=str(MARK), alpha=str(ALPHA),
        columns=COLUMNS, constituent_dimension=DIMENSION, forbidden_weight=LIMIT,
        scope='all column moments zero; all regional boundary states zero',
        full_quiet_event_bounded=False, whole_code_certificate=False,
        failure_probability_lower_bound=False, directed=False,
        all_source_supports_in_chosen_term=True,
        all_fixed_arity_row_sets_in_chosen_term=True,
        source_count=str(comb(8192, GROUPS)), row_set_count_per_group=comb(16, ARITY),
        fixed_row_sets_across_columns=True, row_sets_resampled=False,
        fresh_regional_alpha=1, packet_refresh_keys_already_averaged=True,
        regional_envelope_entry='00', regional_envelope_intercept=1,
        regional_geometry_envelope='max over every feasible full-pair count',
        mask_sum_before_alpha=True, group_mask_sum_before_alpha=True,
        column_degree_allocations_expanded_before_alpha=True,
        unconditional_jensen=True, route_fractionalization_retained=False,
        normalizer_before_alpha=f'255^(-{DIMENSION*ARITY*GROUPS}) per whole word',
        projective_normalizations=1, source_and_rowset_counts_outside_alpha=True,
        old_control_uses_same_fresh_unpowered_envelope=True,
        raw_regional_banks=0, mark_fits=0, slope_fits=0,
        compute_budget_seconds=max_seconds, complete_requested_control=False)
    try:
        check_time()
        regions = forced_nonzero_regions(float(MARK), GROUPS)
        result['raw_regional_banks'] = 1
        check_time()
        envelope = zero_intercept_envelope(regions, GROUPS)
        result['regional_envelope'] = envelope
        if progress:
            progress(dict(stage='fresh-R1-built', elapsed_seconds=perf_counter()-start,
                          contacts=envelope['contacts'], log2_slope=envelope['log2_slope']))
        comparison = compare_contributions(GROUPS, ARITY, envelope['log2_slope'],
            float(ALPHA), check_time=check_time)
        prefix = (log2(comb(8192, GROUPS))+GROUPS*log2(comb(16, ARITY))
                  +float(ALPHA)*(-LIMIT*log2(float(MARK))-log2(255)))
        comparison.update(log2_global_prefix=prefix,
            late_majorant_term_log2=prefix+comparison['late_column_log2'],
            early_majorant_term_log2=prefix+comparison['early_mask_log2'],
            late_minus_early_bits=comparison['late_column_log2']-comparison['early_mask_log2'])
        result['comparison'] = comparison
        result['complete_requested_control'] = True
    except TimeoutError as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-start
    directory = Path(__file__).resolve().parent
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    result['direct_dependency_sha256'] = {name: sha256((directory/name).read_bytes()).hexdigest()
        for name in ('hd_two_sum_closed_kernel.py', 'hd_two_sum_gf16_factorized.py')}
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--max-seconds', type=float, default=30)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen control')
    answer = run(max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(answer.get('comparison', answer.get('incomplete_reason')), allow_nan=False))

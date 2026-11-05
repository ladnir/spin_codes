"""FLOAT all-cohort G16 gate via an affine regional occupancy envelope.

This new scope includes every constituent-row set and every source matching
class in a requested active-group band. It is a proposal, not a certificate.
Only the scalar affine envelope is selected; z and the state potential are
fixed. Dense arities are summed AFTER each group's16-column tail.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize_scalar

from hd_two_sum_closed_kernel import arity_holder_moments
from hd_two_sum_gf16_full_dense import Z, LOGP
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS
from hd_xor_dense_holder_classes import class_log_weights
from hd_xor_multigroup import LIMIT, nonzero_matrices, profile_arrays
from hd_shifted_swap_sparse import sparse_geometry, power_screen, valid_degrees


def logsum(values):
    values = tuple(values)
    return float(np.logaddexp2.reduce(values)) if values else -np.inf


def times(count, log_value):
    return 0. if count == 0 else count*log_value


def forced_nonzero_regions(z, maximum):
    """All R[m,f] through maximum64; actual NZ-amplitude alpha1 caps only.

    Each packet still has four RS[16,4] blocks and32 source pairs. The
    positive down-set s+2f<=maximum is valid independently of the old17 guard.
    There is no import of uniform-amplitude/fractional-power packet caps.
    """
    if type(maximum) is not int or not 1 <= maximum <= 64 or not 0 < z < 1:
        raise ValueError('Expected maximum1..64 and mark in(0,1)')
    _, _, matrices = nonzero_matrices(float(z))
    occupancies, _ = profile_arrays()
    table = {tuple(occ): matrix.ravel() for occ, matrix in zip(occupancies, matrices)}
    base = np.zeros((4, maximum+1, maximum//2+1))
    marker = maximum/4096
    for single, full, count, occupancy in sparse_geometry(maximum):
        base[:, single, full] += count/2**single*table[occupancy]*marker**(single+full)
    powered = power_screen(base, 128, maximum)
    if np.any(~np.isfinite(powered)) or np.any(powered < 0):
        raise ArithmeticError('Invalid FLOAT positive regional coefficients')
    result = {}
    with np.errstate(divide='ignore'):
        for single, full in valid_degrees(maximum):
            normalizer = log2(comb(4096, single)*comb(4096-single, full))+(single+full)*log2(marker)
            result[single+2*full, full] = (np.log2(powered[:, single, full])-normalizer).reshape(2, 2)
    return result


def contraction_envelope(regions, maximum, log_potential):
    """log c_m=max_f,state log((R[m,f]v)_state/v_state)."""
    required = {(m, f) for m in range(maximum+1) for f in range(m//2+1)}
    if set(regions) != required or not np.isfinite(log_potential) or log_potential > 0:
        raise ValueError('Complete ordered pair histograms and0<p<=1 required')
    potential = np.array([0., log_potential])
    result = []
    for m in range(maximum+1):
        caps = []
        for f in range(m//2+1):
            matrix = np.asarray(regions[m, f])
            if matrix.shape != (2, 2) or np.any(np.isnan(matrix)) or np.any(np.isposinf(matrix)):
                raise ValueError('Finite-or-zero log matrices required')
            action = np.logaddexp2.reduce(matrix+potential[None, :], axis=1)-potential
            caps.append(float(np.max(action)))
        result.append(max(caps))
    if any(not np.isfinite(value) for value in result):
        raise ArithmeticError('Positive finite regional contractions required')
    return np.array(result)


@lru_cache(None)
def sparse_pair_terms(columns=16, minimum=9):
    """Exact source-column multiplicities for EACH row's minimum support.

    Key(s,d) counts column words over00,10,01,11 with s singleton and d
    double columns. Both labelled source rows have support>=minimum.
    """
    if type(columns) is not int or type(minimum) is not int or not 1 <= minimum <= columns:
        raise ValueError('Positive individual row floor required')
    result = {}
    for left in range(minimum, columns+1):
        for right in range(minimum, columns+1):
            for double in range(max(0, left+right-columns), min(left, right)+1):
                single = left+right-2*double
                multiplicity = comb(columns, double)*comb(columns-double, left-double)*comb(columns-left, right-double)
                result[single, double] = result.get((single, double), 0)+multiplicity
    return tuple((s, d, count) for (s, d), count in sorted(result.items()))


def positive_binomial_tail(empty_log, nonempty_log, columns=16, minimum=9):
    return logsum(log2(comb(columns, j))+times(columns-j, empty_log)+times(j, nonempty_log)
                  for j in range(minimum, columns+1))


def dense_arity_logs():
    result, pins = {}, {}
    for arity in range(3, 17):
        numerators, denominator = arity_holder_moments(arity, 16, 16, 8)
        result[arity] = class_log_weights(numerators, denominator, 16)
        pins[arity] = sha256(json.dumps([list(numerators), denominator], separators=(',', ':')).encode()).hexdigest()
    return result, pins


def group_bounds(log_b, dense_weights):
    """Log whole-word counts for one active group after regional domination.

    Sparse row choices16/120 and dense row choices inside Holder enter once.
    Sparse individual floors9 and dense nonempty-column floor9 remain exact
    restrictions of their respective positive comparisons.
    """
    if not np.isfinite(log_b) or set(dense_weights) != set(range(3, 17)):
        raise ValueError('Finite occupancy slope and all14 dense arities required')
    logq = log2(255)
    single_action = 16*log_b
    pair_action = logsum((log2(239)+16*log_b, 4+15*log_b))-logq
    one = log2(16)-8*logq+logsum(log2(comb(16, w))+w*(logq+single_action) for w in range(9, 17))
    two = log2(120)-16*logq+logsum(log2(count)+(s+2*d)*logq+s*single_action+d*pair_action
                for s, d, count in sparse_pair_terms())
    dense = {}
    for arity, weights in dense_weights.items():
        if len(weights) != 17 or np.any(np.isnan(weights)) or np.any(np.isposinf(weights)):
            raise ValueError('Complete finite-or-zero cardinality weights required')
        nonempty = logsum(log2(comb(16, w))+weights[w]+w*log_b for w in range(1, 17))
        dense[arity] = positive_binomial_tail(weights[0], nonempty)
    return dict(singleton=one, pair=two, dense_by_arity=dense,
                dense_union=logsum(dense.values()), total=logsum((one, two, *dense.values())))


def all_cohort_band(log_a, group_log, first, last, z, log_potential):
    """All source supports/cohorts via binomial source count and group power."""
    if not 1 <= first <= last <= 64 or not 0 < z < 1 or not log_potential <= 0:
        raise ValueError('Bounded active-group band and valid fixed witness required')
    prefix = 256*log_a-LIMIT*log2(float(z))-log2(255)-log_potential
    rows = [dict(active_groups=r, failure_log2_screen=prefix+log2(comb(8192, r))+r*group_log)
            for r in range(first, last+1)]
    return logsum(row['failure_log2_screen'] for row in rows), rows


def run(*, first=18, last=32, max_seconds=20, progress=None):
    """One inherited z/p; bounded1D affine-envelope selection, no z/p fit."""
    if (type(first) is not int or type(last) is not int or not 18 <= first <= last <= 64
            or not 0 < max_seconds <= 20):
        raise ValueError('New band18..64 and compute budget at most20 seconds required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Bounded factorized FLOAT budget exhausted')
    result = dict(construction=CONSTRUCTION, field_order=256, labels=list(LABELS),
        forbidden_weight=LIMIT, active_group_band=[first, last], total_columns=16,
        constituent_dimension=8, group_size=16, arbitrary_dense_group_count=True,
        all_cohort_splits=True, all_constituent_row_sets=True, all_source_matching_classes=True,
        source_counts='C(8192,r)', cohort_counts='generated by (G1+G2+sum_a GH_a)^r',
        dense_arities_summed_after_group_word_power=True, dense_arities_merged_before_holder=False,
        sparse_individual_row_floor=9, dense_individual_nonempty_column_floor=9,
        fixed_row_sets_across_columns=True, independent_rowset_resampling_assumed=False,
        regional_cap_amplitude_law='independent-uniform-nonzero', regional_alpha=1,
        regional_pair_histograms_complete=True, regional_potential_before_column_mixture=True,
        occupancy_envelope='R[m,f]v <= A*B^m*v for every m<=maximum and f',
        projective_normalizations=1, endpoint_potential_factors=1, pair_swaps=True,
        mark=str(Z), log2_potential=str(LOGP), z_p_fits=0, directed=False,
        old_certificates_transferred=False, whole_code_certificate=False,
        complete_requested_band=False, compute_budget_seconds=max_seconds)
    trials = []
    try:
        regions = forced_nonzero_regions(float(Z), last)
        caps = contraction_envelope(regions, last, float(LOGP))
        weights, pins = dense_arity_logs()
        check_time()
        bank_seconds = perf_counter()-started
        def objective(log_b):
            check_time()
            log_a = float(np.max(caps-log_b*np.arange(last+1)))
            groups = group_bounds(log_b, weights)
            bound, _ = all_cohort_band(log_a, groups['total'], first, last, Z, float(LOGP))
            trials.append(dict(log2_b=float(log_b), log2_a=log_a, failure_log2_screen=bound))
            return bound
        slope_bounds = (-64., 8.)
        optimized = minimize_scalar(objective, bounds=slope_bounds, method='bounded',
                                    options=dict(maxiter=64, xatol=1e-9))
        objective(slope_bounds[0])
        objective(slope_bounds[1])
        best = min(trials, key=lambda row: row['failure_log2_screen'])
        log_a, log_b = best['log2_a'], best['log2_b']
        groups = group_bounds(log_b, weights)
        bound, rows = all_cohort_band(log_a, groups['total'], first, last, Z, float(LOGP))
        slack = log_a+log_b*np.arange(last+1)-caps
        if np.min(slack) < -1e-10:
            raise ArithmeticError('Selected FLOAT envelope failed its complete occupancy check')
        check_time()
        result.update(status='completed', complete_requested_band=True, failure_log2_screen=bound,
            entries=rows, exact_dense_moment_sha256=pins, group_log_bounds=groups,
            regional_log_contractions=caps.tolist(), envelope_minimum_slack=float(np.min(slack)),
            log2_a_exact_float=str(F(log_a)), log2_b_exact_float=str(F(log_b)),
            affine_slope_bounds=list(slope_bounds), affine_trials=trials,
            affine_search_success=bool(optimized.success), affine_evaluations=len(trials),
            winning_slope_hits_boundary=min(abs(log_b-edge) for edge in slope_bounds) < 1e-7,
            raw_regional_banks=1, bank_seconds=bank_seconds)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error), affine_trials=trials)
    result['elapsed_seconds'] = perf_counter()-started
    if progress:
        progress({key: value for key, value in result.items()
                  if key not in ('affine_trials', 'entries', 'exact_dense_moment_sha256', 'regional_log_contractions')})
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=int, default=18)
    parser.add_argument('--last', type=int, default=32)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite FLOAT evidence')
    result = run(first=args.first, last=args.last)
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                     if key not in ('affine_trials', 'entries', 'exact_dense_moment_sha256', 'regional_log_contractions')}, allow_nan=False))

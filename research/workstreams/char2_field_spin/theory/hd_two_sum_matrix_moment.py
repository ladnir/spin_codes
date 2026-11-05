"""Fixed four-atom regional measure; positive whole-column first moment.

No provider or optimizer. The00 entry keeps the frozen atom0/B0 law;
01/10 share ONE2/18 secant atom;11 uses one constant cap. Exact binary
state-path multiplicities expand only the00 mixture. Fixed-A Holder
weights include all row-set choices once; chronological column matrix
powers and the aggregate surrogate input floor are retained.

Scalar mode uses one whole-word Chernoff marker derived from the saved
quiet comparator's adjacent floor coefficients. Polynomial mode retains
the exact positive degree tail. Both are FLOAT, not certificates.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, expm1, log, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_gf16_mean_column import raw_holder_weights
from hd_two_sum_gf16_state_envelope import path_classes, elementary_logs
from hd_two_sum_late_alpha_live import log_matrix_power
from hd_two_sum_late_alpha_moment import quiet_mixture_coefficients
from hd_two_sum_late_alpha_quiet import log_convolve, log_power
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
ARTIFACT = 'hd_two_sum_late_alpha_live_gate_float.json'
ARTIFACT_SHA = '33da1eb9517c7d5665d22d2c8b93f5622199d4e1b70dd2a8e3aa52e28d4e04cf'
SOURCE_PINS = {
    'hd_two_sum_late_alpha_live.py': 'e28bd5473fe60be32a984e67cbd0b64f52078a30aa4f6ff4b5b304497fd4fb59',
    'hd_two_sum_gf16_mean_column.py': 'fa6377002fa52fbfcb980ddd4b72d30df16b6c12c7c6db24eed54554e61a1593',
    'hd_two_sum_gf16_state_envelope.py': '5107b945ab18415191c48a3434a354c101d8fc5aeda66ec83a3fd2a5c3498af7',
    'hd_two_sum_closed_kernel.py': '2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979',
    'hd_two_sum_late_alpha_quiet.py': '8a753036a9be44759206b5beb3baa044f28d80b265a9304da81bdde6e8f11ae3',
    'hd_two_sum_late_alpha_moment.py': 'fe9a8936e99f82a4aec31f903bc55ad6ebfec7fa8a67fdb30f6e37bc4f2ca8e6',
}
GROUPS, ARITY, REGIONS, COLUMNS = 18, 8, 16, 16
DIMENSION, MINIMUM, FLOOR = 8, 9, 1296
MARK = F(691, 700)


def entry_measure(curves, logp, logb, *, cross_endpoints=(2, 18)):
    """One prescribed four-atom positive matrix envelope, no slope search."""
    curves = np.asarray(curves, dtype=float)
    left, right = cross_endpoints
    if (curves.ndim != 2 or curves.shape[0] != 4
            or np.any(np.isnan(curves)) or np.any(np.isposinf(curves))
            or not np.isfinite(logp) or logp > 0 or not np.isfinite(logb)
            or not 1 <= left < right < curves.shape[1]
            or not np.isfinite(curves[1, left]) or not np.isfinite(curves[1, right])):
        raise ValueError('Complete four curves and fixed finite atom parameters required')
    m = np.arange(curves.shape[1])
    zero = logp+logb*m
    zero[0] = 0.
    cross = float((curves[1, right]-curves[1, left])/(right-left))
    scale01 = float(np.max(curves[1]-cross*m))
    scale10 = float(np.max(curves[2]-cross*m))
    scale11 = float(np.max(curves[3]))
    bound = np.stack((zero, scale01+cross*m, scale10+cross*m,
                      np.full(len(m), scale11)))
    if (np.any(~np.isfinite(bound)) or np.any(bound < curves-1e-9)):
        raise ArithmeticError('The prescribed matrix measure does not cover every entry')
    return dict(log2_probability=float(logp), log2_atom00=float(logb),
        log2_atom_cross=cross, log2_scale01=scale01,
        log2_scale10=scale10, log2_scale11=scale11,
        cross_secant_endpoints=list(cross_endpoints), atom00_fitted=False,
        contacts=[np.flatnonzero(row < 1e-9).tolist() for row in bound-curves],
        minimum_slack=float(np.min(bound-curves)),
        exact_float_parameters={key: str(F(value)) for key, value in
            (('p_log', logp), ('b00_log', logb), ('cross_log', cross),
             ('scale01_log', scale01), ('scale10_log', scale10), ('scale11_log', scale11))})


def group_polynomial(weights, counts, active00, measure, arity):
    """One group; zero atoms forbid masks at the remaining00 positions."""
    weights = np.asarray(weights)
    n00, n01, n10, n11 = counts
    if (weights.shape != (sum(counts)+1,) or not 0 <= active00 <= n00
            or not 1 <= arity < len(weights) or np.any(np.isnan(weights))
            or np.any(np.isposinf(weights))):
        raise ValueError('Complete mask weights and feasible path/atom counts required')
    markers = np.array([measure['log2_atom00']]*active00
        +[measure['log2_atom_cross']]*(n01+n10)+[0.]*n11)
    elementary = elementary_logs(markers)
    result = np.full(arity+1, -np.inf)
    for weight in range(len(elementary)):
        degree = min(weight, arity)
        result[degree] = np.logaddexp2(result[degree], weights[weight]+elementary[weight])
    return result


def column_matrix(weights, groups, measure, *, arity, log_marker=None,
                  check_time=lambda: None):
    """Exact ordered path expansion of the stated positive matrix measure.

    log_marker=None returns a2x2 polynomial. Otherwise return its exact
    positive scalar evaluation atX=2**log_marker, not a per-path optimum.
    The cached group polynomial depends only on nonzero atom counts.
    """
    weights = np.asarray(weights)
    size = len(weights)-1
    if (type(groups) is not int or groups < 1 or weights.ndim != 1
            or not 1 <= size <= 16 or not 1 <= arity <= size
            or (log_marker is not None and (not np.isfinite(log_marker) or log_marker < 0))):
        raise ValueError('Valid group/region counts and a fixed X>=1 required')
    degree = groups*arity
    result = np.full((2, 2, degree+1) if log_marker is None else (2, 2), -np.inf)
    logp = measure['log2_probability']
    if not np.isfinite(logp) or logp > 0:
        raise ValueError('00 atom probability must belong to(0,1]')
    logother = -np.inf if logp == 0 else log2(-expm1(logp*log(2)))
    cached = {}
    terms = 0
    for initial, final, counts, multiplicity in path_classes(size):
        check_time()
        n00, n01, n10, n11 = counts
        state_mass = (log2(multiplicity)+n01*measure['log2_scale01']
            +n10*measure['log2_scale10']+n11*measure['log2_scale11'])
        for active00 in range(n00+1):
            if logp == 0 and active00 != n00:
                continue
            atom_mass = (log2(comb(n00, active00))+active00*logp
                         +((n00-active00)*logother if active00 < n00 else 0.))
            key = (active00, n01+n10, n11)
            if key not in cached:
                polynomial = group_polynomial(weights, counts, active00, measure, arity)
                cached[key] = (log_power(polynomial, groups, check_time) if log_marker is None else
                    groups*float(np.logaddexp2.reduce(polynomial+log_marker*np.arange(arity+1))))
            result[initial, final] = np.logaddexp2(result[initial, final],
                                                cached[key]+atom_mass+state_mass)
            terms += 1
    return result, dict(path_classes=len(path_classes(size)), positive_path_atom_terms=terms,
                        distinct_group_polynomials=len(cached))


def matrix_poly_product(left, right, check_time=lambda: None):
    left, right = np.asarray(left), np.asarray(right)
    if (left.ndim != 3 or right.ndim != 3 or left.shape[:2] != (2, 2)
            or right.shape[:2] != (2, 2) or np.any(np.isnan(left))
            or np.any(np.isnan(right)) or np.any(np.isposinf(left))
            or np.any(np.isposinf(right))):
        raise ValueError('Two finite-or-zero2x2 positive log polynomials required')
    result = np.full((2, 2, left.shape[2]+right.shape[2]-1), -np.inf)
    for i in range(2):
        for k in range(2):
            for j in range(2):
                check_time()
                if np.any(np.isfinite(left[i, k])) and np.any(np.isfinite(right[k, j])):
                    term = log_convolve(left[i, k], right[k, j], check_time)
                    np.logaddexp2(result[i, j], term, out=result[i, j])
    return result


def matrix_poly_power(matrix, exponent, check_time=lambda: None):
    if type(exponent) is not int or exponent < 1:
        raise ValueError('Positive integer matrix-polynomial exponent required')
    result = np.array([[[0.], [-np.inf]], [[-np.inf], [0.]]])
    base = np.asarray(matrix)
    while exponent:
        if exponent & 1:
            result = matrix_poly_product(result, base, check_time)
        exponent >>= 1
        if exponent:
            base = matrix_poly_product(base, base, check_time)
    return result


def adjacent_floor_marker(word, floor):
    """ONE frozen rule: max(1,a_floor/a_(floor+1)); no neighbor search."""
    word = np.asarray(word)
    if (word.ndim != 1 or type(floor) is not int or not 0 <= floor < len(word)-1
            or not np.isfinite(word[floor]) or not np.isfinite(word[floor+1])):
        raise ValueError('Two finite adjacent whole-word coefficients required')
    return max(0., float(word[floor]-word[floor+1]))


def load_saved():
    raw = (DIRECTORY/ARTIFACT).read_bytes()
    if sha256(raw).hexdigest() != ARTIFACT_SHA:
        raise ValueError('Saved fresh four-curve bank changed')
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed dependency changed: '+name)
    data = json.loads(raw)
    required = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=1,
        mark=str(MARK), alpha='1', regional_alpha=1, complete_requested_control=True,
        matrix_entry_order=['00', '01', '10', '11'], directed=False, whole_code_certificate=False)
    if any(data.get(key) != value for key, value in required.items()):
        raise ValueError('Exact homogeneous a8/F1 unpowered bank required')
    return data, np.array([[-np.inf if x is None else x for x in row]
                          for row in data['regional_log2_curves']])


def run(*, mode='scalar', quiet=False, max_seconds=30, progress=None):
    if (mode not in ('scalar', 'polynomial', 'both') or type(quiet) is not bool
            or type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30):
        raise ValueError('A fixed scope/mode and explicit budget at most30 seconds required')
    data, curves = load_saved()
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Four-atom saved-bank control exhausted its budget')
    old = data['two_atom_envelope']
    measure = entry_measure(curves, old['log2_probability'], old['log2_nonzero_atom'])
    weights = raw_holder_weights(quiet=quiet)
    prefix = (log2(int(data['source_count']))-LIMIT*log2(float(MARK))-log2(255)
              -DIMENSION*ARITY*GROUPS*log2(255))
    result = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=1,
        mark=str(MARK), alpha='1', scope='quiet-only' if quiet else 'complete homogeneous arity8/F1 class',
        all_source_sets_in_stratum=True, source_count=data['source_count'],
        fixed_rowset_holder_exponent=16, raw_weights_exclude_rs_normalizer=True,
        extra_rowset_factor=False, raw_regional_banks_built=0, regional_alpha=1,
        mode=mode, requested_controls=2 if mode == 'both' else 1,
        whole_word_surrogate_floor=FLOOR, surrogate='min(output_weight,8)',
        individual_row_floors=False, offdiagonal_empty_zeros_retained=False,
        alpha_after_partial_mask_sum=False, source_count_once=True,
        column_matrix_order='chronological', terminal_vector=[1, 1],
        moment_secants={'00': [5, 18], 'cross': [2, 18], '11': 'constant'},
        state_measure=measure, prefix_log2=prefix, parameter_fits=0,
        artifact=ARTIFACT, artifact_sha256=ARTIFACT_SHA, reviewed_source_sha256=SOURCE_PINS,
        compute_budget_seconds=max_seconds, directed=False, whole_code_certificate=False,
        complete_requested_controls=False, controls=[])
    try:
        if mode in ('scalar', 'both'):
            moments = old['log2_probability']+old['log2_nonzero_atom']*np.arange(GROUPS+1)
            moments[0] = 0.
            raw_quiet = quiet_mixture_coefficients(GROUPS, ARITY, moments, check_time=check_time)
            word = log_power(raw_quiet, COLUMNS, check_time)
            marker = adjacent_floor_marker(word, FLOOR)
            matrix, costs = column_matrix(weights, GROUPS, measure, arity=ARITY,
                                          log_marker=marker, check_time=check_time)
            power = log_matrix_power(matrix, COLUMNS)
            scalar = float(np.logaddexp2.reduce(power[0]))-FLOOR*marker
            row = dict(mode='scalar', log2_marker=marker, marker_log_as_exact_float=str(F(marker)),
                marker_rule='max(0,logcoeff1296-logcoeff1297) of frozen quiet moment word',
                adjacent_quiet_coefficients=[float(word[FLOOR]), float(word[FLOOR+1])],
                scalar_column_matrix_log2=matrix.tolist(), tail_upper_log2=scalar,
                log2_screen=prefix+scalar, **costs)
            result['controls'].append(row)
            if progress:
                progress(row)
        if mode in ('polynomial', 'both'):
            matrix, costs = column_matrix(weights, GROUPS, measure, arity=ARITY,
                                          check_time=check_time)
            word = matrix_poly_power(matrix, COLUMNS, check_time)
            tail = float(np.logaddexp2.reduce(word[0, :, FLOOR:].ravel()))
            row = dict(mode='polynomial', column_degree=matrix.shape[2]-1,
                word_degree=word.shape[2]-1, positive_tail_log2=tail,
                log2_screen=prefix+tail, **costs)
            result['controls'].append(row)
            if progress:
                progress(row)
        check_time()
        result['complete_requested_controls'] = True
    except TimeoutError as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('scalar', 'polynomial', 'both'), default='scalar')
    parser.add_argument('--quiet', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite an isolated control')
    answer = run(mode=args.mode, quiet=args.quiet, max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(complete=answer['complete_requested_controls'], elapsed_seconds=answer['elapsed_seconds'])))

"""Saved coupled-bank gate with a common live cone and separate entry moments.

One fixed normalized live mixture is used at every region. Compression
precedes the maximum over full-pair counts. Four prescribed positive
moment laws dominate the resulting entries; their chronological column
matrix polynomial is then raised to sixteen and its complete floor tail
is retained. No regional provider, parameter search, or old receipt is used.

The 01 and 10 atoms have distinct fixed secants. Their intercepts cover
ALL occupancies, including zero. The 00 atom is accepted only if its
computed probability is at most one. These are FLOAT screens only.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, expm1, log, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_coupled_state_compress import (FIXED_LOG_WEIGHTS, compress_bank,
    entrywise_curves, load_saved_bank)
from hd_two_sum_gf16_mean_column import raw_holder_weights
from hd_two_sum_gf16_state_envelope import elementary_logs, path_classes
from hd_two_sum_late_alpha_quiet import log_power
from hd_two_sum_matrix_moment import (entry_measure as old_entry_measure,
    load_saved as load_old_bank, matrix_poly_power)
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
GROUPS, ARITY, COLUMNS, FLOOR = 18, 8, 16, 1296
MARK = F(691, 700)
OLD_GATE = 'hd_two_sum_matrix_moment_gate_float.json'
OLD_GATE_SHA = '0465d355fd56c5c0cdfabac206a377ed74f62bf638b78af354bb4a147a9b5783'
SOURCE_PINS = {
    'hd_coupled_state_compress.py': '3cd4d2de91fa0bb1fa13b267d0157a1bd81f66c86f0f194bf2acd5f66e56d94f',
    'hd_coupled_state_gate.py': '974db96f33705f1d58f83fece7da5e1bd408e082946115bd9e865565eaa4efe0',
    'hd_two_sum_matrix_moment.py': 'bc192ebdafc9abec59f7977e660c5d4741737eb62cf479c920fdbb83a060ce67',
    'hd_two_sum_gf16_mean_column.py': 'fa6377002fa52fbfcb980ddd4b72d30df16b6c12c7c6db24eed54554e61a1593',
    'hd_two_sum_gf16_state_envelope.py': '5107b945ab18415191c48a3434a354c101d8fc5aeda66ec83a3fd2a5c3498af7',
    'hd_two_sum_late_alpha_quiet.py': '8a753036a9be44759206b5beb3baa044f28d80b265a9304da81bdde6e8f11ae3',
    'hd_two_sum_closed_kernel.py': '2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979',
}


class MomentFitError(ValueError):
    """The one requested moment family cannot cover the supplied curve."""


def _logs(values):
    result = np.asarray(values, dtype=float)
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Finite-or-zero logarithms required')
    return result


def entry_measure(curves, *, zero_endpoints=(5, 18), cross_endpoints=(2, 18)):
    """Prescribed 00,01,10 secants; independent all-occupancy intercepts."""
    curves = _logs(curves)
    if curves.ndim != 2 or curves.shape[0] != 4 or curves.shape[1] < 3:
        raise ValueError('Complete four-entry regional curves required')
    for endpoints in (zero_endpoints, cross_endpoints):
        if (len(endpoints) != 2 or any(type(i) is not int for i in endpoints)
                or not 1 <= endpoints[0] < endpoints[1] < curves.shape[1]):
            raise ValueError('Two ordered positive in-scope secant endpoints required')
    if curves[0, 0] > 0:
        raise MomentFitError('Requested 00 law fixes moment zero at one, below the supplied cap')
    slopes = []
    for row, endpoints in ((0, zero_endpoints), (1, cross_endpoints), (2, cross_endpoints)):
        left, right = endpoints
        if not np.isfinite(curves[row, [left, right]]).all():
            raise MomentFitError('Each prescribed secant needs two positive finite coefficients')
        slopes.append(float((curves[row, right]-curves[row, left])/(right-left)))
    degrees = np.arange(curves.shape[1])
    logp = float(np.max(curves[0, 1:]-slopes[0]*degrees[1:]))
    if not np.isfinite(logp) or logp > 0:
        raise MomentFitError('Prescribed 00 secant requires a nonzero-atom probability greater than one')
    scale01 = float(np.max(curves[1]-slopes[1]*degrees))
    scale10 = float(np.max(curves[2]-slopes[2]*degrees))
    scale11 = float(np.max(curves[3]))
    if not np.isfinite([scale01, scale10, scale11]).all():
        raise MomentFitError('Requested matrix measure requires finite positive live entry caps')
    bound = np.stack((logp+slopes[0]*degrees, scale01+slopes[1]*degrees,
                      scale10+slopes[2]*degrees, np.full(len(degrees), scale11)))
    bound[0, 0] = 0.
    if np.any(bound < curves-1e-9):
        raise ArithmeticError('Prescribed moment family misses an occupancy')
    return dict(log2_probability=logp, log2_atom00=slopes[0],
        log2_atom01=slopes[1], log2_atom10=slopes[2], log2_scale01=scale01,
        log2_scale10=scale10, log2_scale11=scale11,
        zero_secant_endpoints=list(zero_endpoints), cross_secant_endpoints=list(cross_endpoints),
        distinct_cross_slopes=True, all_occupancies_checked=True,
        zero_moment_preserved=True, minimum_slack=float(np.min(bound-curves)),
        contacts=[np.flatnonzero(row < 1e-9).tolist() for row in bound-curves],
        exact_float_parameters={key: str(F(value)) for key, value in
            (('probability_log', logp), ('atom00_log', slopes[0]), ('atom01_log', slopes[1]),
             ('atom10_log', slopes[2]), ('scale01_log', scale01),
             ('scale10_log', scale10), ('scale11_log', scale11))})


def from_shared_cross(measure):
    """Embed the pinned old law for a matched arithmetic control only."""
    result = dict(measure)
    result['log2_atom01'] = result['log2_atom10'] = result['log2_atom_cross']
    result['distinct_cross_slopes'] = False
    return result


def group_polynomial(weights, counts, active00, measure, arity):
    """Keep each directed cross atom in the labelled regional multiset."""
    weights = _logs(weights)
    if (len(counts) != 4 or any(type(n) is not int or n < 0 for n in counts)
            or weights.shape != (sum(counts)+1,) or type(active00) is not int
            or not 0 <= active00 <= counts[0] or type(arity) is not int
            or not 1 <= arity < len(weights)):
        raise ValueError('Complete weights and feasible ordered transition counts required')
    n00, n01, n10, n11 = counts
    markers = np.array([measure['log2_atom00']]*active00
        +[measure['log2_atom01']]*n01+[measure['log2_atom10']]*n10+[0.]*n11)
    if np.any(~np.isfinite(markers)):
        raise ValueError('Finite nonzero atom logarithms required')
    elementary = elementary_logs(markers)
    result = np.full(arity+1, -np.inf)
    for weight in range(len(elementary)):
        degree = min(weight, arity)
        result[degree] = np.logaddexp2(result[degree], weights[weight]+elementary[weight])
    return result


def column_matrix(weights, groups, measure, *, arity, check_time=lambda: None):
    """Positive chronological matrix polynomial, with separate 01/10 atoms.

    State-path multiplicity remains outside the group polynomial. Each
    cached key retains n01 and n10 separately. The type and row-set law is
    unchanged across all paths; raw Holder weights contain its union once.
    """
    weights = _logs(weights)
    size = len(weights)-1
    if (weights.ndim != 1 or type(groups) is not int or groups < 1
            or not 1 <= size <= 16 or type(arity) is not int or not 1 <= arity <= size):
        raise ValueError('Positive group and region counts and a valid arity required')
    logp = measure['log2_probability']
    if not np.isfinite(logp) or logp > 0:
        raise ValueError('00 atom probability must belong to (0,1]')
    if any(not np.isfinite(measure[key]) for key in ('log2_scale01', 'log2_scale10', 'log2_scale11')):
        raise ValueError('Finite positive entry scales required')
    logother = -np.inf if logp == 0 else log2(-expm1(logp*log(2)))
    result = np.full((2, 2, groups*arity+1), -np.inf)
    cached, terms = {}, 0
    classes = path_classes(size)
    for initial, final, counts, multiplicity in classes:
        check_time()
        n00, n01, n10, n11 = counts
        state_mass = (log2(multiplicity)+n01*measure['log2_scale01']
            +n10*measure['log2_scale10']+n11*measure['log2_scale11'])
        for active00 in range(n00+1):
            if logp == 0 and active00 != n00:
                continue
            atom_mass = (log2(comb(n00, active00))+active00*logp
                         +((n00-active00)*logother if active00 < n00 else 0.))
            key = active00, n01, n10, n11
            if key not in cached:
                polynomial = group_polynomial(weights, counts, active00, measure, arity)
                cached[key] = log_power(polynomial, groups, check_time)
            result[initial, final] = np.logaddexp2(result[initial, final],
                cached[key]+atom_mass+state_mass)
            terms += 1
    return result, dict(path_classes=len(classes), positive_path_atom_terms=terms,
        distinct_group_polynomials=len(cached), separate_cross_counts_in_cache=True)


def polynomial_control(weights, measure, prefix, *, check_time=lambda: None):
    matrix, costs = column_matrix(weights, GROUPS, measure, arity=ARITY, check_time=check_time)
    word = matrix_poly_power(matrix, COLUMNS, check_time)
    tail = float(np.logaddexp2.reduce(word[0, :, FLOOR:].ravel()))
    return dict(log2_screen=prefix+tail, positive_tail_log2=tail,
        column_degree=matrix.shape[2]-1, word_degree=word.shape[2]-1,
        terminal_vector=[1, 1], extra_terminal_factor_log2=0., **costs)


def _load_old_gate():
    raw = (DIRECTORY/OLD_GATE).read_bytes()
    if sha256(raw).hexdigest() != OLD_GATE_SHA:
        raise ValueError('Pinned old complete polynomial control changed')
    data = json.loads(raw)
    required = dict(active_groups=18, arity=8, source_full_pair_max=1,
        mark='691/700', complete_requested_controls=True, extra_rowset_factor=False,
        source_sha256=SOURCE_PINS['hd_two_sum_matrix_moment.py'], directed=False,
        whole_code_certificate=False)
    if any(data.get(key) != value for key, value in required.items()):
        raise ValueError('Exact old homogeneous r18/a8/F1 control required')
    rows = [row for row in data['controls'] if row.get('mode') == 'polynomial']
    if len(rows) != 1:
        raise ValueError('Exactly one old complete polynomial comparison required')
    return data, rows[0]


def run(*, max_seconds=30, progress=None):
    """One matched old control and one coupled control, using saved banks only."""
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit positive budget at most thirty seconds required')
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source changed: '+name)
    coupled_data, bank = load_saved_bank()
    old_data, old_curves = load_old_bank()
    old_gate, expected_old = _load_old_gate()
    if not coupled_data['source_count'] == old_data['source_count'] == old_gate['source_count']:
        raise ValueError('Matched controls must count the same exact source stratum')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Saved coupled moment control exhausted its budget')
    prefix = (log2(int(coupled_data['source_count']))-LIMIT*log2(float(MARK))
              -log2(255)-1152*log2(255))
    result = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=1,
        source_count=coupled_data['source_count'], mark=str(MARK), alpha=1,
        construction='proposed-g16-coupled-second-order-inner',
        scope='complete homogeneous r18/arity8/Fmax1 class only',
        raw_regional_banks_built=0, regional_alpha=1, requested_controls=2,
        fixed_live_log_weights=list(FIXED_LOG_WEIGHTS), moment_secants={'00': [5, 18],
            '01': [2, 18], '10': [2, 18], '11': 'constant'},
        fixed_rowsets_across_columns=True, raw_holder_exponent=16,
        extra_rowset_factor=False, source_count_once=True, projective_normalizations=1,
        rs_normalizer_exponent=1152, aggregate_surrogate_floor=FLOOR,
        all_outer_masks=True, all_boundary_state_paths=True,
        offdiagonal_empty_zeros_retained=False, parameter_fits=0,
        old_construction_receipts_transferred=False, prefix_log2=prefix,
        reviewed_source_sha256=SOURCE_PINS, old_control_artifact=OLD_GATE,
        old_control_artifact_sha256=OLD_GATE_SHA, compute_budget_seconds=max_seconds,
        directed=False, whole_code_certificate=False, complete_requested_controls=False,
        controls=[])
    try:
        check_time()
        weights = raw_holder_weights()
        old_atoms = old_data['two_atom_envelope']
        old_measure = from_shared_cross(old_entry_measure(old_curves,
            old_atoms['log2_probability'], old_atoms['log2_nonzero_atom']))
        old = polynomial_control(weights, old_measure, prefix, check_time=check_time)
        old.update(name='matched_old_shared_cross_control', state_measure=old_measure,
            saved_log2_screen=expected_old['log2_screen'],
            reproduction_difference_log2=old['log2_screen']-expected_old['log2_screen'])
        if abs(old['reproduction_difference_log2']) > 1e-8:
            raise ArithmeticError('Separate-cross engine failed to reproduce the pinned old control')
        result['controls'].append(old)
        if progress:
            progress(old)
        compressed = compress_bank(bank, FIXED_LOG_WEIGHTS)
        curves = entrywise_curves(compressed['regions'])
        result.update(compression=compressed['metadata'],
            compressed_entry_log2_curves=[[None if value == -np.inf else float(value)
                                          for value in row] for row in curves])
        measure = entry_measure(curves)
        coupled = polynomial_control(weights, measure, prefix, check_time=check_time)
        coupled.update(name='coupled_distinct_entry_moments', state_measure=measure)
        result['controls'].append(coupled)
        check_time()
        result['complete_requested_controls'] = True
        if progress:
            progress(coupled)
    except (TimeoutError, MomentFitError) as error:
        result['incomplete_reason'] = str(error)
        result['incomplete_kind'] = type(error).__name__
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
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
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')

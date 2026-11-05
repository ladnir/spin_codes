"""FLOAT full-branch comparison with fixed row sets and positive Young bounds.

For a fixed arity-eight row set A and its complement O, split the exact
column law K_A=Kq+Kf+Kh. Quiet masks have coefficient Z2(b) and mark X^b.
Full/hole masks have inside coefficients f_b/h_b and the surrogate mark
X^8, not X^b. For positive lambdas summing to one, K_A^r is at most
sum_i lambda_i^(1-r) K_i^r. The hole sum is bounded by
|O|^(r-1) sum_{j in O} K_h,j^r. Both inequalities are pointwise in the
same positive regional path/atom measure.

Outside A, full terms use E(r); a hole term uses E(0) at exactly one
outside position. These are evaluations of that SAME moment measure,
not raw regional caps. An exact degree-one flag sums the hole positions.
The three positive column matrices are added before the ordered word
power. Whole-word group Minkowski then sums fixed A once. No provider,
parameter fit, or production calculation runs on import. This is a FLOAT
comparison for one homogeneous class, not a whole-code certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_coupled_state_moment_gate import entry_measure
from hd_post_refresh_compress import compress_bank, entrywise_curves
from hd_post_refresh_moment_gate import (CONSTRUCTION, MARK, SOURCE_COUNT,
    SOURCE_PINS, check_source_pins, load_saved_bank)
from hd_post_refresh_quiet_rows import (combine_row_bounds, log_matmul,
    log_matrix_power, moment_atoms, quiet_column_matrices, quiet_group_values,
    rowset_gaps, transition_geometry)
from hd_two_sum_closed_kernel import branch_counts
from hd_xor_multigroup import LIMIT


QUIET_SHA = '75e38819d225e9771ade1ed90509b7d4582a2637b24b6522f5b62b0200eb49e7'
LAMBDAS = (F(99, 100), F(1, 200), F(1, 200))  # quiet, full, hole


def _logs(values):
    values = np.asarray(values, dtype=float)
    if np.any(np.isnan(values)) or np.any(np.isposinf(values)):
        raise ValueError('Finite-or-zero log coefficients required')
    return values


def moment_evaluation(measure, occupancy):
    """Evaluate the common positive matrix measure, including its empty slack."""
    moment_atoms(measure)
    if type(occupancy) is not int or occupancy < 0:
        raise ValueError('Nonnegative integer occupancy required')
    zero = (0. if occupancy == 0 else
            measure['log2_probability']+occupancy*measure['log2_atom00'])
    return np.array([[zero, measure['log2_scale01']+occupancy*measure['log2_atom01']],
        [measure['log2_scale10']+occupancy*measure['log2_atom10'], measure['log2_scale11']]])


def gap_powers(base, maximum, *, hole=None):
    """Return [Y^0,Y^1](base+Y*hole)^n in chronological matrix order.

    If hole is omitted, retain only degree zero. Degree one is an exact
    coefficient, not an OR flag: terms with two hole positions are dropped.
    """
    base = _logs(base)
    if base.shape != (2, 2) or type(maximum) is not int or not 0 <= maximum <= 16:
        raise ValueError('Two-state gap matrix and maximum0..16 required')
    if hole is not None:
        hole = _logs(hole)
        if hole.shape != (2, 2):
            raise ValueError('Two-state hole matrix required')
    flags = 1 if hole is None else 2
    result = np.full((maximum+1, flags, 2, 2), -np.inf)
    result[0, 0, 0, 0] = result[0, 0, 1, 1] = 0.
    for length in range(1, maximum+1):
        result[length, 0] = log_matmul(result[length-1, 0], base)
        if flags == 2:
            result[length, 1] = np.logaddexp2(
                log_matmul(result[length-1, 1], base),
                log_matmul(result[length-1, 0], hole))
    return result


def marked_column_matrices(log_counts, groups, measure, *, size=16,
                           one_hole=False, batch_size=128, check_time=lambda: None):
    """Full or one-hole branch for every A, with no activity or Young factors.

    The marked-A count vector is f_b or h_b. Its elementary polynomial has
    NO X^b factor. E(groups) occupies every gap site, except the exact one
    E(0) chosen by the hole flag. Five marked-edge counts share the frozen
    quiet helper's atom geometry. Flags and matrices preserve chronology.
    """
    log_counts = _logs(log_counts)
    if log_counts.ndim != 1 or type(one_hole) is not bool:
        raise ValueError('A count vector and boolean one-hole selection required')
    rows, gaps = rowset_gaps(size, len(log_counts)-1)
    outside = size-len(log_counts)+1
    if (type(groups) is not int or groups < 1 or type(batch_size) is not int
            or not 1 <= batch_size <= 512 or outside < 1):
        raise ValueError('Positive groups/outside count and batch size1..512 required')
    masses, _ = moment_atoms(measure)
    values = groups*quiet_group_values(log_counts, measure, 0.)
    counts, layers = transition_geometry(len(log_counts)-1)
    active = moment_evaluation(measure, groups)
    empty = moment_evaluation(measure, 0)
    powers = gap_powers(active, outside, hole=empty if one_hole else None)
    flags = 2 if one_hole else 1
    result = np.full((len(rows), 2, 2), -np.inf)
    endpoints = ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1))
    flag_products = ((0, 0, 0), (1, 0, 1), (0, 1, 1)) if one_hole else ((0, 0, 0),)
    for start in range(0, len(rows), batch_size):
        check_time()
        stop = min(len(rows), start+batch_size)
        block = gaps[start:stop]
        current = powers[block[:, 0]].transpose(1, 0, 2, 3)[:, None]
        for position, (next_size, targets) in enumerate(layers):
            check_time()
            following = np.full((flags, next_size, stop-start, 2, 2), -np.inf)
            after = powers[block[:, position+1]].transpose(1, 0, 2, 3)
            for label, ((incoming, outgoing), mass) in enumerate(zip(endpoints, masses)):
                if mass == -np.inf:
                    continue
                target = targets[label]
                for old_flag, gap_flag, new_flag in flag_products:
                    term = (current[old_flag, ..., incoming, None]
                        +after[gap_flag, None, :, None, outgoing, :]+mass)
                    following[new_flag, target] = np.logaddexp2(following[new_flag, target], term)
            current = following
        result[start:stop] = np.logaddexp2.reduce(
            current[int(one_hole)]+values[:, None, None, None], axis=0)
    return dict(rowsets=rows, column_log2_matrices=result,
        gap_log2_matrix=active, hole_log2_matrix=empty if one_hole else None,
        one_hole=one_hole, costs=dict(rowset_count=len(rows), final_atom_count_classes=len(counts),
            batch_size=batch_size, hole_flag_degree=int(one_hole),
            maximum_dp_layer_bytes=flags*len(counts)*min(batch_size, len(rows))*4*8))


def checked_lambdas(lambdas):
    values = tuple(F(value) for value in lambdas)
    if len(values) != 3 or any(value <= 0 for value in values) or sum(values) != 1:
        raise ValueError('Three positive exact Young weights summing to one required')
    return values


def combine_columns(quiet, full, holes, groups, arity, outside, log2_x, *, lambdas=LAMBDAS):
    """Add positive Young branch matrices BEFORE the whole-word power."""
    quiet, full, holes = map(_logs, (quiet, full, holes))
    if (quiet.shape != full.shape or quiet.shape != holes.shape or quiet.shape[-2:] != (2, 2)
            or any(type(value) is not int or value < 1 for value in (groups, arity, outside))
            or not np.isfinite(log2_x) or log2_x < 0):
        raise ValueError('Matched branch matrices and positive group/arity/outside scope required')
    weights = checked_lambdas(lambdas)
    shifts = np.array([(1-groups)*log2(weights[0]),
        (1-groups)*log2(weights[1])+arity*groups*log2_x,
        (1-groups)*log2(weights[2])+(groups-1)*log2(outside)+arity*groups*log2_x])
    result = np.logaddexp2(np.logaddexp2(quiet+shifts[0], full+shifts[1]), holes+shifts[2])
    return result, shifts


def full_row_bounds(branch_logs, groups, measure, log2_x, *, size=16, columns=16,
                    lambdas=LAMBDAS, batch_size=128, check_time=lambda: None, progress=None):
    """Complete Young comparison, then fixed-A whole-word group Minkowski."""
    logs = _logs(branch_logs)
    if (logs.ndim != 2 or logs.shape[0] != 3 or type(columns) is not int or columns < 1
            or not np.isfinite(log2_x) or log2_x < 0):
        raise ValueError('Quiet/full/hole counts, positive columns and X>=1 required')
    checked_lambdas(lambdas)
    arity, outside = logs.shape[1]-1, size-logs.shape[1]+1
    quiet = quiet_column_matrices(logs[0], groups, measure, log2_x, size=size,
        batch_size=batch_size, check_time=check_time)
    if progress:
        progress('quiet_complete')
    full = marked_column_matrices(logs[1], groups, measure, size=size, batch_size=batch_size,
        check_time=check_time)
    if progress:
        progress('full_complete')
    holes = marked_column_matrices(logs[2], groups, measure, size=size, one_hole=True,
        batch_size=batch_size, check_time=check_time)
    if progress:
        progress('hole_complete')
    combined, shifts = combine_columns(quiet['column_log2_matrices'], full['column_log2_matrices'],
        holes['column_log2_matrices'], groups, arity, outside, log2_x, lambdas=lambdas)
    check_time()
    word = log_matrix_power(combined, columns)
    values = np.logaddexp2.reduce(word[:, 0, :], axis=1)
    return dict(rowsets=quiet['rowsets'], per_rowset_word_log2=values,
        combined_word_log2=combine_row_bounds(values, groups), column_log2_matrices=combined,
        branch_column_log2_matrices={key: value['column_log2_matrices'] for key, value in
            (('quiet', quiet), ('full', full), ('holes', holes))},
        costs={key: value['costs'] for key, value in (('quiet', quiet), ('full', full), ('holes', holes))},
        young_log2_scales=shifts, terminal_vector=[1, 1], extra_rowset_factor=False)


def run(bank_path, bank_sha256, *, log2_x, max_seconds=30, batch_size=128, progress=None):
    """One frozen-X saved-bank comparison; no numerical provider or fit."""
    if (type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30
            or not np.isfinite(log2_x) or log2_x < 0
            or type(batch_size) is not int or not 1 <= batch_size <= 512):
        raise ValueError('X>=1, batch1..512 and a positive budget at most30seconds required')
    directory = Path(__file__).resolve().parent
    check_source_pins()
    if sha256((directory/'hd_post_refresh_quiet_rows.py').read_bytes()).hexdigest() != QUIET_SHA:
        raise ValueError('Frozen quiet helper changed')
    _, bank = load_saved_bank(bank_path, bank_sha256)
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Full fixed-row-set comparison exhausted its budget')
    result = dict(construction=CONSTRUCTION, active_groups=18, arity=8, columns=16,
        scope='complete homogeneous r18/arity8/Fmax exactly1 class only',
        source_count=str(SOURCE_COUNT), mark=str(MARK), alpha=1, regional_alpha=1,
        log2_activity_marker=float(log2_x), aggregate_surrogate_floor=1296,
        rs_normalizer_exponent=1152, projective_normalizations=1,
        rowsets_fixed_across_columns=True, rowset_count=comb(16, 8),
        rowset_count_inside_wholeword_Minkowski=True, extra_rowset_factor=False,
        common_positive_path_atom_measure=True, exact_empty_replacement=False,
        young_weights=[str(value) for value in LAMBDAS], young_weights_selected_before_data=True,
        one_hole_summed_exactly_once=True, hole_jensen_log2_factor=17*log2(8),
        branch_sum_before_word_power=True, empty_offdiagonal_relaxation_retained=True,
        all_boundary_state_paths=True, all_outer_branches=True, terminal_vector=[1, 1],
        raw_regional_banks_built=0, parameter_fits=0, raw_bank_artifact=str(bank_path),
        raw_bank_artifact_sha256=bank_sha256, reviewed_source_sha256=dict(SOURCE_PINS,
            **{'hd_post_refresh_quiet_rows.py': QUIET_SHA}),
        compute_budget_seconds=max_seconds, completed_stages=[], complete_requested_control=False,
        directed=False, whole_code_certificate=False)
    def stage(name):
        result['completed_stages'].append(name)
        if progress:
            progress(dict(stage=name, elapsed_seconds=perf_counter()-started))
    try:
        compressed = compress_bank(bank)
        measure = entry_measure(entrywise_curves(compressed['regions']))
        result.update(state_measure=measure, compression=compressed['metadata'])
        counts = [branch_counts(8, b, 16) for b in range(9)]
        logs = np.array([[log2(row[index]) if row[index] else -np.inf
                          for row in counts] for index in (0, 2, 1)])
        answer = full_row_bounds(logs, 18, measure, log2_x, batch_size=batch_size,
            check_time=check_time, progress=stage)
        prefix = (log2(SOURCE_COUNT)-LIMIT*log2(float(MARK))-log2(255)
                  -1152*log2(255)-1296*log2_x)
        result.update(prefix_log2=prefix, combined_word_log2=answer['combined_word_log2'],
            log2_screen=prefix+answer['combined_word_log2'], costs=answer['costs'],
            young_log2_scales=answer['young_log2_scales'].tolist(),
            per_rowset=[dict(positions=list(row), word_log2=float(value), root_log2=float(value/18))
                for row, value in zip(answer['rowsets'], answer['per_rowset_word_log2'])])
        check_time()
        result['complete_requested_control'] = True
    except TimeoutError as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bank', type=Path, required=True)
    parser.add_argument('--bank-sha256', required=True)
    parser.add_argument('--log2-x', type=float, required=True)
    parser.add_argument('--batch-size', type=int, default=128)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen comparison')
    answer = run(args.bank, args.bank_sha256, log2_x=args.log2_x,
        max_seconds=args.max_seconds, batch_size=args.batch_size,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')

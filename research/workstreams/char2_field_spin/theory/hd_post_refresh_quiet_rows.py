"""FLOAT quiet-only comparison retaining each row set across the whole word.

Fix a positive activity marker X and the common regional moment measure.
For a row set A, K_A(t,X)=sum_b Z2(b)e_b(t_A)X^b counts its quiet masks.
Q_A is the ordered column matrix with K_A**groups inside the SAME positive
path/atom measure. Write f_A=e0 Q_A**columns 1. Holder on that whole-word
measure gives the all-row-set comparison (sum_A f_A**(1/groups))**groups.
The row sets are fixed, not resampled between columns. Their multiplicity
is already in this sum; there is no additional binomial factor.

Empty positions use S(0) from the COMMON moment measure, including its
artificial off-diagonal entries. Replacing these gaps by a smaller exact
empty matrix separately for each A would invalidate the mixed-A bound.
Only quiet outer branches are included. This is not a full-class bound,
a directed calculation, or a distance certificate. No provider runs here.
"""
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import comb, expm1, log, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_coupled_state_moment_gate import entry_measure
from hd_post_refresh_compress import compress_bank, entrywise_curves
from hd_post_refresh_moment_gate import (CONSTRUCTION, MARK, SOURCE_COUNT,
    check_source_pins, load_saved_bank)
from hd_two_sum_closed_kernel import zero_two_moments
from hd_two_sum_gf16_state_envelope import elementary_logs
from hd_xor_multigroup import LIMIT


def _logs(values):
    result = np.asarray(values, dtype=float)
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Finite-or-zero logarithms required')
    return result


def log_matmul(left, right):
    """Ordered positive 2x2 product, permitting leading batch dimensions."""
    left, right = _logs(left), _logs(right)
    if left.shape[-2:] != (2, 2) or right.shape[-2:] != (2, 2):
        raise ValueError('Two-state matrices required')
    return np.logaddexp2(left[..., :, 0, None]+right[..., None, 0, :],
                        left[..., :, 1, None]+right[..., None, 1, :])


def log_matrix_power(matrix, exponent):
    matrix = _logs(matrix)
    if matrix.shape[-2:] != (2, 2) or type(exponent) is not int or exponent < 0:
        raise ValueError('Two-state matrix and nonnegative integer exponent required')
    answer = np.full(matrix.shape, -np.inf)
    answer[..., 0, 0] = answer[..., 1, 1] = 0.
    while exponent:
        if exponent & 1:
            answer = log_matmul(answer, matrix)
        exponent >>= 1
        if exponent:
            matrix = log_matmul(matrix, matrix)
    return answer


@lru_cache(None)
def rowset_gaps(size, arity):
    """All labelled A and their ordered gaps, with no orbit reduction."""
    if (type(size) is not int or type(arity) is not int
            or not 1 <= arity <= min(size, 8) or not size <= 16):
        raise ValueError('Scope 1<=arity<=8 and arity<=size<=16 required')
    rows = tuple(combinations(range(size), arity))
    gaps = np.array([[row[0], *[row[j]-row[j-1]-1 for j in range(1, arity)],
                     size-row[-1]-1] for row in rows], dtype=np.int8)
    gaps.setflags(write=False)
    return rows, gaps


@lru_cache(None)
def transition_geometry(arity):
    """Five marked-edge counts: 00/atom0, 00/atomB, 01, 10, 11."""
    if type(arity) is not int or not 1 <= arity <= 8:
        raise ValueError('One to eight marked positions required')
    counts = ((0, 0, 0, 0, 0),)
    layers = []
    for _ in range(arity):
        following = sorted({tuple(n+int(j == label) for j, n in enumerate(row))
                            for row in counts for label in range(5)})
        index = {row: j for j, row in enumerate(following)}
        targets = tuple(np.array([index[tuple(n+int(j == label)
            for j, n in enumerate(row))] for row in counts], dtype=np.int32)
            for label in range(5))
        for target in targets:
            target.setflags(write=False)
        layers.append((len(following), targets))
        counts = tuple(following)
    return counts, tuple(layers)


def moment_atoms(measure):
    keys = ('log2_probability', 'log2_atom00', 'log2_atom01', 'log2_atom10',
            'log2_scale01', 'log2_scale10', 'log2_scale11')
    if any(key not in measure or not np.isfinite(measure[key]) for key in keys):
        raise ValueError('Seven finite moment parameters required')
    logp = measure['log2_probability']
    if logp > 0:
        raise ValueError('The 00 nonzero-atom probability must not exceed one')
    other = -np.inf if logp == 0 else log2(-expm1(logp*log(2)))
    masses = (other, logp, measure['log2_scale01'], measure['log2_scale10'],
              measure['log2_scale11'])
    empty = np.array([[0., masses[2]], [masses[3], masses[4]]])
    return masses, empty


def quiet_group_values(log_counts, measure, log2_x):
    """One log K_A for each final atom-count tuple; no row-set factor."""
    log_counts = _logs(log_counts)
    if (log_counts.ndim != 1 or not 2 <= len(log_counts) <= 9
            or not np.isfinite(log2_x) or log2_x < 0):
        raise ValueError('Quiet counts through arity<=8 and a finite X>=1 required')
    moment_atoms(measure)
    counts, _ = transition_geometry(len(log_counts)-1)
    answer = []
    for _, n00, n01, n10, n11 in counts:
        markers = ([measure['log2_atom00']]*n00+[measure['log2_atom01']]*n01
                   +[measure['log2_atom10']]*n10+[0.]*n11)
        elementary = elementary_logs(np.asarray(markers))
        answer.append(float(np.logaddexp2.reduce(log_counts[:len(elementary)]
            +elementary+log2_x*np.arange(len(elementary)))))
    return np.asarray(answer)


def quiet_column_matrices(log_counts, groups, measure, log2_x, *, size=16,
                          batch_size=512, check_time=lambda: None):
    """Return every log Q_A using bounded-memory vectorized path-count DP.

    The gap matrix is S(0), never an A-dependent replacement by exact R(0).
    Each DP edge preserves its direction. Only the symmetric quiet mask
    polynomial is grouped by atom counts. All group factors precede the
    column/path summation, and all column products remain ordered.
    """
    log_counts = _logs(log_counts)
    if log_counts.ndim != 1:
        raise ValueError('One quiet count per support weight required')
    arity = len(log_counts)-1
    rows, gaps = rowset_gaps(size, arity)
    if (type(groups) is not int or groups < 1 or type(batch_size) is not int
            or not 1 <= batch_size <= 4096):
        raise ValueError('Positive group count and explicit batch size1..4096 required')
    masses, empty = moment_atoms(measure)
    values = groups*quiet_group_values(log_counts, measure, log2_x)
    counts, layers = transition_geometry(arity)
    gap_powers = np.stack([log_matrix_power(empty, n) for n in range(size-arity+1)])
    result = np.full((len(rows), 2, 2), -np.inf)
    endpoints = ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1))
    for start in range(0, len(rows), batch_size):
        check_time()
        stop = min(len(rows), start+batch_size)
        block = gaps[start:stop]
        current = gap_powers[block[:, 0]][None, ...]
        for position, (next_size, targets) in enumerate(layers):
            check_time()
            following = np.full((next_size, stop-start, 2, 2), -np.inf)
            after = gap_powers[block[:, position+1]]
            for label, ((incoming, outgoing), mass) in enumerate(zip(endpoints, masses)):
                if mass == -np.inf:
                    continue
                # For one label, target indices are injective: no scatter sum is lost.
                term = (current[..., incoming, None]
                        +after[None, :, None, outgoing, :]+mass)
                target = targets[label]
                following[target] = np.logaddexp2(following[target], term)
            current = following
        result[start:stop] = np.logaddexp2.reduce(current+values[:, None, None, None], axis=0)
    return dict(rowsets=rows, column_log2_matrices=result,
        common_empty_log2_matrix=empty,
        costs=dict(rowset_count=len(rows), final_atom_count_classes=len(counts),
            batch_size=batch_size, maximum_dp_layer_bytes=len(counts)*min(batch_size, len(rows))*4*8),
        exact_empty_replacement=False, common_positive_measure=True)


def combine_row_bounds(log_word_values, groups):
    """Whole-word group Holder/Minkowski; row sets are summed exactly once."""
    values = _logs(log_word_values)
    if values.ndim != 1 or not len(values) or type(groups) is not int or groups < 1:
        raise ValueError('Per-row-set word values and a positive group count required')
    return float(groups*np.logaddexp2.reduce(values/groups))


def quiet_row_bounds(log_counts, groups, measure, log2_x, *, size=16, columns=16,
                     batch_size=512, check_time=lambda: None):
    if type(columns) is not int or columns < 1:
        raise ValueError('Positive whole-word column count required')
    answer = quiet_column_matrices(log_counts, groups, measure, log2_x, size=size,
        batch_size=batch_size, check_time=check_time)
    check_time()
    word = log_matrix_power(answer['column_log2_matrices'], columns)
    values = np.logaddexp2.reduce(word[:, 0, :], axis=1)
    answer.update(per_rowset_word_log2=values, combined_word_log2=combine_row_bounds(values, groups),
        columns=columns, terminal_vector=[1, 1], extra_rowset_factor=False)
    return answer


def run(bank_path, bank_sha256, *, log2_x, max_seconds=30, batch_size=512):
    """One saved-bank, fixed-X quiet comparison; no provider or optimizer."""
    if (type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30
            or not np.isfinite(log2_x) or log2_x < 0):
        raise ValueError('One finite X>=1 and an explicit budget at most30seconds required')
    check_source_pins()
    _, bank = load_saved_bank(bank_path, bank_sha256)
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Quiet fixed-row-set comparison exhausted its budget')
    result = dict(construction=CONSTRUCTION, active_groups=18, arity=8, columns=16,
        scope='all-quiet outer component of homogeneous r18/arity8/Fmax exactly1 only',
        source_count=str(SOURCE_COUNT), mark=str(MARK), alpha=1, regional_alpha=1,
        log2_activity_marker=float(log2_x), aggregate_input_floor=1296,
        rs_normalizer_exponent=1152, projective_normalizations=1,
        rowsets_fixed_across_columns=True, rowset_count=comb(16, 8),
        rowset_count_inside_wholeword_Minkowski=True, extra_rowset_factor=False,
        common_positive_path_atom_measure=True, exact_empty_replacement=False,
        empty_offdiagonal_relaxation_retained=True, all_boundary_state_paths=True,
        all_outer_branches=False, terminal_vector=[1, 1], raw_regional_banks_built=0,
        parameter_fits=0, raw_bank_artifact=str(bank_path), raw_bank_artifact_sha256=bank_sha256,
        compute_budget_seconds=max_seconds, complete_requested_control=False,
        directed=False, whole_code_certificate=False)
    try:
        compressed = compress_bank(bank)
        curves = entrywise_curves(compressed['regions'])
        measure = entry_measure(curves)
        result.update(state_measure=measure, compression=compressed['metadata'])
        quiet = np.array([log2(n) if n else -np.inf
                          for n in (zero_two_moments(b) for b in range(9))])
        answer = quiet_row_bounds(quiet, 18, measure, log2_x, batch_size=batch_size,
                                  check_time=check_time)
        prefix = (log2(SOURCE_COUNT)-LIMIT*log2(float(MARK))-log2(255)
                  -1152*log2(255)-1296*log2_x)
        result.update(prefix_log2=prefix, combined_word_log2=answer['combined_word_log2'],
            log2_screen=prefix+answer['combined_word_log2'], costs=answer['costs'],
            common_empty_log2_matrix=answer['common_empty_log2_matrix'].tolist(),
            per_rowset=[dict(positions=list(row), word_log2=float(value),
                            root_log2=float(value/18)) for row, value in
                        zip(answer['rowsets'], answer['per_rowset_word_log2'])])
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
    parser.add_argument('--batch-size', type=int, default=512)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen comparison')
    answer = run(args.bank, args.bank_sha256, log2_x=args.log2_x,
                 max_seconds=args.max_seconds, batch_size=args.batch_size)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')

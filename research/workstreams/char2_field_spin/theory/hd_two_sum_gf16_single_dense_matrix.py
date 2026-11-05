"""FLOAT r1 dense gate retaining full2x2 matrices through all columns.

The common-potential control is preserved in hd_two_sum_gf16_single_dense.
Here Q_w sums ordered regional matrix products for masks of weightw.
Fixed-arity Holder weights form matrix coefficients B0 and B1; a positive
matrix polynomial computes(B0+Y B1)^16. Start in state0 and sum both final
states, retaining Y>=9. There is no potential or endpoint-potential factor.
Only two previously used r1 marks are tested, with no new optimization.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_gf16_arity_dense import arity_weights
from hd_two_sum_gf16_full_dense import _check_deadline
from hd_two_sum_gf16_single_dense import single_source_tables
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


MARKS = (('old-r1-inherited', F('499451537/500000000')),
         ('old-r1-rational', F('1199/1200')))


def log_matmul(left, right):
    """Positive2x2 matrix multiplication, broadcasting leading dimensions."""
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if (left.shape[-2:] != (2, 2) or right.shape[-2:] != (2, 2)
            or any(np.any(np.isnan(a)) or np.any(np.isposinf(a)) for a in (left, right))):
        raise ValueError('Finite-or-zero2x2 log matrices required')
    return np.logaddexp2(left[..., :, 0, None]+right[..., None, 0, :],
                        left[..., :, 1, None]+right[..., None, 1, :])


def weight_matrices(regions):
    """Q[w]=sum_|T|=w R0(1_T(0))...R_(L-1)(1_T(L-1))."""
    tables = np.asarray(regions, dtype=float)
    if (tables.ndim != 4 or tables.shape[1:] != (2, 2, 2) or not 1 <= len(tables) <= 16
            or np.any(np.isnan(tables)) or np.any(np.isposinf(tables))):
        raise ValueError('One to16 ordered occupancy0/1 regional log matrices required')
    values = np.array([[[0., -np.inf], [-np.inf, 0.]]])
    for table in tables:
        following = np.full((len(values)+1, 2, 2), -np.inf)
        for bit in (0, 1):
            target = following[bit:bit+len(values)]
            np.logaddexp2(target, log_matmul(values, table[bit]), out=target)
        values = following
    return values


def matrix_branches(matrices, weights):
    """Apply per-mask Holder weights before splitting empty/nonempty columns."""
    matrices, weights = np.asarray(matrices, dtype=float), np.asarray(weights, dtype=float)
    if (weights.ndim != 1 or len(weights) < 2 or matrices.shape != (len(weights), 2, 2)
            or any(np.any(np.isnan(a)) or np.any(np.isposinf(a)) for a in (matrices, weights))):
        raise ValueError('Complete matrix weight coefficients and finite-or-zero Holder weights required')
    return np.array([matrices[0]+weights[0],
                     np.logaddexp2.reduce(matrices[1:]+weights[1:, None, None], axis=0)])


def matrix_tail(branches, columns=16, y_floor=9):
    """Exact positive matrix-coefficient recurrence in FLOAT log arithmetic."""
    branches = np.asarray(branches, dtype=float)
    if (branches.shape != (2, 2, 2) or type(columns) is not int or columns < 1
            or type(y_floor) is not int or not 0 <= y_floor <= columns
            or np.any(np.isnan(branches)) or np.any(np.isposinf(branches))):
        raise ValueError('Two matrix coefficients and valid column/floor bounds required')
    values = np.array([[[0., -np.inf], [-np.inf, 0.]]])
    for _ in range(columns):
        following = np.full((len(values)+1, 2, 2), -np.inf)
        for bit in (0, 1):
            target = following[bit:bit+len(values)]
            np.logaddexp2(target, log_matmul(values, branches[bit]), out=target)
        values = following
    # The initial state is0; the terminal row-vector action sums both states.
    degrees = np.logaddexp2.reduce(values[y_floor:, 0, :], axis=1)
    return float(np.logaddexp2.reduce(degrees)), [
        dict(nonempty_dense_columns=degree, log2_mass=float(value))
        for degree, value in enumerate(degrees, y_floor)]


def run(*, max_seconds=20, progress=None):
    if (isinstance(max_seconds, bool) or not isinstance(max_seconds, (int, float))
            or not np.isfinite(max_seconds) or not 0 < max_seconds <= 20):
        raise ValueError('A positive compute budget at most20 seconds is required')
    started, deadline = perf_counter(), perf_counter()+max_seconds
    result = dict(construction=CONSTRUCTION, labels=list(LABELS), field_order=256,
        field_polynomial='0x11b', active_groups=1, cohort_counts=[0, 0, 1],
        all_source_supports=True, source_supports=8192, source_colorings=1,
        source_full_pairs=0, forbidden_weight=LIMIT, total_columns=16,
        constituent_dimension=8, sparse_groups=0, dense_nonempty_column_floor=9,
        holder_exponent=16, arities=list(range(3, 17)), fixed_dense_row_sets=65399,
        arities_kept_separate=True, dense_rowset_multiplicity_inside_holder=True,
        extra_dense_rowset_factor=False, fixed_rows_across_columns=True,
        fixed_source_geometry=True, pair_swaps=True, ordered_regional_matrix_products=True,
        ordered_column_matrix_polynomial=True, initial_state=0, terminal_vector=[1, 1],
        common_potential_used=False, endpoint_potential_factors=0,
        projective_normalizations=1, alpha=1, marks_requested=2, fitted_new_marks=0,
        mark_selection_per_arity=True, extra_mark_bank_factor=False,
        directed=False, whole_code_certificate=False, max_seconds=max_seconds, banks=[])
    try:
        _check_deadline(deadline)
        weights = {a: arity_weights(a) for a in range(3, 17)}
        for name, mark in MARKS:
            _check_deadline(deadline)
            raw = pair_regions_screen(float(mark), 1)
            matrices = weight_matrices(single_source_tables(raw))
            prefix = log2(8192)-log2(255)-LIMIT*log2(float(mark))
            rows = []
            for arity in range(3, 17):
                _check_deadline(deadline)
                tail, degrees = matrix_tail(matrix_branches(matrices, weights[arity][0]))
                peak = max(degrees, key=lambda row: row['log2_mass'])
                rows.append(dict(active_constituents=arity, fixed_row_sets=comb(16, arity),
                    exact_moments_sha256=weights[arity][1], failure_log2_screen=prefix+tail,
                    largest_dense_degree=peak['nonempty_dense_columns']))
            bank = dict(witness_label=name, mark=str(mark), entries=rows,
                failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in rows])))
            result['banks'].append(bank)
            if progress:
                progress({key: value for key, value in bank.items() if key != 'entries'})
        _check_deadline(deadline)
        selected = []
        for index, arity in enumerate(range(3, 17)):
            choice = min(range(2), key=lambda j: result['banks'][j]['entries'][index]['failure_log2_screen'])
            bank = result['banks'][choice]
            selected.append(dict(bank['entries'][index], witness_label=bank['witness_label'], mark=bank['mark']))
        result.update(status='complete-r1-dense-matrix-screen', complete_selected_cohort=True,
            entries=selected, failing_arities=[row['active_constituents'] for row in selected
                                             if row['failure_log2_screen'] > -40],
            failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in selected])))
    except TimeoutError:
        result.update(status='compute-budget-exceeded', complete_selected_cohort=False)
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=20)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded FLOAT artifact')
    result = run(max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, 'hd_two_sum_gf16_single_dense.py',
                     'hd_two_sum_closed_kernel.py', 'hd_two_sum_gf16_arity_dense.py', 'hd_xor_pair_route.py')}
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'banks'}, allow_nan=False))

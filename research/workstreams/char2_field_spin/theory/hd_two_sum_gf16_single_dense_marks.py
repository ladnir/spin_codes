"""Four fixed rational r1 marks using the unchanged full-matrix interface.

This is a bounded diagnostic bank, not an optimizer or a directed proof.
The prior common-potential and two-old-mark artifacts remain unchanged.
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
from hd_two_sum_gf16_single_dense_matrix import weight_matrices, matrix_branches, matrix_tail
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


MARKS = tuple(map(F, ('1999/2000', '2499/2500', '9997/10000', '4999/5000')))


def run(*, max_seconds=20, progress=None):
    if (isinstance(max_seconds, bool) or not isinstance(max_seconds, (int, float))
            or not np.isfinite(max_seconds) or not 0 < max_seconds <= 20):
        raise ValueError('A positive compute budget at most20 seconds is required')
    started = perf_counter()
    deadline = started+max_seconds
    result = dict(construction=CONSTRUCTION, labels=list(LABELS), field_order=256,
        field_polynomial='0x11b', active_groups=1, cohort_counts=[0, 0, 1],
        all_source_supports=True, source_supports=8192, source_colorings=1,
        forbidden_weight=LIMIT, total_columns=16, constituent_dimension=8,
        dense_nonempty_column_floor=9, holder_exponent=16, arities=list(range(3, 17)),
        fixed_dense_row_sets=65399, arities_kept_separate=True,
        dense_rowset_multiplicity_inside_holder=True, extra_dense_rowset_factor=False,
        fixed_rows_across_columns=True, fixed_source_geometry=True, pair_swaps=True,
        ordered_regional_matrix_products=True, ordered_column_matrix_polynomial=True,
        initial_state=0, terminal_vector=[1, 1], common_potential_used=False,
        endpoint_potential_factors=0, projective_normalizations=1, alpha=1,
        frozen_marks=[str(z) for z in MARKS], requested_banks=4, optimizer_used=False,
        mark_selection_per_arity=True, extra_mark_bank_factor=False,
        directed=False, whole_code_certificate=False, max_seconds=max_seconds, banks=[])
    try:
        _check_deadline(deadline)
        weights = {a: arity_weights(a) for a in range(3, 17)}
        for mark in MARKS:
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
            bank = dict(mark=str(mark), entries=rows,
                failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in rows])))
            result['banks'].append(bank)
            if progress:
                progress({key: value for key, value in bank.items() if key != 'entries'})
        _check_deadline(deadline)
        selected = []
        for index, arity in enumerate(range(3, 17)):
            choice = min(range(4), key=lambda j: result['banks'][j]['entries'][index]['failure_log2_screen'])
            bank = result['banks'][choice]
            selected.append(dict(bank['entries'][index], mark=bank['mark']))
        result.update(status='complete-r1-dense-four-mark-screen', complete_selected_cohort=True,
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
        for name in (Path(__file__).name, 'hd_two_sum_gf16_single_dense_matrix.py',
                     'hd_two_sum_gf16_single_dense.py', 'hd_two_sum_closed_kernel.py',
                     'hd_two_sum_gf16_arity_dense.py', 'hd_xor_pair_route.py')}
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'banks'}, allow_nan=False))

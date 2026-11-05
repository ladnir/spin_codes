"""Bounded FLOAT control for labelled middle-arity Hlambda columns.

This uses one reviewed common positive matrix/path measure throughout.
For each fixed original row set A, the exact fibre factors C,Z satisfy
K-1 <= (C-1)Z+(Z-1) at atoms in [0,1]. Minkowski and Holder therefore give
||K-1||_r <= ||C-1||_(2r)||Z||_(2r)+||Z-1||_r, entry by entry.
The empty branch retains its full sixteen-region matrix moment J_empty.
After entrywise group powers, the chronological activity polynomial keeps
degree >=9r. This drops individual RS row floors but excludes zero words.

Only prescribed fixed-A diagnostics are returned. No representative is
multiplied by a row-set count, and no all-arity certificate is claimed.
"""
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_lambda_fiber_outer import column_fiber_counts
from hd_pair_identity_kernel import mask_table_from_entries
from hd_pair_scheduled_control import load_measure
from hd_pair_state_kernel import common_moment_atoms
from hd_pair_state_tail_general import tail_from_polynomial
from hd_scheduled_column_control import column_word_polynomial

DIRECTORY = Path(__file__).resolve().parent
SOURCES = ((0,1,2,3,4), (0,1,2,4,8), tuple(range(8)), (0,1,2,3,4,5,8,12))


def moment_table(measure, exponent):
    """Moments above r are only used inside the SAME positive measure."""
    if type(exponent) is not int or not 1 <= exponent <= 36:
        raise ValueError('A moment exponent in1..36 is required')
    atoms, matrices = common_moment_atoms(measure)
    if np.any(atoms > 0):
        raise ValueError('Fibre domination requires atoms at most one')
    return mask_table_from_entries(np.logaddexp2.reduce(matrices, axis=0),
        np.logaddexp2.reduce(matrices+exponent*atoms[:,None,None], axis=0))


def support_norm(counts, table, exponent, *, omit_zero=False):
    """Positive Minkowski sum; exact zero support is omitted, never subtracted."""
    terms = [log2(count)+table[mask]/exponent for mask,count in counts.items()
             if count and (mask or not omit_zero)]
    return np.logaddexp2.reduce(terms, axis=0) if terms else np.full((2,2),-np.inf)


def active_norm(data, table_r, table_2r, groups):
    outside, zero = data['outside_counts'], data['zero_fiber_counts']
    if outside.get(0) != 1 or zero.get(0) != 1:
        raise ValueError('Both exact zero coefficients must equal one')
    first = support_norm(outside,table_2r,2*groups,omit_zero=True)
    first += support_norm(zero,table_2r,2*groups)
    second = support_norm(zero,table_r,groups,omit_zero=True)
    return np.logaddexp2(first,second)


def fixed_source(source, table_r, table_2r, prefix, *, groups=18, check_time=None):
    columns, ranks, masses = [], [], []
    for j in range(16):
        if check_time: check_time()
        data = column_fiber_counts(source,j)
        upper = active_norm(data,table_r,table_2r,groups)
        columns.append(np.stack((table_r[0],groups*upper)))
        ranks.append(data['cross_rank'])
        masses.append([len(data['outside_counts']),len(data['zero_fiber_counts'])])
    polynomial = column_word_polynomial(columns,np.zeros(2),np.array([0,1]),groups,
                                       check_time=check_time)
    arity = len(source)
    tail = float(tail_from_polynomial(polynomial,9*groups))
    root = tail/groups-8*arity*log2(255)
    return dict(original_source=list(source),arity=arity,cross_ranks=ranks,
        positive_support_counts=masses,normalized_fixed_A_root_log2=root,
        fixed_A_repeated_r_bound_log2=prefix+groups*root,
        rowset_sum_included=False,minimum_nonempty_columns_per_word=9,
        relaxed_total_activity_floor=9*groups)


def run(*,sources=SOURCES,max_seconds=120,progress=None):
    start = perf_counter()
    def check_time():
        if perf_counter()-start > max_seconds:
            raise TimeoutError('Bounded middle-fibre control exceeded its budget')
    saved,measure,prefix = load_measure()
    table_r,table_2r = moment_table(measure,18),moment_table(measure,36)
    rows = []
    for source in sources:
        row = fixed_source(source,table_r,table_2r,prefix,check_time=check_time)
        rows.append(row)
        if progress: progress(dict(**row,elapsed_seconds=perf_counter()-start))
    names = ('hd_lambda_fiber_control.py','hd_lambda_fiber_outer.py',
        'hd_pair_scheduled_outer.py','hd_bit_syndrome_flat_geometry.py',
        'hd_pair_identity_kernel.py','hd_pair_state_kernel.py',
        'hd_post_refresh_quiet_rows.py','hd_pair_scheduled_control.py',
        'hd_pair_scalar_control.py','hd_scheduled_column_control.py',
        'hd_pair_state_tail_general.py','hd_pair_state_tail.py')
    return dict(construction='scheduled-Hlambda-bc-no-cross',schedule='input-inverse-output-xor',
        groups=18,source_class='Fmax<=1',mark=saved['mark'],source_count=saved['source_count'],
        common_measure=measure,common_prefix_log2=prefix,full_independent_D1=True,
        empty_branch_is_full_matrix_moment=True,coefficient_subtraction=False,
        representative_extrapolation=False,original_row_floors_dropped=True,
        whole_code_certificate=False,directed=False,controls=rows,
        elapsed_seconds=perf_counter()-start,
        source_pins={name:sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in names})


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite evidence')
    result=run(progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')

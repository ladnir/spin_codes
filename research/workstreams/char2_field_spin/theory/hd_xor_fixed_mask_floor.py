"""Exact fixed-A comparison masks for one dense XOR2 group, FLOAT only.

For each fixed active-row set A, uniform input on GF(q)^A includes zero
coordinates. Its exact output-mask counts N_A(T) sum to q^|A|, and its empty
mask count is one. Independent output diagonals make amplitudes uniform
nonzero conditional on the mask. Hence Q_A=sum_T N_A(T)/q^|A| Q_T uses
the forced-nonzero regional matrices without an additional amplitude factor.

This comparison law is NOT the law of actual nonzero RS rows. Their expected
measure is pointwise dominated by beta_a=q^(32a)/(q-1)^(16a) times this
column-product law. We keep A fixed across all32 columns. Invertibility of
D2 H D1 preserves empty columns, so every actual nonzero group has at least
17 nonempty columns. Restrict the comparison to that event using a positive
matrix polynomial before summing its marked mass. No independent redraw of
A, row-set multiplicity, product-density c_a(p), or endpoint potential occurs.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_column_probe import uniform_column_support_counts
from hd_xor_dense_column_floor import (R1_MARKS, column_polynomial, log_identity,
    matrix_product, tail_log_mass)
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT


def beta_log2(active, n=32, k=16, q=256):
    if not 1 <= active <= 8 or not 0 < k <= n or q < 2:
        raise ValueError('Invalid active count, code dimensions, or field order')
    return n*active*log2(q)-(n-k)*active*log2(q-1)


def mask_operators(regions):
    """Chronological products for every labelled output support mask."""
    regions = np.asarray(regions, dtype=float)
    if (regions.ndim != 4 or regions.shape[1:] != (2, 2, 2)
            or not 2 <= len(regions) <= 8 or len(regions) % 2
            or np.any(np.isnan(regions)) or np.any(np.isposinf(regions))):
        raise ValueError('Expected two occupancy matrices per labelled region of an even group')
    masks = np.arange(1 << len(regions))
    result = np.tile(log_identity(), (len(masks), 1, 1))
    for position, table in enumerate(regions):
        result = matrix_product(result, table[(masks >> position) & 1])
    return result


@lru_cache(None)
def comparison_counts(active_mask, group=8, q=256):
    if type(active_mask) is not int or not 0 < active_mask < 1 << group:
        raise ValueError('A nonempty fixed active-row mask is required')
    counts = tuple(uniform_column_support_counts(q, group, active_mask))
    if counts[0] != 1 or sum(counts) != q**active_mask.bit_count():
        raise ArithmeticError('Uniform-input mask normalization failed')
    return counts


def fixed_column_split(operators, active_mask, group=8, q=256):
    operators = np.asarray(operators, dtype=float)
    if (operators.shape != (1 << group, 2, 2)
            or np.any(np.isnan(operators)) or np.any(np.isposinf(operators))):
        raise ValueError('A complete labelled mask-operator table is required')
    counts = comparison_counts(active_mask, group, q)
    normalizer = active_mask.bit_count()*log2(q)
    empty = operators[0]-normalizer
    nonempty = np.full((2, 2), -np.inf)
    for mask, count in enumerate(counts[1:], 1):
        if count:
            nonempty = np.logaddexp2(nonempty, operators[mask]+log2(count)-normalizer)
    return empty, nonempty


def evaluate(mark):
    started = perf_counter()
    z = F(mark)
    if not 0 < z < 1:
        raise ValueError('A rational weight mark in(0,1) is required')
    raw = pair_regions_screen(float(z), 1)
    operators = mask_operators(np.array([[raw[0, 0], raw[1, 0]]]*8))
    prefix = log2(REGION_LENGTH)-log2(Q-1)-LIMIT*log2(float(z))
    rows, arities = [], []
    for active_mask in range(1, 256):
        active = active_mask.bit_count()
        if active < 3:
            continue
        empty, nonempty = fixed_column_split(operators, active_mask)
        polynomial = column_polynomial(empty, nonempty)
        by_degree = np.logaddexp2(polynomial[:, 0, 0], polynomial[:, 0, 1])
        beta = beta_log2(active)
        rows.append(dict(active_row_mask=active_mask, active_constituents=active,
            beta_log2=beta, column_probability_normalizer=256**active,
            unrestricted_log2_screen=prefix+beta+tail_log_mass(polynomial, 0),
            failure_log2_screen=prefix+beta+tail_log_mass(polynomial, 17),
            largest_retained_degree=int(17+np.argmax(by_degree[17:])),
            nonempty_column_degree_log2_masses=by_degree.tolist()))
    for active in range(3, 9):
        selected = [row for row in rows if row['active_constituents'] == active]
        largest = max(selected, key=lambda row: row['failure_log2_screen'])
        arities.append(dict(active_constituents=active, fixed_row_sets=len(selected),
            failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in selected])),
            unrestricted_log2_screen=float(np.logaddexp2.reduce([row['unrestricted_log2_screen'] for row in selected])),
            largest_row_mask=largest['active_row_mask'], largest_retained_degree=largest['largest_retained_degree']))
    if len(rows) != 219:
        raise ArithmeticError('Incomplete active-row-set coverage')
    return dict(mark=str(z), fixed_row_sets=len(rows), entries=rows, arity_summaries=arities,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in rows])),
        unrestricted_log2_screen=float(np.logaddexp2.reduce([row['unrestricted_log2_screen'] for row in rows])),
        elapsed_seconds=perf_counter()-started)


def run(progress=None):
    witnesses = []
    for source, mark in R1_MARKS:
        row = dict(witness_source=source, **evaluate(mark))
        witnesses.append(row)
        if progress:
            progress({key: value for key, value in row.items() if key != 'entries'})
    return dict(construction=CONSTRUCTION, active_groups=1,
        active_constituent_counts=[3, 4, 5, 6, 7, 8], fixed_row_sets=219,
        exact_uniform_input_mask_mixture=True, uniform_input_includes_zero=True,
        conditional_nonzero_amplitude_operator=True, active_row_set_fixed_across_columns=True,
        comparison_is_actual_rs_law=False, beta_formula='q^(32a)/(q-1)^(16a)',
        source_supports=8192, source_matching_counts=[0, 0], field_order=256,
        forbidden_weight=LIMIT, nonempty_column_floor=17, total_columns=32,
        floor_unit='eight-coordinate-outer-columns', alpha=1,
        endpoint_potential_used=False, product_density_comparison_used=False,
        matrix_state_preserved_across_columns=True, positive_empty_nonempty_split=True,
        projective_normalizations=1, extra_row_set_binomial_factor=False,
        directed=False, proposal_requires_fresh_directed_replay=True,
        whole_code_certificate=False, fitted_new_witnesses=0,
        witnesses=witnesses, best_log2_screen=min(row['failure_log2_screen'] for row in witnesses))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a previous diagnostic artifact')
    result = run(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(output=str(args.output.resolve()), best_log2_screen=result['best_log2_screen'])), flush=True)

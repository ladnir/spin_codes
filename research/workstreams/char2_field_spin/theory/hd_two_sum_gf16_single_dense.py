"""Bounded FLOAT r1 dense screen, keeping all14 arities separate.

One active group has one fixed dense row set. In every region its output
support is either empty or one location, so f=0 under both matchings.
The existing weight DP with zero sparse groups is exact for the supplied
regional envelopes. Per-arity Holder weights retain the fixed row-set
count across16 columns. Three frozen marks are tested, with deterministic
selection per arity; no new witness is fitted and no bank factor is paid.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_gf16_arity_dense import weight_vectors, arity_weights
from hd_two_sum_gf16_full_dense import _check_deadline
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS
from hd_xor_dense_holder_classes import mask_branch_caps
from hd_xor_mixed_mask_floor import bivariate_log_tail
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


WITNESSES = (
    ('sparse-inherited', F('996248217/1000000000'), F('-127/2')),
    ('dense-inherited', F('8970988817957597/9007199254740992'), F('-2266954815/33554432')),
    ('h2-inherited', F('1116892707587883/1125899906842624'), F(-40)),
)


def single_source_tables(raw):
    """No source-geometry maximum is needed for occupancies0/1."""
    tables = np.array([raw[0, 0], raw[1, 0]], dtype=float)
    if tables.shape != (2, 2, 2) or np.any(np.isnan(tables)) or np.any(np.isposinf(tables)):
        raise ValueError('Finite-or-zero occupancy0/1 log matrices required')
    return np.repeat(tables[None], 16, axis=0)


def arity_tail(vectors, logp, weights):
    """Scalar empty/nonempty column power, retaining at least9 of16 columns."""
    empty, nonempty = mask_branch_caps(vectors, logp, weights)
    if empty.shape != (1,) or nonempty.shape != (1,):
        raise ValueError('The single-dense screen has exactly one empty sparse category')
    return bivariate_log_tail(empty, nonempty, columns=16, x_floor=0, y_floor=9)


def run(*, max_seconds=20, progress=None):
    if (isinstance(max_seconds, bool) or not isinstance(max_seconds, (int, float))
            or not np.isfinite(max_seconds) or not 0 < max_seconds <= 20):
        raise ValueError('A positive compute budget at most20 seconds is required')
    started = perf_counter()
    deadline = started+max_seconds
    result = dict(construction=CONSTRUCTION, labels=list(LABELS), field_order=256,
        field_polynomial='0x11b', active_groups=1, cohort_counts=[0, 0, 1],
        all_source_supports=True, source_supports=8192, source_colorings=1,
        source_full_pairs=0, forbidden_weight=LIMIT, total_columns=16,
        constituent_dimension=8, sparse_groups=0, sparse_support_floor=0,
        sparse_normalizer_exponent=0, dense_nonempty_column_floor=9,
        holder_exponent=16, arities=list(range(3, 17)), fixed_dense_row_sets=65399,
        arities_kept_separate=True, dense_rowset_multiplicity_inside_holder=True,
        extra_dense_rowset_factor=False, fixed_rows_across_columns=True,
        fixed_source_geometry=True, pair_swaps=True,
        potential_quotient_after_dense_mask_mixture=True, alpha=1,
        projective_normalizations=1, endpoint_potential_factors=1,
        witness_banks_requested=3, fitted_new_witnesses=0,
        witness_selection_per_arity=True, extra_witness_bank_factor=False,
        directed=False, whole_code_certificate=False, max_seconds=max_seconds, banks=[])
    try:
        _check_deadline(deadline)
        weights = {a: arity_weights(a) for a in range(3, 17)}
        for name, z, logp in WITNESSES:
            _check_deadline(deadline)
            raw = pair_regions_screen(float(z), 1)
            tables = single_source_tables(raw)
            vectors = weight_vectors(tables, float(logp), sparse_groups=0, deadline=deadline)
            prefix = log2(8192)-log2(255)-LIMIT*log2(float(z))-float(logp)
            rows = []
            for arity in range(3, 17):
                tail, degrees = arity_tail(vectors, float(logp), weights[arity][0])
                peak = max(degrees, key=lambda row: row['log2_mass'])
                rows.append(dict(active_constituents=arity, fixed_row_sets=comb(16, arity),
                    exact_moments_sha256=weights[arity][1], failure_log2_screen=prefix+tail,
                    largest_dense_degree=peak['nonempty_dense_columns']))
            bank = dict(witness_label=name, mark=str(z), log2_potential=str(logp), entries=rows,
                failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in rows])))
            result['banks'].append(bank)
            if progress:
                progress({key: value for key, value in bank.items() if key != 'entries'})
        _check_deadline(deadline)
        selected = []
        for index, arity in enumerate(range(3, 17)):
            choice = min(range(3), key=lambda j: result['banks'][j]['entries'][index]['failure_log2_screen'])
            bank = result['banks'][choice]
            selected.append(dict(bank['entries'][index], witness_label=bank['witness_label'],
                                 mark=bank['mark'], log2_potential=bank['log2_potential']))
        result.update(status='complete-r1-dense-arity-screen', complete_selected_cohort=True,
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
        for name in (Path(__file__).name, 'hd_two_sum_closed_kernel.py', 'hd_two_sum_gf16_arity_dense.py',
                     'hd_xor_dense_holder_classes.py', 'hd_xor_mixed_mask_floor.py', 'hd_xor_pair_route.py')}
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'banks'}, allow_nan=False))

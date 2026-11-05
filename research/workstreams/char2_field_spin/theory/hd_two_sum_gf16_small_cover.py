"""Bounded FLOAT all-source G16 screen for r1..17 and zero/one dense group.

One entrywise regional envelope covers EVERY fixed source support. The
source count is C(8192,r), not a selected matching class. Dense row sets
of every arity3..16 are combined inside one exact Holder moment. Neither
the source nor the dense row set is resampled between outer columns.
This is a new conservative comparison, not a transferred certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_affine_outer import column_action_float
from hd_two_sum_closed_kernel import arity_holder_moments
from hd_two_sum_gf16_arity_dense import weight_vectors
from hd_two_sum_gf16_full_dense import Z, LOGP
from hd_two_sum_gf16_sparse import CONSTRUCTION, MARK, LOG_POTENTIAL
from hd_xor_dense_holder_classes import class_log_weights, mask_branch_caps
from hd_xor_joint_columns import categories
from hd_xor_mixed_gate import mixed_column_log_coefficients, sparse_log_count
from hd_xor_mixed_mask_floor import bivariate_log_tail
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


WITNESSES = (
    ('sparse-inherited', MARK, LOG_POTENTIAL),
    ('dense-inherited', Z, LOGP),
    ('h2-inherited', F('1116892707587883/1125899906842624'), F(-40)),
)


def dense_moments():
    """Sum all fixed-A moments with one common denominator, exactly."""
    denominator = 255**128
    result = [0]*17
    for arity in range(3, 17):
        nums, den = arity_holder_moments(arity, 16, 16, 8)
        factor, remainder = divmod(denominator, den)
        if remainder:
            raise ArithmeticError('Incompatible exact moment denominators')
        for weight, value in enumerate(nums):
            result[weight] += factor*value
    return tuple(result), denominator


def source_envelope(raw, maximum=17):
    """Elementwise bound for every regional subset of every source support."""
    if type(maximum) is not int or not 1 <= maximum <= 17:
        raise ValueError('Regional maximum1..17 required')
    rows = [np.maximum.reduce([raw[m, full] for full in range(m//2+1)])
            for m in range(maximum+1)]
    return np.repeat(np.asarray(rows)[None], 16, axis=0)


def support_polynomial(caps, n1, n2):
    """Fixed sparse row support choices, including255 per active symbol."""
    poly = mixed_column_log_coefficients(caps, n1, n2)
    return poly+np.arange(len(poly))*log2(255)


def cohort_prefix(r, n1, n2, dense, z, logp):
    if (any(type(x) is not int or x < 0 for x in (r, n1, n2, dense))
            or not 1 <= r <= 17 or dense not in (0, 1) or n1+n2+dense != r):
        raise ValueError('A nonempty bounded zero/one-dense cohort is required')
    colors = comb(r, dense)*comb(r-dense, n1)
    return (log2(comb(8192, r))+log2(colors)+4*n1+n2*log2(120)
            -log2(255)-LIMIT*log2(float(z))-float(logp))


def one_witness(maximum, z, logp, *, deadline=None, progress=None):
    raw = pair_regions_screen(float(z), maximum)
    tables = source_envelope(raw, maximum)
    sparse_caps = {kind: float(np.max(column_action_float(tables, float(logp), *kind)
                                      -[0., float(logp)]))
                   for kind in categories(maximum)}
    vectors = weight_vectors(tables, float(logp), maximum-1, deadline=deadline)
    weights = class_log_weights(*dense_moments(), columns=16)
    branches = mask_branch_caps(vectors, float(logp), weights)
    dense_caps = [dict(zip(categories(maximum-1), branch)) for branch in branches]
    entries = []
    for r in range(1, maximum+1):
        if deadline is not None and perf_counter() > deadline:
            raise TimeoutError('Bounded small-coverage screen exhausted its budget')
        current = []
        for dense in (0, 1):
            for n1 in range(r-dense+1):
                n2 = r-dense-n1
                prefix = cohort_prefix(r, n1, n2, dense, z, logp)
                if dense:
                    tail, _ = bivariate_log_tail(
                        *(support_polynomial(caps, n1, n2) for caps in dense_caps),
                        columns=16, x_floor=9*(n1+2*n2), y_floor=9)
                    tail -= 8*(n1+2*n2)*log2(255)
                else:
                    tail = sparse_log_count(sparse_caps, n1, n2, n=16, k=8)
                row = dict(active_groups=r, cohort_counts=[n1, n2, dense],
                    source_supports=comb(8192, r), source_colorings=comb(r, dense)*comb(r-dense, n1),
                    sparse_support_floor=9*(n1+2*n2), dense_nonempty_column_floor=9*dense,
                    failure_log2_screen=prefix+tail)
                current.append(row)
        entries.extend(current)
        if progress:
            progress(dict(active_groups=r, cohorts=len(current),
                union_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in current]))))
    return entries


def run(maximum=17, *, max_seconds=300, progress=None):
    if type(maximum) is not int or not 1 <= maximum <= 17 or not 0 < max_seconds <= 600:
        raise ValueError('Maximum1..17 and positive bounded compute time required')
    started = perf_counter()
    deadline = started+max_seconds
    banks = []
    for name, z, logp in WITNESSES:
        rows = one_witness(maximum, z, logp, deadline=deadline,
            progress=(lambda row: progress(dict(witness=name, **row))) if progress else None)
        banks.append(dict(name=name, mark=str(z), log2_potential=str(logp), entries=rows))
    chosen = []
    for candidates in zip(*(bank['entries'] for bank in banks)):
        index = min(range(len(banks)), key=lambda i: candidates[i]['failure_log2_screen'])
        chosen.append(dict(candidates[index], witness=banks[index]['name']))
    worst = max(chosen, key=lambda row: row['failure_log2_screen'])
    return dict(construction=CONSTRUCTION, active_group_range=[1, maximum],
        covered_dense_counts=[0, 1], all_source_supports=True,
        source_geometry_envelope=True, fixed_rows_across_columns=True,
        all_dense_rowsets_inside_holder=True, fixed_dense_rowsets=65399,
        extra_dense_arity_factor=False, pair_swaps=True, directed=False,
        whole_code_certificate=False, fitted_new_witnesses=0,
        selected_cells=len(chosen), source_banks=len(banks),
        failure_log2_screen=float(np.logaddexp2.reduce([r['failure_log2_screen'] for r in chosen])),
        worst_cell=worst, entries=chosen, witnesses=banks, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum', type=int, default=17)
    parser.add_argument('--max-seconds', type=float, default=300)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a research artifact')
    result = run(args.maximum, max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, 'hd_two_sum_closed_kernel.py', 'hd_two_sum_gf16_arity_dense.py',
                     'hd_xor_pair_route.py', 'hd_xor_mixed_gate.py', 'hd_xor_mixed_mask_floor.py')}
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('entries', 'witnesses')}, allow_nan=False))

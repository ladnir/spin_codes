"""Selected r18/Fmax<=1 screen with independent single-column factors.

There is no added outer crossed layer. The global paired route and regional
measure are unchanged. Coefficientwise Holder is applied before composing
columns, not after grouping adjacent independent columns. The original RS
row set stays fixed in all sixteen columns. FLOAT only, not a certificate.
"""
from hashlib import sha256
from itertools import combinations
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_pair_scheduled_control import load_measure
from hd_pair_scheduled_outer import column_permutations, permute_mask
from hd_pair_identity_kernel import mask_log_matrices
from hd_pair_state_tail_general import (branch_polynomial, entrywise_power,
    matrix_polynomial_product, tail_from_polynomial)

DIRECTORY = Path(__file__).resolve().parent
CASES = (('lambda','input-inverse'), ('Hbin','input-inverse'), ('lambda','static'))


def column_branches(source, schedule, column, mixer, *, lambda_value=0xbc):
    if mixer == 'Hbin':
        from hd_pair_scheduled_outer import column_branches as binary
        return binary(source, schedule, column)
    if mixer != 'lambda':
        raise ValueError('Known mixer required')
    from hd_lambda_outer import column_branches as lam
    before, after = column_permutations(schedule, column)
    return tuple((w, permute_mask(mask, after), count) for w, mask, count in
        lam(tuple(before[i] for i in source), lambda_value=lambda_value))


def column_word_polynomial(columns, log_counts, input_weights, groups, *, check_time=None):
    """Chronological independent column moments, up to sixteen factors."""
    if not hasattr(columns, '__len__') or not 1 <= len(columns) <= 16:
        raise ValueError('One to sixteen columns required')
    counts = np.asarray(log_counts, dtype=float)
    word = np.full(counts.shape[:-1]+(1,2,2), -np.inf)
    word[...,0,0,0] = word[...,0,1,1] = 0.
    for column in columns:
        if check_time is not None:
            check_time()
        polynomial = branch_polynomial(column, counts, input_weights, groups)
        factor = entrywise_power(polynomial, groups, check_time=check_time)
        word = matrix_polynomial_product(word, factor, check_time=check_time)
    return word


def data_for_sources(sources, mixer, schedule, table, *, lambda_value=0xbc):
    columns, reference = [], None
    for j in range(16):
        rows = [column_branches(source, schedule, j, mixer,
                               lambda_value=lambda_value) for source in sources]
        counts = np.array([[log2(n) for _,_,n in row] for row in rows])
        weights = np.array([[w for w,_,_ in row] for row in rows], dtype=np.uint8)
        masks = np.array([[mask for _,mask,_ in row] for row in rows], dtype=np.uint32)
        if reference is None:
            reference = counts, weights
        else:
            np.testing.assert_array_equal(counts, reference[0])
            np.testing.assert_array_equal(weights, reference[1])
        columns.append(table[masks])
    return columns, *reference


def run(*, cases=CASES, max_seconds=120, progress=None):
    started = perf_counter()
    def check_time():
        if perf_counter()-started > max_seconds:
            raise TimeoutError('Selected single-column screen exceeded its budget')
    saved, measure, prefix = load_measure()
    table = mask_log_matrices(measure, 18)
    results = []
    for mixer, schedule in cases:
        arities = []
        for arity in (1,2):
            sources = tuple(combinations(range(16), arity))
            columns, counts, weights = data_for_sources(sources, mixer, schedule, table)
            polynomial = column_word_polynomial(columns, counts, weights, 18,
                                                check_time=check_time)
            roots = tail_from_polynomial(polynomial, 9*arity*18)/18-8*arity*log2(255)
            result = dict(arity=arity, labelled_rowsets=len(sources),
                bound_log2=prefix+18*float(np.logaddexp2.reduce(roots)),
                largest_fixed_A=list(sources[int(np.argmax(roots))]),
                sources=[list(source) for source in sources],
                normalized_per_A_roots_log2=roots.tolist())
            arities.append(result)
            if progress:
                progress(dict(mixer=mixer, schedule=schedule, arity=arity,
                    bound_log2=result['bound_log2'], elapsed_seconds=perf_counter()-started))
        mixed = prefix+18*float(np.logaddexp2.reduce(
            [(row['bound_log2']-prefix)/18 for row in arities]))
        results.append(dict(mixer=mixer,schedule=schedule,arities=arities,
                            mixed_arities_1_2_log2=mixed))
    names = ('hd_scheduled_column_control.py','hd_pair_scheduled_control.py',
        'hd_pair_scheduled_outer.py','hd_lambda_outer.py','hd_pair_identity_kernel.py',
        'hd_pair_state_kernel.py','hd_pair_state_tail_general.py','hd_pair_state_tail.py',
        'hd_pair_scalar_control.py','hd_binary_syndrome_outer.py')
    return dict(outer_pair_layer=False, global_paired_route_unchanged=True,
        lambda_value=0xbc, groups=18, source_class='Fmax<=1', mark=saved['mark'],
        source_count=saved['source_count'], common_measure=measure,
        common_prefix_log2=prefix, same_original_A_across_all_columns=True,
        full_independent_D1=True, column_factors=16, aggregate_input_floor='9*a*18',
        directed=False, whole_code_certificate=False, native_implementation_changed=False,
        controls=results, elapsed_seconds=perf_counter()-started,
        source_pins={name:sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in names})


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = run(progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(dict(elapsed_seconds=result['elapsed_seconds'],
        mixed=[(c['mixer'],c['schedule'],c['mixed_arities_1_2_log2']) for c in result['controls']])),flush=True)

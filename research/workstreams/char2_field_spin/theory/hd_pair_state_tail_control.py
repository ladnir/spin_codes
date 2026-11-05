"""Selected low-arity FLOAT control with an exact positive input-weight tail.

This replaces the finite Chernoff floor in hd_pair_state_control, whose
singleton result was dominated by the excluded zero word. The tensor
kernel, common measure, and source class remain unchanged. Neither this
screen nor a successful selected-class result is a whole-code certificate.
"""
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_pair_scalar_control import load_bank
from hd_pair_state_control import branch_data
from hd_pair_state_kernel import (EXACT_RATIO, SYMMETRIC_MAJORANT,
                                  leaves_from_measure)
from hd_pair_state_kernel import branch_log_matrices
from hd_pair_state_tail import word_tail

DIRECTORY = Path(__file__).resolve().parent


def run(*, mixer='H0', max_seconds=120, progress=None):
    started = perf_counter()

    def check_time():
        if perf_counter() - started > max_seconds:
            raise TimeoutError('Exact-tail selected control exceeded its budget')

    saved, curves, prefix = load_bank()
    # Keep numerical parameters only; the pinned bank's exact dyadic
    # strings and contact metadata would describe the pre-lift measure.
    measure = {key: value for key, value in saved['state_measure'].items()
               if key.startswith('log2_')}
    lift = 2. ** -36
    for key in ('log2_probability', 'log2_scale01', 'log2_scale10', 'log2_scale11'):
        measure[key] += lift
    degrees = np.arange(19)
    ceiling = np.stack((measure['log2_probability'] + degrees * measure['log2_atom00'],
        measure['log2_scale01'] + degrees * measure['log2_atom01'],
        measure['log2_scale10'] + degrees * measure['log2_atom10'],
        np.full(19, measure['log2_scale11'])))
    ceiling[0, 0] = 0.
    if measure['log2_probability'] >= 0 or np.any(ceiling < curves):
        raise ArithmeticError('Lifted FLOAT common measure is invalid')
    results = []
    for mode in (EXACT_RATIO, SYMMETRIC_MAJORANT):
        leaves = leaves_from_measure(measure, 18, mode=mode)
        arities = []
        for arity in (1, 2):
            data = branch_data(arity, mixer=mixer)
            blocks = []
            for delta in range(8, 16):
                check_time()
                matrices = branch_log_matrices(data['first'], data['second'], delta, leaves)
                blocks.append(matrices.reshape(len(data['sources']), data['branch_count'], 2, 2))
            tail = word_tail(blocks, data['log_counts'], data['input_weights'],
                             18, 9 * arity * 18, check_time=check_time)
            roots = tail / 18 - 8 * arity * log2(255)
            row = dict(arity=arity, labelled_rowsets=len(data['sources']),
                minimum_total_input_weight=9 * arity * 18,
                bound_log2=prefix + 18 * float(np.logaddexp2.reduce(roots)),
                largest_fixed_A=list(data['sources'][int(np.argmax(roots))]),
                largest_fixed_A_without_rowset_sum_log2=prefix + 18 * float(np.max(roots)))
            arities.append(row)
            if progress:
                progress(dict(mode=mode, mixer=mixer, **row,
                              elapsed_seconds=perf_counter() - started))
        mixed = prefix + 18 * float(np.logaddexp2.reduce(
            [(row['bound_log2'] - prefix) / 18 for row in arities]))
        results.append(dict(mode=mode, arities=arities, mixed_arities_1_2_log2=mixed))
    names = ('hd_pair_state_tail_control.py', 'hd_pair_state_tail.py',
             'hd_pair_state_control.py', 'hd_pair_state_kernel.py',
             'hd_pair_scalar_control.py', 'hd_parity_outer.py')
    return dict(mixer=mixer, groups=18, source_class='Fmax<=1', mark=saved['mark'],
        source_count=saved['source_count'], common_prefix_log2=prefix,
        common_measure=measure, float_measure_log2_lift=lift,
        exact_positive_total_input_weight_tail=True,
        separate_original_column_activity_floor_retained=False,
        same_fixed_A_across_all_eight_blocks=True, all_labelled_arities=[1, 2],
        full_unnormalized_D1=True, controls=results,
        directed=False, whole_code_certificate=False,
        regional_provider_runs=0, native_implementation_changed=False,
        source_pins={name: sha256((DIRECTORY / name).read_bytes()).hexdigest() for name in names},
        elapsed_seconds=perf_counter() - started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mixer', choices=('H0', 'Hbin'), default='H0')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = run(mixer=args.mixer, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(mixer=args.mixer, elapsed_seconds=result['elapsed_seconds'],
        mixed=[(row['mode'], row['mixed_arities_1_2_log2']) for row in result['controls']])), flush=True)

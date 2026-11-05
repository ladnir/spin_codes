"""Prescribed scheduled-outer comparisons, selected r18/Fmax<=1 only.

No stochastic schedule search, provider generation, native implementation,
or certificate promotion. All bounds are FLOAT under one positive measure.
For singleton rows, compare the earlier aggregate-weight tail against a
whole-word Minkowski bound that retains >=9 original active columns first.
Neither relaxation is assumed to dominate the other.
"""
from hashlib import sha256
from itertools import combinations, product
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_pair_scalar_control import load_bank
from hd_pair_scheduled_outer import (SCHEDULES, branch_data as binary_data,
                                     column_permutations, permute_mask)
from hd_pair_state_kernel import SYMMETRIC_MAJORANT, leaves_from_measure, branch_log_matrices
from hd_pair_state_tail import word_tail
from hd_pair_singleton_word import whole_word_roots

DIRECTORY = Path(__file__).resolve().parent


def load_measure():
    saved, curves, prefix = load_bank()
    measure = {key: value for key, value in saved['state_measure'].items() if key.startswith('log2_')}
    for key in ('log2_probability','log2_scale01','log2_scale10','log2_scale11'):
        measure[key] += 2.**-36
    degrees = np.arange(19)
    ceiling = np.stack((measure['log2_probability']+degrees*measure['log2_atom00'],
        measure['log2_scale01']+degrees*measure['log2_atom01'],
        measure['log2_scale10']+degrees*measure['log2_atom10'],
        np.full(19,measure['log2_scale11'])))
    ceiling[0,0] = 0.
    if measure['log2_probability'] >= 0 or np.any(ceiling < curves):
        raise ArithmeticError('Common FLOAT measure fails the regional comparison')
    return saved, measure, prefix


def prepare_data(arity, schedule, mixer):
    if mixer == 'Hbin':
        return binary_data(arity, schedule)
    if mixer != 'lambda':
        raise ValueError('Unknown mixer')
    from hd_lambda_outer import column_branches
    sources = tuple(combinations(range(16),arity))
    first, second, logs, weights = [], [], [], []
    branch_count = (2 if arity == 1 else 7)**2
    for block in range(8):
        block_first, block_second = [], []
        for source in sources:
            columns = []
            for j in (2*block,2*block+1):
                before, after = column_permutations(schedule,j)
                branches = column_branches(tuple(before[i] for i in source), lambda_value=2,q=256)
                columns.append([(w,permute_mask(mask,after),n) for w,mask,n in branches])
            row_logs, row_weights = [], []
            for left,right in product(*columns):
                block_first.append(left[1]); block_second.append(right[1])
                row_logs.append(log2(left[2])+log2(right[2]))
                row_weights.append(left[0]+right[0])
            if block == 0:
                logs.append(row_logs); weights.append(row_weights)
            elif row_logs != logs[len(block_first)//branch_count-1] or row_weights != weights[len(block_first)//branch_count-1]:
                raise AssertionError('Schedule changed branch category weights')
        first.append(block_first); second.append(block_second)
    return dict(sources=sources,branch_count=branch_count,
        first=np.asarray(first,dtype=np.uint32),second=np.asarray(second,dtype=np.uint32),
        log_counts=np.asarray(logs),input_weights=np.asarray(weights))


def screen(data, leaves, arity, prefix, check_time):
    blocks = []
    for b in range(8):
        check_time()
        matrix = branch_log_matrices(data['first'][b],data['second'][b],8+b,leaves)
        blocks.append(matrix.reshape(len(data['sources']),data['branch_count'],2,2))
    roots = {'aggregate_input_weight_tail':word_tail(blocks,data['log_counts'],data['input_weights'],
                                                  18,9*arity*18,check_time=check_time)/18}
    if arity == 1:
        roots['whole_word_activity_before_minkowski'] = whole_word_roots(blocks,18,check_time=check_time)
    def summarize(norms):
        values = norms-8*arity*log2(255)
        return dict(bound_log2=prefix+18*float(np.logaddexp2.reduce(values)),
            largest_fixed_A=list(data['sources'][int(np.argmax(values))]),
            largest_fixed_A_without_rowset_sum_log2=prefix+18*float(np.max(values)))
    combined = np.minimum.reduce(list(roots.values()))
    return dict(arity=arity,labelled_rowsets=len(data['sources']),
                methods={method:summarize(norms) for method,norms in roots.items()},
                best_per_A=summarize(combined))


def run(*,mixer='Hbin',schedules=SCHEDULES,max_seconds=120,progress=None):
    started = perf_counter()
    def check_time():
        if perf_counter()-started > max_seconds:
            raise TimeoutError('Prescribed scheduled screen exceeded its budget')
    saved,measure,prefix = load_measure()
    leaves = leaves_from_measure(measure,18,mode=SYMMETRIC_MAJORANT)
    results = []
    for schedule in schedules:
        arities = []
        for arity in (1,2):
            data = prepare_data(arity,schedule,mixer)
            result = screen(data,leaves,arity,prefix,check_time)
            arities.append(result)
            if progress:
                progress(dict(mixer=mixer,schedule=schedule,**result,elapsed_seconds=perf_counter()-started))
        mixed = prefix+18*float(np.logaddexp2.reduce(
            [(row['best_per_A']['bound_log2']-prefix)/18 for row in arities]))
        results.append(dict(schedule=schedule,arities=arities,mixed_arities_1_2_log2=mixed))
    names = ['hd_pair_scheduled_control.py','hd_pair_scheduled_outer.py','hd_pair_singleton_word.py',
             'hd_pair_state_tail.py','hd_pair_state_kernel.py','hd_pair_scalar_control.py',
             'hd_binary_syndrome_outer.py']
    if mixer == 'lambda':
        names.append('hd_lambda_outer.py')
    return dict(mixer=mixer,lambda_value=2 if mixer=='lambda' else None,
        pair_mode=SYMMETRIC_MAJORANT,groups=18,source_class='Fmax<=1',mark=saved['mark'],
        source_count=saved['source_count'],common_measure=measure,common_prefix_log2=prefix,
        schedules=list(schedules),same_original_A_across_all_columns=True,
        full_independent_D1=True,exact_positive_tail=True,
        minimum_between_valid_per_A_norm_bounds=True,
        directed=False,whole_code_certificate=False,native_implementation_changed=False,
        controls=results,elapsed_seconds=perf_counter()-started,
        source_pins={name:sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in names})


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mixer',choices=('Hbin','lambda'),default='Hbin')
    parser.add_argument('--schedules',nargs='+',choices=SCHEDULES,default=list(SCHEDULES))
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = run(mixer=args.mixer,schedules=tuple(args.schedules),
        progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(dict(mixer=args.mixer,elapsed_seconds=result['elapsed_seconds'],
        mixed=[(r['schedule'],r['mixed_arities_1_2_log2']) for r in result['controls']])),flush=True)

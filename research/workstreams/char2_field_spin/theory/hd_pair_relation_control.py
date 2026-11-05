"""Complete labelled arity3/4 selected-class FLOAT comparison in chunks.

The original A is fixed across all columns. Relation kernels provide
exact raw input counts, not resampled messages. Regional and source scope
is r18,Fmax<=1 at the pinned mark; no full-code certificate is claimed.
"""
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_pair_scheduled_control import load_measure
from hd_pair_scheduled_outer import SCHEDULES
from hd_pair_relation_data import iter_branch_data, source_sets
from hd_pair_state_kernel import SYMMETRIC_MAJORANT, leaves_from_measure, branch_log_matrices
from hd_pair_state_tail_general import word_tail

DIRECTORY = Path(__file__).resolve().parent


def deduplicated_pair_moments(first, second, delta, leaves):
    """Kernel depends only on two labelled masks, not their source labels."""
    keys = (first.astype(np.uint64) << np.uint64(16)) | second
    distinct, inverse = np.unique(keys, return_inverse=True)
    a = (distinct >> np.uint64(16)).astype(np.uint32)
    b = (distinct & np.uint64(65535)).astype(np.uint32)
    values = branch_log_matrices(a,b,delta,leaves)
    return values[inverse],len(distinct)


def run(*,arity=3,schedule='input-xor',pair_layer=True,max_seconds=300,progress=None):
    started = perf_counter()
    def check_time():
        if perf_counter()-started > max_seconds:
            raise TimeoutError('Complete selected arity comparison exceeded its budget')
    saved,measure,prefix = load_measure()
    if pair_layer:
        leaves = leaves_from_measure(measure,18,mode=SYMMETRIC_MAJORANT)
    else:
        from hd_pair_identity_kernel import mask_log_matrices
        table = mask_log_matrices(measure,18)
    observed,roots = [],[]
    distinct_pairs = 0
    chunks = 0
    for data in iter_branch_data(arity,schedule):
        check_time()
        blocks = []
        for b in range(8):
            check_time()
            if pair_layer:
                matrix,count = deduplicated_pair_moments(data['first'][b],data['second'][b],8+b,leaves)
                distinct_pairs += count
            else:
                from hd_pair_identity_kernel import branch_log_matrices as identity_matrices
                matrix = identity_matrices(data['first'][b],data['second'][b],table)
            blocks.append(matrix.reshape(len(data['sources']),data['branch_count'],2,2))
        tail = word_tail(blocks,data['log_counts'],data['input_weights'],18,9*arity*18,
                         check_time=check_time)
        roots.extend((tail/18-8*arity*log2(255)).tolist())
        observed.extend(data['sources'])
        chunks += 1
        if progress:
            progress(dict(arity=arity,schedule=schedule,pair_layer=pair_layer,
                completed_rowsets=len(observed),total_rowsets=len(source_sets(arity)),
                partial_sum_log2=prefix+18*float(np.logaddexp2.reduce(roots)),
                elapsed_seconds=perf_counter()-started))
    if tuple(observed) != source_sets(arity):
        raise AssertionError('Labelled source coverage or order changed')
    value = prefix+18*float(np.logaddexp2.reduce(roots))
    names = ['hd_pair_relation_control.py','hd_pair_relation_data.py','hd_binary_outer_relations.py',
        'hd_pair_scheduled_outer.py','hd_binary_syndrome_outer.py','hd_pair_state_tail_general.py',
        'hd_pair_state_tail.py','hd_pair_scheduled_control.py','hd_pair_scalar_control.py',
        'hd_pair_state_kernel.py']
    if not pair_layer:
        names.append('hd_pair_identity_kernel.py')
    return dict(mixer='Hbin',arity=arity,schedule=schedule,pair_layer=pair_layer,
        pair_mode=SYMMETRIC_MAJORANT if pair_layer else None,
        groups=18,source_class='Fmax<=1',mark=saved['mark'],source_count=saved['source_count'],
        common_measure=measure,common_prefix_log2=prefix,
        labelled_rowset_count=len(observed),all_labelled_rowsets=True,
        original_A_fixed_across_all_columns=True,full_independent_D1=True,
        minimum_total_original_input_weight=9*arity*18,
        exact_positive_tail=True,separate_column_activity_floor=False,
        bound_log2=value,largest_fixed_A=list(observed[int(np.argmax(roots))]),
        largest_fixed_A_without_rowset_sum_log2=prefix+18*float(np.max(roots)),
        rowsets=[list(row) for row in observed],normalized_roots_log2=roots,
        chunks=chunks,unique_pair_evaluations_across_chunks=distinct_pairs,
        elapsed_seconds=perf_counter()-started,
        directed=False,whole_code_certificate=False,native_implementation_changed=False,
        source_pins={name:sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in names})


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arity',type=int,choices=(3,4),default=3)
    parser.add_argument('--schedule',choices=SCHEDULES,default='input-xor')
    parser.add_argument('--no-pair-layer',action='store_true')
    parser.add_argument('--max-seconds',type=float,default=300)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result=run(arity=args.arity,schedule=args.schedule,pair_layer=not args.no_pair_layer,
        max_seconds=args.max_seconds,progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({key:result[key] for key in ('arity','schedule','pair_layer','bound_log2',
        'largest_fixed_A','elapsed_seconds')}),flush=True)

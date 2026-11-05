"""Complete labelled low-arity FLOAT gate for scheduled Hlambda, no cross.

The GF4 coefficient subfield is explicitly embedded in AES-GF256. Counts
are raw category masses. Optional positive integer scores retain weighted
consequences of each original RS row's weight floor; their sum is at most4
so the positive polynomial degree remains at most1152. A score stays fixed
across all columns and group replicas. Different score bounds may be
minimized per original A, before summing roots.
"""
from hashlib import sha256
from itertools import combinations
import json
from math import log2
from pathlib import Path
from time import perf_counter
import numpy as np

from hd_lambda_scheduled_relations import category_data, column_masks
from hd_pair_scheduled_control import load_measure
from hd_pair_scheduled_outer import SCHEDULES
from hd_pair_identity_kernel import mask_log_matrices
from hd_scheduled_column_control import column_word_polynomial
from hd_pair_state_tail_general import tail_from_polynomial

DIRECTORY=Path(__file__).resolve().parent


def scored_weights(input_masks, score):
    score=tuple(score)
    if (not 1<=len(score)<=4 or any(type(c) is not int or c<1 for c in score)
            or sum(score)>4):
        raise ValueError('Positive integer row scores with sum at most4 required')
    if any(type(mask) is not int or not 0<=mask<1<<len(score) for mask in input_masks):
        raise ValueError('Original input activity mask out of range')
    return np.array([sum(c for i,c in enumerate(score) if mask>>i&1)
                     for mask in input_masks],dtype=np.uint8)


def run(*,arity=3,schedule='input-inverse-output-xor',scores=None,
        chunk_size=32,max_seconds=360,progress=None):
    if type(arity) is not int or not 1<=arity<=4:
        raise ValueError('Arity1..4 required')
    if schedule not in SCHEDULES or type(chunk_size) is not int or not 1<=chunk_size<=128:
        raise ValueError('Known schedule and bounded positive chunk size required')
    scores=(tuple([1]*arity),) if scores is None else tuple(tuple(s) for s in scores)
    if not scores or any(len(s)!=arity for s in scores):
        raise ValueError('Nonempty scores with one entry per original row required')
    categories=category_data(arity)
    degree_sets=[scored_weights(categories['category_input_masks'],s) for s in scores]
    raw_logs=np.array([log2(n) for n in categories['category_counts']])
    saved,measure,prefix=load_measure()
    started=perf_counter()
    def check_time():
        if perf_counter()-started>max_seconds:
            raise TimeoutError('Complete single-column gate exceeded its budget')
    table=mask_log_matrices(measure,18)
    sources=tuple(combinations(range(16),arity))
    roots=[];by_score=[[] for _ in scores];observed=[]
    for offset in range(0,len(sources),chunk_size):
        chunk=sources[offset:offset+chunk_size]
        columns=[]
        for j in range(16):
            check_time()
            masks=np.array([column_masks(a,schedule,j,lambda_value=0xbc) for a in chunk],
                           dtype=np.uint32)
            columns.append(table[masks])
        counts=np.broadcast_to(raw_logs,(len(chunk),len(raw_logs)))
        chunk_roots=[]
        for index,(score,degrees) in enumerate(zip(scores,degree_sets)):
            weights=np.broadcast_to(degrees,counts.shape)
            polynomial=column_word_polynomial(columns,counts,weights,18,check_time=check_time)
            values=tail_from_polynomial(polynomial,9*sum(score)*18)/18-8*arity*log2(255)
            by_score[index].extend(values.tolist());chunk_roots.append(values)
        roots.extend(np.minimum.reduce(chunk_roots).tolist());observed.extend(chunk)
        if progress:
            progress(dict(arity=arity,schedule=schedule,completed_rowsets=len(observed),
                total_rowsets=len(sources),partial_sum_log2=prefix+18*float(np.logaddexp2.reduce(roots)),
                elapsed_seconds=perf_counter()-started))
    if tuple(observed)!=sources:
        raise AssertionError('Complete original row-set coverage required')
    names=('hd_lambda_column_control.py','hd_lambda_scheduled_relations.py',
        'hd_quaternary_outer_relations.py','hd_lambda_outer.py','hd_pair_scheduled_outer.py',
        'hd_pair_scheduled_control.py','hd_pair_scalar_control.py','hd_pair_identity_kernel.py',
        'hd_pair_state_kernel.py','hd_scheduled_column_control.py','hd_pair_state_tail_general.py',
        'hd_pair_state_tail.py')
    return dict(mixer='Hlambda',lambda_aes_byte=0xbc,coefficient_field_order=4,
        arity=arity,schedule=schedule,outer_crossed_layer=False,global_paired_route_unchanged=True,
        groups=18,source_class='Fmax<=1',mark=saved['mark'],source_count=saved['source_count'],
        common_measure=measure,common_prefix_log2=prefix,
        labelled_rowset_count=len(sources),all_labelled_rowsets=True,
        original_A_fixed_across_all_columns=True,full_independent_D1=True,
        row_scores=[list(s) for s in scores],row_scores_fixed_across_columns_and_replicas=True,
        exact_positive_tail=True,bound_log2=prefix+18*float(np.logaddexp2.reduce(roots)),
        score_bounds_log2=[prefix+18*float(np.logaddexp2.reduce(v)) for v in by_score],
        largest_fixed_A=list(sources[int(np.argmax(roots))]),
        rowsets=[list(a) for a in sources],normalized_roots_log2=roots,
        normalized_roots_by_score_log2=by_score,elapsed_seconds=perf_counter()-started,
        directed=False,whole_code_certificate=False,native_implementation_changed=False,
        source_pins={name:sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in names})


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arity',type=int,choices=range(1,5),default=3)
    parser.add_argument('--schedule',choices=SCHEDULES,default='input-inverse-output-xor')
    parser.add_argument('--scores',nargs='+',help='Per-row positive integer scores, e.g. 11 21 12')
    parser.add_argument('--max-seconds',type=float,default=360)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Refuse to overwrite evidence')
    scores=None if args.scores is None else [tuple(map(int,s)) for s in args.scores]
    result=run(arity=args.arity,schedule=args.schedule,scores=scores,max_seconds=args.max_seconds,
        progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({key:result[key] for key in ('arity','schedule','bound_log2','largest_fixed_A',
                                                'elapsed_seconds')}),flush=True)

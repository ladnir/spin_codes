"""Selected dense-arity FLOAT gate using the exact dual column code.

For bare orthogonal M_j=Q_j Hlambda P_j, the dual of M_j F^A is M_j
F^(A-complement). The GF4 categories enumerate raw dual words, not primal
input weights. This controller DROPS all original RS weight/activity floors.
The full-D1 comparison and255^(-8a) factor remain. No directed certificate.

First apply the regional majorant to the original nonnegative primal
support polynomial. Only then use MacWilliams under that common positive
measure. Substituting moment bounds into an alternating dual expansion
would not be valid.
"""
from hashlib import sha256
from itertools import combinations
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_lambda_scheduled_relations import column_masks, category_data
from hd_pair_scheduled_control import load_measure
from hd_pair_scheduled_outer import SCHEDULES
from hd_post_refresh_quiet_rows import log_matmul
from hd_lambda_dual_kernel import mask_log_matrices, column_log_bound

DIRECTORY=Path(__file__).resolve().parent


def run(*,arity=13,schedule='input-inverse-output-xor',chunk_size=32,
        max_seconds=180,progress=None):
    if type(arity) is not int or not 12<=arity<=16:
        raise ValueError('Dense arity12..16 required')
    if schedule not in SCHEDULES or type(chunk_size) is not int or not 1<=chunk_size<=128:
        raise ValueError('Prescribed schedule and bounded chunk required')
    saved,measure,prefix=load_measure()
    started=perf_counter()
    def check_time():
        if perf_counter()-started>max_seconds:
            raise TimeoutError('Dense dual gate exceeded its budget')
    b=16-arity
    complements=tuple(combinations(range(16),b))
    raw_counts=(1,) if b==0 else category_data(b)['category_counts']
    table=mask_log_matrices(measure,18)
    roots=[];observed=[]
    for offset in range(0,len(complements),chunk_size):
        chunk=complements[offset:offset+chunk_size]
        word=np.full((len(chunk),2,2),-np.inf)
        word[:,0,0]=word[:,1,1]=0.
        for j in range(16):
            check_time()
            masks=(np.zeros((len(chunk),1),dtype=np.uint32) if b==0 else
                   np.array([column_masks(source,schedule,j,lambda_value=0xbc)
                             for source in chunk],dtype=np.uint32))
            factor=column_log_bound(masks,raw_counts,b,18,table)
            word=log_matmul(word,factor)
        norms=np.logaddexp2(word[:,0,0],word[:,0,1])/18-8*arity*log2(255)
        roots.extend(norms.tolist());observed.extend(chunk)
        if progress:
            progress(dict(arity=arity,completed_rowsets=len(observed),total_rowsets=len(complements),
                partial_sum_log2=prefix+18*float(np.logaddexp2.reduce(roots)),
                elapsed_seconds=perf_counter()-started))
    if tuple(observed)!=complements:raise AssertionError('Complement coverage changed')
    sources=[sorted(set(range(16))-set(row)) for row in complements]
    names=('hd_lambda_dual_control.py','hd_lambda_dual_kernel.py',
        'hd_lambda_scheduled_relations.py','hd_quaternary_outer_relations.py','hd_lambda_outer.py',
        'hd_pair_scheduled_outer.py','hd_pair_scheduled_control.py','hd_pair_scalar_control.py',
        'hd_pair_state_kernel.py','hd_pair_identity_kernel.py','hd_post_refresh_quiet_rows.py')
    return dict(mixer='Hlambda',lambda_aes_byte=0xbc,coefficient_field_order=4,
        arity=arity,dual_arity=b,schedule=schedule,outer_crossed_layer=False,
        global_paired_route_unchanged=True,groups=18,source_class='Fmax<=1',
        mark=saved['mark'],source_count=saved['source_count'],common_measure=measure,
        common_prefix_log2=prefix,full_independent_D1=True,original_RS_floors_retained=False,
        primal_positive_majorant_before_dual_rewrite=True,
        raw_dual_zero_word_included=True,all_labelled_rowsets=True,labelled_rowset_count=len(sources),
        original_complement_fixed_across_columns=True,dual_normalization_inside_group_power=True,
        bound_log2=prefix+18*float(np.logaddexp2.reduce(roots)),
        largest_fixed_A=sources[int(np.argmax(roots))],rowsets=sources,
        complements=[list(row) for row in complements],normalized_roots_log2=roots,
        elapsed_seconds=perf_counter()-started,directed=False,whole_code_certificate=False,
        native_implementation_changed=False,
        source_pins={name:sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in names})


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arity',type=int,choices=range(12,17),default=13)
    parser.add_argument('--schedule',choices=SCHEDULES,default='input-inverse-output-xor')
    parser.add_argument('--max-seconds',type=float,default=180)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Refuse to overwrite evidence')
    result=run(arity=args.arity,schedule=args.schedule,max_seconds=args.max_seconds,
        progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({key:result[key] for key in ('arity','bound_log2','largest_fixed_A',
                                                'elapsed_seconds')}),flush=True)

"""Prescribed middle-row diagnostics with positive RS deficit truncation.

Exact5 uses full GF4 relation categories, including determined inside output
coordinates and true original-input weights. Marked-fibre uses the reviewed
quiet-exact/nonquiet-inflated upper-tail comparison. Both retain the same
fixed original A, regional measure, matrix chronology and normalizers.
No row-set extrapolation or whole-code certificate is supplied.
"""
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter
import numpy as np

from hd_column_deficit_tail import column_deficit_tail
from hd_lambda_arity5_relations import joint_counts,category_data
from hd_lambda_fiber_control import moment_table
from hd_lambda_marked_fiber_control import (SOURCES,prepare_source,optimize_prepared,
                                          UnsupportedQuietFiber)
from hd_pair_scheduled_control import load_measure

DIRECTORY=Path(__file__).resolve().parent
EXACT_SOURCES=SOURCES[:2]


def exact_columns(source,table,groups=18,*,check_time=None):
    if len(source)!=5:
        raise ValueError('This exact extension is bounded to arity5')
    columns,counts=[],[]
    for j in range(16):
        if check_time:check_time()
        branches=joint_counts(source,'input-inverse-output-xor',j)
        if sum(branches.values())!=256**5:
            raise ArithmeticError('Raw GF256 input mass changed')
        column=np.full((6,2,2),-np.inf)
        for (weight,mask),count in branches.items():
            column[weight]=np.logaddexp2(column[weight],log2(count)+table[mask]/groups)
        columns.append(column);counts.append(len(branches))
    return columns,dict(exact_merged_branch_counts=counts,
        kernel_categories=len(category_data(5)['category_counts']),
        raw_column_mass=256**5,original_input_weights_exact=True)


def run(mode,*,max_seconds=120,progress=None):
    if mode not in ('exact5','marked-fibre') or not 0<max_seconds<=120:
        raise ValueError('Prescribed exact5/marked-fibre mode and at most120 seconds required')
    start=perf_counter()
    def check_time():
        if perf_counter()-start>max_seconds:
            raise TimeoutError('Bounded middle control exceeded its budget')
    saved,measure,prefix=load_measure()
    table,table2=moment_table(measure,18),moment_table(measure,36)
    sources=EXACT_SOURCES if mode=='exact5' else SOURCES
    rows=[];complete=True;error=None
    try:
        for source in sources:
            try:
                if mode=='exact5':
                    columns,diagnostics=exact_columns(source,table,check_time=check_time)
                else:
                    columns,diagnostics=prepare_source(source,table,table2,check_time=check_time)
                cher=optimize_prepared(np.asarray(columns),len(source),prefix,check_time=check_time)
                tail=float(column_deficit_tail(columns,len(source),18,check_time=check_time))
                root=tail/18-8*len(source)*log2(255)
                value=prefix+18*root
                if value>cher['best']['fixed_A_repeated_r_bound_log2']+1e-8:
                    raise ArithmeticError('Positive tail exceeded a Chernoff upper bound')
                row=dict(original_source=list(source),arity=len(source),status='complete',
                    fixed_A_repeated_r_bound_log2=value,normalized_fixed_A_root_log2=root,
                    chernoff=cher,diagnostics=diagnostics,retained_deficit=7*len(source)*18)
            except UnsupportedQuietFiber as exception:
                row=dict(original_source=list(source),arity=len(source),status='unsupported_zero_fiber',
                         error=str(exception))
            rows.append(row)
            if progress:
                progress({key:value for key,value in row.items() if key not in ('diagnostics','chernoff')})
    except TimeoutError as exception:
        complete=False;error=str(exception)
    names=('hd_lambda_middle_control.py','hd_column_deficit_tail.py',
        'hd_lambda_arity5_relations.py','hd_quaternary_outer_relations.py',
        'hd_lambda_marked_fiber_control.py','hd_lambda_fiber_control.py',
        'hd_lambda_fiber_outer.py','hd_lambda_outer.py','hd_outer_pair_layers.py',
        'hd_bit_syndrome_outer.py','hd_bit_syndrome_flat_geometry.py',
        'hd_pair_scheduled_outer.py','hd_pair_scheduled_control.py','hd_pair_scalar_control.py',
        'hd_inverse_column_control.py','hd_pair_identity_kernel.py','hd_pair_state_kernel.py',
        'hd_post_refresh_quiet_rows.py','hd_scheduled_column_control.py',
        'hd_pair_state_tail_general.py','hd_pair_state_tail.py')
    return dict(construction='scheduled-Hlambda-bc-no-cross',schedule='input-inverse-output-xor',
        mode=mode,groups=18,source_class='Fmax<=1',mark=saved['mark'],source_count=saved['source_count'],
        common_measure=measure,common_prefix_log2=prefix,full_independent_D1=True,
        original_A_fixed_across_columns=True,minimum_aggregate_input_activity='9*a*r',
        rowset_sum_included=False,representative_extrapolation=False,directed=False,
        whole_code_certificate=False,complete=complete,error=error,controls=rows,
        elapsed_seconds=perf_counter()-start,
        source_pins={name:sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in names})


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('exact5','marked-fibre'),required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Refuse to overwrite evidence')
    result=run(args.mode,progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(dict(complete=result['complete'],elapsed_seconds=result['elapsed_seconds'],
                         error=result['error'])),flush=True)

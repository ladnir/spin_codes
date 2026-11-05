"""Bounded complete low-arity matrix-state screen for H0 plus one pair layer.

The source scope is r18,F<=1 at one pinned output mark, not a whole-code
certificate. All labelled singleton and two-row sets are retained. Each
set persists across all eight column pairs. Positive branch Minkowski
precedes the chronological block product; whole-word Minkowski sums row
sets only after that product. Original column activity and input weight
are marked separately. Both pair constructions use full independent D1.
"""
from hashlib import sha256
from itertools import combinations, product
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_pair_scalar_control import load_bank
from hd_pair_state_kernel import (EXACT_RATIO, SYMMETRIC_MAJORANT, leaves_from_measure,
                                  branch_log_matrices, metadata)
from hd_parity_outer import support_counts
from hd_post_refresh_quiet_rows import log_matmul

DIRECTORY = Path(__file__).resolve().parent
MARKERS = ((1,1),(2,1),(4,1),(16,1),(256,1),(1,2),(1,4),(1,16),(2,2),(4,4))


def joint_branches(source, *, mixer='H0'):
    """Exact (input weight, output support, multiplicity), including zero."""
    source = tuple(source)
    if source != tuple(sorted(set(source))) or len(source) not in (1,2) or not all(0<=i<16 for i in source):
        raise ValueError('One or two distinct sorted row labels required')
    if mixer == 'H0':
        column = lambda i: 65535 ^ (1 << i)
    elif mixer == 'Hbin':
        column = lambda i: sum(((int(i==j)^1^((i&j).bit_count()&1)) & 1) << j for j in range(16))
    else:
        raise ValueError('Unknown binary first mixer')
    rows = [(0,0,1),(1,column(source[0]),255)]
    if len(source)==2:
        a,b=map(column,source)
        rows += [(1,b,255),(2,a^b,255),(2,a|b,255*254)]
    if mixer=='H0':
        merged={}
        for _,mask,count in rows:
            merged[mask]=merged.get(mask,0)+count
        if merged!=support_counts(sum(1<<i for i in source)):
            raise AssertionError('Independent parity enumerator disagrees')
    if sum(count for _,_,count in rows)!=256**len(source):
        raise AssertionError('Binary two-source mass changed')
    return tuple(rows)


def branch_data(arity, *, mixer='H0'):
    sources=tuple(combinations(range(16),arity))
    rows=[joint_branches(source,mixer=mixer) for source in sources]
    count=len(rows[0])**2
    first,second,weights,activities,input_weights=[],[],[],[],[]
    for branches in rows:
        for left,right in product(branches,repeat=2):
            first.append(left[1]);second.append(right[1])
            weights.append(log2(left[2])+log2(right[2]))
            activities.append(int(left[0]>0)+int(right[0]>0))
            input_weights.append(left[0]+right[0])
    return dict(sources=sources,branch_count=count,
        first=np.asarray(first,dtype=np.uint32),second=np.asarray(second,dtype=np.uint32),
        log_counts=np.asarray(weights).reshape(len(sources),count),
        activities=np.asarray(activities).reshape(len(sources),count),
        input_weights=np.asarray(input_weights).reshape(len(sources),count))


def combine_branches(log_matrices, log_counts, activities, input_weights, groups, x, y):
    """Entrywise Minkowski over one block's positive polynomial branches."""
    matrices=np.asarray(log_matrices,dtype=float)
    counts=np.asarray(log_counts,dtype=float)
    activities=np.asarray(activities);input_weights=np.asarray(input_weights)
    if (matrices.shape!=counts.shape+(2,2) or activities.shape!=counts.shape
            or input_weights.shape!=counts.shape or counts.ndim!=2
            or np.any(np.isnan(matrices)) or np.any(np.isposinf(matrices))
            or not np.all(np.isfinite(counts)) or np.any(activities<0) or np.any(input_weights<0)
            or type(groups) is not int or groups<1 or not 1<=x<1e100 or not 1<=y<1e100):
        raise ValueError('Compatible positive block data and markers>=1 required')
    terms=(counts+activities*log2(y)+input_weights*log2(x))[...,None,None]+matrices/groups
    return groups*np.logaddexp2.reduce(terms,axis=1)


def run(*, mixer='H0', max_seconds=60, progress=None):
    started=perf_counter()
    def check_time():
        if perf_counter()-started>max_seconds:
            raise TimeoutError('Bounded low-arity state control exceeded its budget')
    saved,curves,prefix=load_bank()
    measure=dict(saved['state_measure'])
    # Repair the recorded tiny FLOAT deficit conservatively before reuse.
    # This remains a FLOAT screen, not interval-verified evidence.
    lift=2.**-36
    for key in ('log2_probability','log2_scale01','log2_scale10','log2_scale11'):
        measure[key]+=lift
    if measure['log2_probability']>=0:
        raise ValueError('Invalid lifted nonzero-atom probability')
    degrees=np.arange(19)
    ceiling=np.stack((measure['log2_probability']+degrees*measure['log2_atom00'],
        measure['log2_scale01']+degrees*measure['log2_atom01'],
        measure['log2_scale10']+degrees*measure['log2_atom10'],
        np.full(19,measure['log2_scale11'])))
    ceiling[0,0]=0.
    if np.any(ceiling<curves):
        raise ArithmeticError('Lifted common measure misses a regional entry')
    results=[]
    for mode in (EXACT_RATIO,SYMMETRIC_MAJORANT):
        leaves=leaves_from_measure(measure,18,mode=mode)
        arities=[]
        for arity in (1,2):
            data=branch_data(arity,mixer=mixer)
            blocks=[]
            for delta in range(8,16):
                check_time()
                matrices=branch_log_matrices(data['first'],data['second'],delta,leaves)
                blocks.append(matrices.reshape(len(data['sources']),data['branch_count'],2,2))
            controls=[]
            for x,y in MARKERS:
                check_time()
                word=np.full((len(data['sources']),2,2),-np.inf)
                word[:,0,0]=word[:,1,1]=0.
                for block in blocks:
                    current=combine_branches(block,data['log_counts'],data['activities'],data['input_weights'],18,x,y)
                    word=log_matmul(word,current)
                row_words=np.logaddexp2.reduce(word[:,0,:],axis=1)
                fixed_offset=-8*arity*log2(255)-9*arity*log2(x)-9*log2(y)
                roots=row_words/18+fixed_offset
                value=prefix+18*float(np.logaddexp2.reduce(roots))
                controls.append(dict(input_weight_marker=x,column_activity_marker=y,bound_log2=value,
                    largest_fixed_A=list(data['sources'][int(np.argmax(roots))]),
                    largest_fixed_A_without_rowset_sum_log2=prefix+18*float(np.max(roots))))
            best=min(controls,key=lambda row:row['bound_log2'])
            arities.append(dict(arity=arity,labelled_rowsets=len(data['sources']),controls=controls,best=best))
            if progress:
                progress(dict(mode=mode,mixer=mixer,arity=arity,**best,elapsed_seconds=perf_counter()-started))
        mixed=prefix+18*float(np.logaddexp2.reduce([(row['best']['bound_log2']-prefix)/18 for row in arities]))
        results.append(dict(mode=mode,arities=arities,mixed_arities_1_2_log2=mixed,kernel_metadata=metadata(mode=mode)))
    pins={name:sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in
          ('hd_pair_state_control.py','hd_pair_state_kernel.py','hd_parity_outer.py','hd_pair_scalar_control.py')}
    return dict(mixer=mixer,groups=18,source_class='Fmax<=1',mark=saved['mark'],source_count=saved['source_count'],
        common_prefix_log2=prefix,common_measure=measure,float_measure_log2_lift=lift,
        original_activity_floor=9,original_input_weight_floor='9*arity',
        same_fixed_A_across_all_eight_blocks=True,all_labelled_rowsets_for_arities=[1,2],
        no_noncommutative_orbit_quotient=True,full_unnormalized_D1=True,
        controls=results,regional_provider_runs=0,native_implementation_changed=False,
        directed=False,whole_code_certificate=False,source_pins=pins,elapsed_seconds=perf_counter()-started)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--mixer',choices=('H0','Hbin'),default='H0')
    args=parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result=run(mixer=args.mixer,progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(dict(mixer=args.mixer,elapsed_seconds=result['elapsed_seconds'],
        mixed=[(row['mode'],row['mixed_arities_1_2_log2']) for row in result['controls']])),flush=True)

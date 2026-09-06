"""Engineering Q1 projections anchored only to complete exact spectra.

The model removes the minimum-shell multiplicity and row-position count,
fits the remaining margin against distance, then restores those costs at
the projected size. Window variation and largest-size holdouts expose its
sensitivity. These estimates are not full-distance bounds or certificates.
"""
import csv
import hashlib
import json
import math

import numpy as np

import read_grid_receipts
import run_occupation_grid as source


def rm_minimum(block):
    m=block.bit_length()-1
    if 1<<m!=block or m%2!=1 or m<3:raise ValueError('rate-half RM length required')
    r=(m-1)//2;numerator=denominator=1
    for j in range(r):numerator*=(1<<(m-j))-1;denominator*=(1<<(r-j))-1
    if numerator%denominator:raise ArithmeticError('nonintegral affine-subspace count')
    return 1<<(m-r),(1<<r)*(numerator//denominator)


def linear(x,y,target):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
    if len(x)<2 or np.ptp(x)==0:raise ValueError('two distinct anchors required')
    center=float(x.mean());denominator=float(np.square(x-center).sum())
    slope=float((x-center)@(y-y.mean()))/denominator
    return float(y.mean()+slope*(target-center))


def project(anchors, family, block, exponent, window):
    anchors=sorted(anchors,key=lambda r:r['block_bits'])[-window:]
    if len(anchors)!=window or window<2:raise ValueError('missing projection anchors')
    if family=='rm':
        distance,multiplicity=rm_minimum(block);log_count=math.log2(multiplicity)
    elif family=='bch':
        distance=2**linear([math.log2(r['block_bits']) for r in anchors],
                            [math.log2(r['distance']) for r in anchors],math.log2(block))
        log_count=linear([r['distance'] for r in anchors],[r['log_minimum_count'] for r in anchors],distance)
    else:raise ValueError('random codes are computed rather than extrapolated')
    normalized=[r['margin_bits']+r['log_minimum_count']+exponent-math.log2(r['block_bits']//2) for r in anchors]
    margin=linear([r['distance'] for r in anchors],normalized,distance)-log_count-exponent+math.log2(block//2)
    return dict(margin_bits=margin,distance_proxy=float(distance),log_minimum_count_proxy=float(log_count),
                anchor_blocks=[r['block_bits'] for r in anchors],window=window)


def write_csv(path,rows):
    with path.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    _,observations,_=read_grid_receipts.snapshot();spectra={};groups={};excluded=[]
    for row in observations.values():
        if row['outer_model']=='random-ensemble-average':continue
        family='rm' if row['series'].startswith('RM(') else 'bch'
        if row['series'] not in spectra:
            counts=source.exact_counts(row);d=min(counts)
            spectra[row['series']]=(d,math.log2(counts[d]),counts[d])
            if family=='rm' and rm_minimum(int(row['block_bits']))!=(d,counts[d]):
                raise ValueError('RM minimum-shell formula disagrees with exact spectrum')
        d,log_count,minimum_count=spectra[row['series']]
        if int(row['dominant_weight'])!=d or int(row['dominant_witness_at_grid_edge']):
            excluded.append(list(source.grid.key(row)));continue
        key=(family,int(row['step_bits']),int(row['state_bits']),int(row['message_exponent']))
        groups.setdefault(key,[]).append(dict(block_bits=int(row['block_bits']),distance=d,
            log_minimum_count=log_count,minimum_count=minimum_count,
            margin_bits=float(row['margin_bits']),map_tag=row['map_tag']))
    validation=[];projections=[];training=[]
    for (family,t,s,e),anchors in sorted(groups.items()):
        anchors.sort(key=lambda r:r['block_bits'])
        expected=[8,32,128,512] if family=='rm' else [8,32,64,128]
        if [r['block_bits'] for r in anchors]!=expected:continue
        training.extend(dict(family=family,step_bits=t,state_bits=s,message_exponent=e,**anchor) for anchor in anchors)
        actual=anchors[-1]
        for window in (2,3):
            prediction=project(anchors[:-1],family,actual['block_bits'],e,window)
            validation.append(dict(family=family,step_bits=t,state_bits=s,message_exponent=e,window=window,
                held_out_block=actual['block_bits'],observed_q1_margin_bits=actual['margin_bits'],
                projected_q1_margin_bits=prediction['margin_bits'],error_bits=prediction['margin_bits']-actual['margin_bits']))
        for block in ([2048] if family=='rm' else [256,512]):
            length=(1<<e)//(block//2)
            if length<t or length%t:continue
            candidates=[project(anchors,family,block,e,w) for w in (2,3,4)]
            margins=[c['margin_bits'] for c in candidates]
            projections.append(dict(family=family,block_bits=block,dimension=block//2,step_bits=t,state_bits=s,
                message_exponent=e,outer_rows=length,median_q1_margin_estimate=float(np.median(margins)),
                minimum_window_estimate=min(margins),maximum_window_estimate=max(margins),
                median_distance_proxy=float(np.median([c['distance_proxy'] for c in candidates])),
                median_log_minimum_count_proxy=float(np.median([c['log_minimum_count_proxy'] for c in candidates])),
                map_tag=anchors[0]['map_tag'],evidence='Q1 engineering estimate; unknown higher shells and occupations; no certificate'))
    summaries={}
    for family in ('bch','rm'):
        errors=np.array([r['error_bits'] for r in validation if r['family']==family])
        summaries[family]=dict(heldout_predictions=len(errors),median_error_bits=float(np.median(errors)),
            median_absolute_error_bits=float(np.median(abs(errors))),maximum_absolute_error_bits=float(max(abs(errors))),
            percentile95_absolute_error_bits=float(np.quantile(abs(errors),.95)))
    here=source.grid.HERE
    write_csv(here/'exact_family_q1_projections.csv',projections)
    write_csv(here/'exact_family_projection_holdouts.csv',validation)
    write_csv(here/'exact_family_projection_training.csv',training)
    payload=dict(status='ENGINEERING_Q1_ESTIMATES_ONLY',primary_exact_blocks=dict(bch=[8,32,64,128],rm=[8,32,128,512]),
        projection_rows=len(projections),excluded_boundary_or_nonminimum_rows=len(excluded),
        excluded_parameter_keys=excluded,holdout_summary=summaries,
        projection_csv_sha256=hashlib.sha256((here/'exact_family_q1_projections.csv').read_bytes()).hexdigest(),
        holdout_csv_sha256=hashlib.sha256((here/'exact_family_projection_holdouts.csv').read_bytes()).hexdigest(),
        training_csv_sha256=hashlib.sha256((here/'exact_family_projection_training.csv').read_bytes()).hexdigest(),
        producer_sha256=source.grid.pilot.sha(here/'extrapolate_exact_families.py'),
        limitations=['Only Q1 is projected; these estimates do not establish a full distance bound.',
            'BCH future distance and minimum multiplicity are fitted proxies, not properties of a specified code.',
            'RM minimum-shell parameters are exact, but its remaining spectrum is not supplied.',
            'Window spread and holdout errors describe model sensitivity; they are not confidence intervals.',
            'BCH-256 bounded-spectrum evidence is excluded from every fit; random spectra are not extrapolated.'])
    (here/'exact_family_projection_summary.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({k:v for k,v in payload.items() if k!='excluded_parameter_keys'},indent=2))


if __name__=='__main__':main()

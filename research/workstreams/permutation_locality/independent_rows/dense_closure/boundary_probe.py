"""Selected boundary compositions: diagnostic bounds, never a full cover."""
import argparse
from fractions import Fraction as Q
import json
from math import lgamma
from pathlib import Path
from time import monotonic

import numpy as np
from flint import arb,ctx

import kernel
import selected
import density_tangent
import shuffle_comparison
from mixture import actual_components,G,N,REGIONS,PACKETS,capped_density_loss,logq,log_power

CASES={
    'light65':{'0003':65},
    'light80':{'0003':80},
    'light128':{'0003':128},
    'central_one65':{'0002':65},
    'central_two65':{'0022':65},
    'central_four65':{'2222':65},
    'light_with_full':{'0003':64,'1111':1},
    'mixed65':{'0003':32,'0034':33},
    'full65':{'1111':65},
    'mixed_light65':{'0003':32,'0002':33},
    'near_light65':{'0003':64,'0002':1},
    'three_laws65':{'0003':30,'0002':30,'0034':5},
    'mixed_full65':{'0003':32,'0034':32,'1111':1},
    'overlap_full65':{'2222':64,'1111':1},
}
SCHEMA='four-bit-boundary-composition-probes-1'


def composition(components,case):
    lookup={row[0]:i for i,row in enumerate(components)}
    if any(name not in lookup or name=='0000' or type(n) is not int or n<0 for name,n in case.items()):
        raise ValueError('known active types and nonnegative integer counts required')
    if not 1<=sum(case.values())<=G:raise ValueError('active-group count in 1..2048 required')
    counts=[0]*len(components)
    for name,n in case.items():counts[lookup[name]]=n
    counts[lookup['0000']]=G-sum(counts)
    return counts


def parameters(result):
    lam,*logs=result['tilt']
    return [Q(round(lam*10**9),10**9),*(Q(round(float(np.exp(x))*10**12),10**12) for x in logs)]


def breakdown(data,components,counts,threshold,witness,anchors=None,exact_shuffle=False,posterior_shuffle=False):
    """Floating decomposition of one fixed witness, not proof arithmetic."""
    ids=[i for i,c in enumerate(counts) if c]
    n=np.array([counts[i] for i in ids])
    laws=np.array([list(map(float,components[i][2])) for i in ids])
    lam,*tilts=map(float,witness);z=laws@np.array([1.,*tilts])
    weights=(n/G/z)@laws
    caps=tuple(sum(counts[i] for i in ids if components[i][2][j]>0) for j in range(5))
    if exact_shuffle:loss=shuffle_comparison.floating(shuffle_comparison.structure(components,counts),[1.,*tilts])
    elif posterior_shuffle:
        loss=min(logq(capped_density_loss(G,caps)),
                 shuffle_comparison.posterior_floating(shuffle_comparison.posterior_structure(components,counts),[1.,*tilts]))
    elif anchors is None:loss=logq(capped_density_loss(G,caps))
    else:
        loss,logtau=density_tangent.floating(G,anchors)
        weights*=np.exp(logtau)
    scale=weights.sum()
    pieces=dict(
        outer=sum(counts[i]*logq(components[i][1]) for i in ids)+lgamma(G+1)-sum(lgamma(c+1) for c in counts),
        shuffle_loss=REGIONS*loss,
        input_reweighting=REGIONS*float(n@np.log(z))+PACKETS*np.log(scale),
        inner_moment=log_power(kernel.floating(data,weights/scale,lam)),
        cutoff=lam*threshold)
    return {key:float(value/np.log(2)) for key,value in pieces.items()}


def refine_density(data,components,counts,threshold,baseline):
    best=baseline;best_anchors=None;seen=set()
    witness=parameters(baseline)
    weights=selected.reference_weights(components,counts,witness)
    raw=tuple(int(round(G*w/weights.sum())) for w in weights)
    anchors=raw
    for _ in range(5):
        if anchors in seen:break
        seen.add(anchors)
        candidate=selected.selected(data,components,counts,threshold,density_anchors=anchors)
        if candidate['log2_upper']<best['log2_upper']:best,best_anchors=candidate,anchors
        anchors=selected.posterior_anchors(data,components,counts,parameters(candidate),anchors)
    # Very rare categories can disappear under the output tilt. A zero
    # anchor is valid even if the actual category count is positive.
    for anchor in list(seen):
        zeroed=tuple(0 if x<=4 else x for x in anchor)
        if zeroed in seen:continue
        seen.add(zeroed)
        candidate=selected.selected(data,components,counts,threshold,density_anchors=zeroed)
        if candidate['log2_upper']<best['log2_upper']:best,best_anchors=candidate,zeroed
    return best,best_anchors


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',nargs='+',choices=tuple(CASES),default=list(CASES))
    parser.add_argument('--distances',nargs='+',default=['1/100','1/50','1/20'])
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--density-tangents',action='store_true')
    parser.add_argument('--exact-shuffle',action='store_true',help='Special bound for one stochastic type plus disjoint constant categories')
    parser.add_argument('--posterior-shuffle',action='store_true',help='General bound retaining posterior packet-type count randomness')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--replay',type=Path)
    args=parser.parse_args()
    if sum((args.density_tangents,args.exact_shuffle,args.posterior_shuffle))>1:parser.error('choose one shuffle refinement')
    distances=list(map(Q,args.distances))
    if args.precision<128 or any(not 0<d<1 for d in distances):parser.error('precision >=128 and distances in (0,1) required')
    components=actual_components()
    if args.exact_shuffle and not args.replay:
        for name in args.cases:
            try:shuffle_comparison.structure(components,composition(components,CASES[name]))
            except ValueError as error:parser.error(f'{name}: {error}; select compatible --cases')
    data=kernel.actual();ctx.prec=args.precision
    if args.replay:
        record=json.loads(args.replay.read_text())
        if record.get('schema')!=SCHEMA:parser.error('unsupported record')
        for row in record['rows']:
            counts=composition(components,row['composition'])
            if row['threshold']!=int(Q(row['distance'])*N):raise ValueError('cutoff mismatch')
            upper=selected.outward(data,components,counts,row['threshold'],row['parameters'],
                                   density_anchors=row.get('density_anchors'),exact_shuffle=row.get('exact_shuffle',False),
                                   posterior_shuffle=row.get('posterior_shuffle',False))
            print('REPLAY',row['case'],row['distance'],'log2 upper',upper.log()/arb(2).log(),flush=True)
        print('Selected comparison compositions only; not a full-code certificate.',flush=True)
        return
    record=dict(schema=SCHEMA,precision=args.precision,rows=[])
    for name in args.cases:
        counts=composition(components,CASES[name])
        for distance in distances:
            start=monotonic();threshold=int(distance*N)
            proposal=selected.selected(data,components,counts,threshold,exact_shuffle=args.exact_shuffle,posterior_shuffle=args.posterior_shuffle)
            anchors=None;baseline=proposal['log2_upper']
            if args.density_tangents:proposal,anchors=refine_density(data,components,counts,threshold,proposal)
            witness=parameters(proposal)
            upper=selected.outward(data,components,counts,threshold,witness,density_anchors=anchors,exact_shuffle=args.exact_shuffle,posterior_shuffle=args.posterior_shuffle)
            logupper=upper.log()/arb(2).log()
            row=dict(case=name,composition=CASES[name],distance=str(distance),threshold=threshold,
                     parameters=list(map(str,witness)),outward_log2=str(logupper),
                     log2_upper=float(logupper),proposal=proposal['log2_upper'],baseline=baseline,
                     density_anchors=anchors,
                     exact_shuffle=args.exact_shuffle,
                     posterior_shuffle=args.posterior_shuffle,
                     terms=breakdown(data,components,counts,threshold,witness,anchors,args.exact_shuffle,args.posterior_shuffle),seconds=monotonic()-start)
            record['rows'].append(row)
            if args.output:
                args.output.parent.mkdir(parents=True,exist_ok=True)
                args.output.write_text(json.dumps(record,indent=2)+'\n')
            print('SELECTED',name,str(distance),'log2 upper',logupper,'seconds',round(row['seconds'],2),flush=True)
    print('Selected comparison compositions only; not a full-code certificate.',flush=True)


if __name__=='__main__':main()

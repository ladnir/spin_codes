"""Floating diagnostic: dominate lazy feedback by a uniform density.

Not a certificate. The candidate uses the same encoder; outward transfer
and full-domain replay are deliberately not implemented in this probe.
"""
import argparse
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
import numpy as np
import rank_return
import trimmed_return
import scalar_cover as sc


def prepare(data):
    W=data['windows'];denoms=[comb(W,j)*15**j for j in range(W+1)]
    budgets=[int(D*max(Q(a),Q(z))) for D,a,z in
             zip(denoms,data['nonzero_atom_caps'],data['zero_probabilities'])]
    outputs=[trimmed_return.output_counts(h) for h in data['histograms']]
    values=[[trimmed_return.lightest(row,budgets[j]) for j,row in enumerate(rows)] for rows in outputs]
    denominator=np.array(denoms,dtype=float)[None,:,None]
    return dict(data,any_fiber_float=np.array(values,dtype=float)/denominator,
                output_float=np.array(outputs,dtype=float)/denominator)


def density(data,p,z):
    W=data['windows'];L=(1<<data['bits'])-1
    powers=z**np.arange(4*W+1);factors=(1+z)**(4*np.arange(W+1)-data['bits'])
    values=np.minimum((data['rank_float']@powers)*factors[None,:],data['any_fiber_float']@powers)
    weights=np.array([comb(W,j)*p**j*(1-p)**(W-j) for j in range(W+1)])
    return float(weights@values.max(axis=0)),float(weights@(data['histogram_multiplicities']@values/L))


def candidates(data,p,z,original):
    maximum,mean=density(data,p,z);L=(1<<data['bits'])-1;alpha=2.**-data['updates']
    result={'original':original}
    for m,u in ((True,False),(False,True),(True,True)):
        candidate=original.copy()
        if m:candidate[1,1]=0.;candidate[1,2]+=alpha*L*maximum
        if u:candidate[2,1]=0.;candidate[2,2]+=alpha*L*mean
        result[f'convert_m{int(m)}_u{int(u)}']=candidate
    W=data['windows'];powers=z**np.arange(4*W+1)
    factors=(1+z)**(4*np.arange(W+1)-data['bits'])
    fibers=np.minimum((data['rank_float']@powers)*factors[None,:],data['any_fiber_float']@powers)
    moments=data['output_float']@powers
    weights=np.array([comb(W,j)*p**j*(1-p)**(W-j) for j in range(W+1)])
    for threshold in sorted(set([1,2,4,6,8,10,12,16,W])):
        if threshold>W:continue
        remaining=weights.copy();remaining[threshold:]=0
        converting=weights-remaining
        kept= moments@remaining
        maximum=float(converting@fibers.max(axis=0))
        mean=float(converting@(data['histogram_multiplicities']@fibers/L))
        for m,u in ((True,False),(False,True),(True,True)):
            candidate=original.copy()
            if m:candidate[1,1]=alpha*kept.max();candidate[1,2]+=alpha*L*maximum
            if u:candidate[2,1]=alpha*(data['histogram_multiplicities']@kept/L);candidate[2,2]+=alpha*L*mean
            result[f'from_j{threshold}_m{int(m)}_u{int(u)}']=candidate
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path,help='rank-return frontier probe output')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    source=json.loads(args.input.read_text());data=prepare(rank_return.actual());rows=[]
    for row in source['rows']:
        p=float(Q(row['activity']));lam=float(Q(row['witness']['parameters'][0]));z=np.exp(-lam)
        original=rank_return.floating(data,sc.probabilities(Q(row['activity'])),lam)
        base=sc.log_power(original)
        scores={name:float(row['proposal']+(sc.log_power(matrix)-base)/np.log(2))
                for name,matrix in candidates(data,p,z,original).items()}
        result=dict(x=row['x'],scores=scores);rows.append(result);print(json.dumps(result),flush=True)
    if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,rows=rows),indent=2)+'\n')


if __name__=='__main__':main()

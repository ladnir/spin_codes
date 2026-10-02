"""Diagnostic: retain a density cap for births at fixed packet occupancy.

Local caps use outward arithmetic. Candidate selection and any printed
whole-length score use floating arithmetic and are not certificates.
No existing verifier or production encoder uses this module.
"""
import argparse
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
import numpy as np
from flint import arb,arb_poly,ctx
import occupancy_birth_classes as base
import scalar_cover as sc


def feedback_caps(data,z):
    """Cap E[z^wt(X) 1{CX=b} | J=j] for every nonzero b."""
    if not 0<z<=1:raise ValueError('output weight in (0,1] required')
    W=data['windows'];S=1<<data['bits']
    values=[((1+z)**(4-r)*(1-z)**r-1)/15 for r in range(5)]
    powers=[[arb_poly([1,value])**j for j in range(W+1)] for value in values]
    polynomials=[]
    for profile in data['records']:
        polynomial=arb_poly([1])
        for row,n in zip(powers,profile):polynomial*=row[int(n)]
        polynomials.append(polynomial)
    result=[];packet=((1+z)**4-1)/15
    for j in range(W+1):
        transform=[polynomial[j]/comb(W,j) for polynomial in polynomials]
        # Any exact centering constant is valid at nonzero feedback.
        # Floating ordering merely proposes such a constant.
        order=sorted(range(len(transform)),key=lambda i:float(transform[i]))
        seen=0;center=arb(0)
        for i in order:
            seen+=int(data['multiplicities'][i])
            if 2*seen>=S:center=arb(float(transform[i]));break
        cap=base.up(sum((int(n)*abs(value-center) for n,value in
                        zip(data['multiplicities'],transform)),arb(0))/S)
        result.append(min(cap,base.up(packet**j)))
    return result


def candidate(local,caps,bits,cutoff):
    """Replace high-occupancy births by a uniform-density upper envelope."""
    W=len(local)-1;L=(1<<bits)-1
    if len(caps)!=W+1 or type(cutoff) is not int or not 1<=cutoff<=W+1:
        raise ValueError('matching caps and a valid positive cutoff required')
    result=[]
    for j,matrix in enumerate(local):
        value=matrix*arb(1)
        if j>=cutoff:
            for k in range(1,value.ncols()):value[0,k]=arb(0)
            value[0,2]=base.up(L*caps[j])
        result.append(value)
    return result


def arrays(local):
    size=local[0].nrows()
    return np.array([[[float(m[i,j]) for j in range(size)] for i in range(size)] for m in local])


def iid_scores(local,caps,bits,p,cutoffs):
    W=len(local)-1
    if not 0<=p<=1:raise ValueError('activity probability required')
    weights=np.array([comb(W,j)*p**j*(1-p)**(W-j) for j in range(W+1)])
    original=sc.log_power(np.tensordot(weights,arrays(local),axes=1))/np.log(2)
    rows=[]
    for cutoff in cutoffs:
        matrices=candidate(local,caps,bits,cutoff)
        score=sc.log_power(np.tensordot(weights,arrays(matrices),axes=1))/np.log(2)
        rows.append(dict(cutoff=cutoff,iid_log2=float(score),gain_bits=float(original-score)))
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--updates',type=int,choices=(2,3),default=3)
    parser.add_argument('--activities',nargs='+',default=['.175','.2'])
    parser.add_argument('--tilts',nargs='+',default=['.08','.10','.12','.14'])
    parser.add_argument('--cutoffs',nargs='+',type=int,default=[1,2,4,6,8,10,12,16,33])
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.precision<128 or any(Q(t)<=0 for t in args.tilts):parser.error('precision >=128 and positive tilts required')
    data=base.actual(args.updates);ctx.prec=args.precision;rows=[]
    for tilt in map(Q,args.tilts):
        z=(-base.aq(tilt)).exp();local=base.outward_at_z(data,z);caps=feedback_caps(data,z)
        for activity in map(Q,args.activities):
            scores=iid_scores(local,caps,data['bits'],float(activity),args.cutoffs)
            row=dict(activity=str(activity),tilt=str(tilt),scores=scores)
            rows.append(row);print(json.dumps(row),flush=True)
            if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,
                updates=args.updates,precision=args.precision,rows=rows),indent=2)+'\n')
    print('IID comparison only; no regional or whole-code certificate.',flush=True)


if __name__=='__main__':main()

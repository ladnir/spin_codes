"""Diagnostic: retain a density envelope through selected lazy updates.

The local changes have outward bounds. Floating iid comparisons only
select candidates; no certificate verifier uses this module.
"""
import argparse
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
import numpy as np
from flint import ctx
import occupancy_birth_classes as base
import scalar_cover as sc
from lazy_density import density_caps,candidate


def arrays(local):
    size=local[0].nrows()
    return np.array([[[float(m[i,j]) for j in range(size)] for i in range(size)] for m in local])


def iid_scores(data,local,z,caps,p,cutoffs):
    W=data['windows']
    if not 0<=p<=1:raise ValueError('activity probability required')
    weights=np.array([comb(W,j)*p**j*(1-p)**(W-j) for j in range(W+1)])
    original=sc.log_power(np.tensordot(weights,arrays(local),axes=1))/np.log(2)
    rows=[]
    for cutoff in cutoffs:
        matrices=candidate(data,local,z,caps,cutoff)
        score=sc.log_power(np.tensordot(weights,arrays(matrices),axes=1))/np.log(2)
        rows.append(dict(through=cutoff,iid_log2=float(score),gain_bits=float(original-score)))
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--updates',type=int,choices=(2,3),default=3)
    parser.add_argument('--activities',nargs='+',default=['.175','.2'])
    parser.add_argument('--tilts',nargs='+',default=['.10','.12','.14'])
    parser.add_argument('--through',nargs='+',type=int,default=[0,1,2,4,6,8,12,32])
    parser.add_argument('--one-packet',action='store_true',help='Sharpen the density cap with the checked exact one-packet census')
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.precision<128 or any(Q(t)<=0 for t in args.tilts):parser.error('precision >=128 and positive tilts required')
    data=base.actual(args.updates);ctx.prec=args.precision;rows=[]
    if args.one_packet:
        import single_packet
        images,columns,_=sc.kernel.maps()
        if columns!=data['columns']:raise ArithmeticError('feedback map mismatch')
    for tilt in map(Q,args.tilts):
        z=(-base.aq(tilt)).exp();local=base.outward_at_z(data,z)
        single=single_packet.census(images,columns,z) if args.one_packet else None
        caps=density_caps(data,z,single)
        for activity in map(Q,args.activities):
            row=dict(activity=str(activity),tilt=str(tilt),scores=iid_scores(data,local,z,caps,float(activity),args.through))
            rows.append(row);print(json.dumps(row),flush=True)
            if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,
                updates=args.updates,precision=args.precision,one_packet=args.one_packet,rows=rows),indent=2)+'\n')
    print('IID comparison only; no regional or whole-code certificate.',flush=True)


if __name__=='__main__':main()

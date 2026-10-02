"""Diagnostic: exact small-occupancy cancellation moments in class operators.

Integer counts and dyadic local bounds are checked, but the comparison
scores are floating and do not establish a whole-code certificate.
"""
import argparse
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
import numpy as np
from flint import arb,ctx
import occupancy_birth_classes as base
from occupancy_birth_density_probe import arrays
import return_moment
import scalar_cover as sc


refine=return_moment.refine_class_returns


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--updates',type=int,choices=(2,3),default=3)
    parser.add_argument('--through',type=int,choices=(1,2,3,4),default=3)
    parser.add_argument('--activities',nargs='+',default=['.175','.2'])
    parser.add_argument('--tilts',nargs='+',default=['.10','.12','.14'])
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.precision<128 or any(Q(t)<=0 for t in args.tilts):parser.error('precision >=128 and positive tilts required')
    data=base.actual(args.updates);images,columns,_=sc.kernel.maps();ctx.prec=args.precision
    if columns!=data['columns']:raise ArithmeticError('feedback map mismatch')
    census=return_moment.census(images,columns,data['bits'],args.through,verbose=True)
    print('EXACT JOINT RETURN census checked through',args.through,flush=True)
    W=data['windows'];rows=[]
    for tilt in map(Q,args.tilts):
        z=(-base.aq(tilt)).exp();local=base.outward_at_z(data,z);improved=refine(data,local,census,z)
        for p in map(Q,args.activities):
            if not 0<=p<=1:parser.error('activities must be probabilities')
            weights=np.array([float(comb(W,j)*p**j*(1-p)**(W-j)) for j in range(W+1)])
            old=sc.log_power(np.tensordot(weights,arrays(local),axes=1))/np.log(2)
            new=sc.log_power(np.tensordot(weights,arrays(improved),axes=1))/np.log(2)
            row=dict(activity=str(p),tilt=str(tilt),baseline_log2=old,refined_log2=new,gain_bits=old-new)
            rows.append(row);print(json.dumps(row),flush=True)
            if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,updates=args.updates,
                precision=args.precision,through=args.through,rows=rows),indent=2,allow_nan=False)+'\n')
    print('IID comparison only; not a whole-code certificate.',flush=True)


if __name__=='__main__':main()

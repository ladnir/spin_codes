"""Choose fixed per-occupancy convex witnesses, then replay outward.

The optimizer sees binary64 coefficients only to propose rational mixing
fractions. Rebuilding each complete coupled column with those fractions
preserves its local inequality for every input shape and incoming measure.
The final replay covers a selected homogeneous support vector, not a full
occupancy. Nothing from this module enters the production operator cache.
"""
import argparse
from fractions import Fraction as Q

import numpy as np
from scipy.optimize import minimize

from mass_density_screen import as_array, baseline, score
from mass_density import blend
from local_sensitivity import placement_tape, sensitivity
from occupancy_memory import Z, C
from local_family import build as build_families


def family(base, zero, density, fractions):
    fractions=np.asarray(fractions,dtype=float)
    if (fractions.shape!=(2,len(base)-1) or not np.isfinite(fractions).all()
            or np.any(fractions<0) or np.any(fractions>1)):
        raise ValueError('two finite convex fractions per positive local occupancy required')
    result=base.copy()
    for index,target,alternative in ((0,Z,zero),(1,C,density)):
        f=fractions[index,:,None]
        result[1:,:,target]=(1-f)*base[1:,:,target]+f*alternative[1:,:,target]
    return result


def objective(base, zero, density, fractions, q, p, terminal, **geometry):
    ops=family(base,zero,density,fractions)
    value,_,gradient=sensitivity(ops,q,p,terminal,return_gradient=True,**geometry)
    direction=np.array([np.sum(gradient[1:,:,target]*(alternative[1:,:,target]-base[1:,:,target]),axis=1)
                        for target,alternative in ((Z,zero),(C,density))])
    return value,direction


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=80)
    parser.add_argument('--support',type=int,default=200)
    parser.add_argument('--tilt',default='.064')
    parser.add_argument('--penalty',default='.9')
    parser.add_argument('--maximum',type=int,default=16)
    parser.add_argument('--iterations',type=int,default=40)
    parser.add_argument('--passes',type=int,default=2)
    parser.add_argument('--epoch-cache',help='Optional source-bound exact memo shared with mass_verify')
    feedback_options=parser.add_mutually_exclusive_group()
    feedback_options.add_argument('--spectral-feedback',type=int,default=0,choices=[0,7,8,9,10])
    feedback_options.add_argument('--exact-feedback',type=int,default=0,choices=[0,7,8,9,10])
    parser.add_argument('--joint-four',action='store_true')
    parser.add_argument('--density-through',type=int,default=0,choices=[0,7,8,9,10])
    args=parser.parse_args()
    if (not 1<=args.groups<=2048 or not 38<=args.support<256
            or Q(args.tilt)<=0 or not 0<Q(args.penalty)<=1
            or not 6<=args.maximum<=32 or args.iterations<1 or args.passes<1):
        parser.error('invalid parameters')
    print('FIXED-WITNESS SEARCH; only final selected-point replay is outward',flush=True)
    base,coefficients=build_families([args.tilt],args.penalty,maximum=max(16,args.maximum),
        exact_feedback=args.exact_feedback,spectral_feedback=args.spectral_feedback,
        joint_four=args.joint_four,density_through=args.density_through,directory=args.epoch_cache)[args.tilt]
    all_ones={j:Q(1) for j in range(1,args.maximum+1)}
    zero=blend(base,coefficients,all_ones,target=Z)
    density=blend(base,coefficients,all_ones,target=C)
    arrays=tuple(map(as_array,(base,zero,density)))
    caps=baseline.authenticated_caps()
    cdf=baseline.integer_cdf(baseline.weighted_cdf_upper(caps,1<<128,full_weight=1/Q(args.penalty)))
    shells=baseline.weighted_union_shells(caps,full_weight=1/Q(args.penalty)); shells[0]-=1
    count=min(Q(cdf[args.support]),shells[args.support])
    fractions=np.zeros((2,len(base)-1)); fractions[:,5:min(12,args.maximum)]=1
    initial=family(*arrays,fractions)
    best,p=score(placement_tape(initial,args.groups)[-1],args.groups,args.support,count,args.tilt)
    print('INITIAL mass-ZC-6-12 log2 score',best,'p',p,flush=True)
    bounds=[(0.,1.) if j<args.maximum else (0.,0.) for _ in range(2) for j in range(len(base)-1)]
    for pass_index in range(args.passes):
        def evaluate(vector):
            value,gradient=objective(*arrays,vector.reshape(fractions.shape),args.groups,p,baseline.TAIL_TERMINAL)
            return value,gradient.ravel()
        fit=minimize(evaluate,fractions.ravel(),jac=True,bounds=bounds,method='L-BFGS-B',
                     options={'maxiter':args.iterations,'ftol':1e-12,'gtol':1e-8})
        proposal=np.clip(fit.x.reshape(fractions.shape),0,1)
        candidate,new_p=score(placement_tape(family(*arrays,proposal),args.groups)[-1],
                              args.groups,args.support,count,args.tilt)
        if candidate<best:
            fractions,best,p=proposal,candidate,new_p
        print('PASS',pass_index+1,'log2 score',best,'p',p,'iterations',fit.nit,'status',fit.message,flush=True)
    rational=[{j+1:Q(round(float(f)*10**6),10**6) for j,f in enumerate(row)
               if j<args.maximum and round(float(f)*10**6)>0} for row in fractions]
    print('FIXED RATIONAL WITNESSES zero/density',
          [{j:str(f) for j,f in row.items()} for row in rational],flush=True)
    exact_ops=blend(base,coefficients,rational[0],target=Z)
    exact_ops=blend(exact_ops,coefficients,rational[1],target=C)
    exact=baseline.placement(exact_ops,rounding=baseline.rounded,maximum_groups=args.groups)
    label=args.penalty+':optimized-mass'
    if args.spectral_feedback:
        label+=f':spectral-{args.spectral_feedback}'
    if args.exact_feedback:
        label+=f':exact-{args.exact_feedback}'
    if args.joint_four:
        label+=':joint-4'
    if args.density_through:
        label+=f':density-{args.density_through}'
    args.probe_supports=[args.support]; args.target_bits=52
    baseline.shell_points(args,{(args.tilt,label):(exact,as_array(exact))},{label:cdf},{label:shells})


if __name__=='__main__':
    main()

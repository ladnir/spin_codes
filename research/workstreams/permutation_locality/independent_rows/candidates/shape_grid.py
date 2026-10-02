"""Fast selected-point grid over tilt, all-one penalty, and update count.

Exact local shape caps are regenerated or loaded from the source-bound memo.
All rescaling and global evaluation here are binary64 DIAGNOSTICS ONLY.
"""
import argparse
from fractions import Fraction as Q
import numpy as np
import attack_cache
import shape_potential as shape
import mass_density_screen as screen
from occupancy_memory import Z,F,M,C,U


def retarget(families,spectrum,tilt,rounds,old_rounds=2):
    if rounds is not None and (type(rounds) is not int or rounds<old_rounds):
        raise ValueError('cannot reduce update count')
    a=2.**(-old_rounds);b=0. if rounds is None else 2.**(-rounds)
    lazy=b/a;refresh=(1-b)/(1-a);sizes=np.array([spectrum[v] for v in sorted(spectrum)])
    result=[]
    for occupancy,old in enumerate(families):
        new=old.copy();new[:,1:,:]*=lazy;new[:,1:,U:U+5]=old[:,1:,U:U+5]*refresh
        if occupancy:
            for i in (M,9,10):new[:,i,Z]=old[:,i,Z]*refresh
            for i in (F,*range(U,U+5)):
                upper=np.min(old[:,i,U:U+5]/sizes,axis=1)
                new[:,i,Z]=lazy*old[:,i,Z]+(refresh-lazy)*upper
        else:
            for i,v in enumerate(sorted(spectrum),U):
                new[:,i,i]=refresh*old[:,i,i]-(refresh*a-b)*np.exp(-float(tilt)*v)
        if not np.isfinite(new).all() or (new<0).any():raise ValueError('invalid retargeted family')
        result.append(new)
    return result


def penalize(families,penalty,old_penalty='.9'):
    ratio=float(Q(penalty)/Q(old_penalty));result=[]
    for j,family in enumerate(families):
        if j and len(family)>1:
            catalog=sorted(shape.expected_shapes(j))
            if len(catalog)!=len(family):raise ValueError('incomplete shape family')
            scale=np.array([ratio**s.count(4) for s in catalog])[:,None,None]
        else:scale=max(1.,ratio**j)
        result.append(family*scale)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts',nargs='+',default=['.056','.072','.088'])
    parser.add_argument('--penalties',nargs='+',default=['.7','.8','.9','.95','1'])
    parser.add_argument('--rounds',type=int,nargs='+',default=[2,3,4,8])
    parser.add_argument('--groups',type=int,nargs='+',default=[96,128])
    parser.add_argument('--directory',default='tmp/four-bit-attack-cache')
    args=parser.parse_args()
    if any(not 0<Q(rho)<=1 for rho in args.penalties) or any(r<2 for r in args.rounds):
        parser.error('penalties in (0,1] and update counts >=2 required')
    data=attack_cache.build(args.tilts,directory=args.directory)
    caps=screen.baseline.authenticated_caps();counts={};best={}
    for rho in args.penalties:
        full=1/Q(rho)
        cdf=screen.baseline.integer_cdf(screen.baseline.weighted_cdf_upper(caps,1<<128,full_weight=full))
        counts[rho]=min(Q(cdf[200]),screen.baseline.weighted_union_shells(caps,full_weight=full)[200])
    for tilt,(base,feedback,details) in data.items():
        shaped=shape.shape_matrices(base,feedback,details,tilt,'.9',2)
        families=shape.as_families(base,shaped)
        for rounds in args.rounds:
            mixed=retarget(families,details['spectrum'],tilt,rounds)
            for rho in args.penalties:
                operators=np.array([f.max(axis=0) for f in penalize(mixed,rho)])
                region=screen.float_placement(operators,max(args.groups))
                for q in args.groups:
                    value,p=screen.score(region,q,200,counts[rho],tilt,cutoff=193986)
                    print('SHAPE GRID q/tilt/rho/R',q,tilt,rho,rounds,'log2',value,'p',p,flush=True)
                    key=q,rounds
                    if key not in best or value<best[key][0]:best[key]=(value,tilt,rho,p)
    print('BEST SELECTED HOMOGENEOUS PROPOSALS',best,flush=True)
    print('No complete support cover or global outward replay.')


if __name__=='__main__':main()

"""Stronger-mixing diagnostics and transition ablations for four-bit SPIN.

None denotes the ideal uniform-refresh limit, not an implemented encoder.
Deletion ablations are counterfactuals, never valid probability bounds.
"""
import argparse
from math import log

import numpy as np
from flint import arb
import attack_cache
import mass_density_screen as screen
from group_rank_one_verify import up
from occupancy_memory import Z,F,M,C,U


def retarget_matrix(old,occupancy,spectrum,tilt,old_rounds,new_rounds):
    """Retarget an UNSPLIT envelope with the baseline lazy/refresh semantics.

In particular M/L->Z must contain refresh only. A coupled zero-column
replacement invalidates this precondition even if M/L->C remains zero.
"""
    if (type(old_rounds) is not int or old_rounds<1
        or new_rounds is not None and (type(new_rounds) is not int or new_rounds<old_rounds)
        or old.nrows()!=11 or old.ncols()!=11):
        raise ValueError('eleven coordinates and nondecreasing integer update count required')
    if any(old[i,C]!=0 for i in (M,9,10)):
        raise ValueError('retarget before coupled-column blending; unsplit envelope required')
    if new_rounds==old_rounds:return old*1
    a=arb(2)**(-old_rounds);b=arb(0) if new_rounds is None else arb(2)**(-new_rounds)
    lazy=b/a;refresh=(1-b)/(1-a);levels=sorted(spectrum)
    new=old*1
    for i in range(1,11):
        for j in range(11):new[i,j]=up(old[i,j]*(refresh if U<=j<U+5 else lazy))
        if occupancy:
            if i in (M,9,10):new[i,Z]=up(refresh*old[i,Z])
            elif i==F or U<=i<U+5:
                refreshed=min(up(old[i,U+k]/spectrum[v]) for k,v in enumerate(levels))
                new[i,Z]=up(lazy*old[i,Z]+(refresh-lazy)*refreshed)
        elif U<=i<U+5:
            new[i,i]=up(refresh*old[i,i]-(refresh*a-b)*(-arb(tilt)*levels[i-U]).exp())
    if any(new[i,j]<0 for i in range(11) for j in range(11)):
        raise ValueError('negative retargeted coefficient')
    return new


def retarget(base,spectrum,tilt,old_rounds,new_rounds):
    return [retarget_matrix(t,j,spectrum,tilt,old_rounds,new_rounds) for j,t in enumerate(base)]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts',nargs='+',default=['.056','.072','.088'])
    parser.add_argument('--directory',default='tmp/four-bit-attack-cache')
    parser.add_argument('--groups',type=int,nargs='+',default=[96,128])
    parser.add_argument('--supports',type=int,nargs='+',default=[200])
    parser.add_argument('--rounds',nargs='+',default=['2','3','4','8','ideal'])
    parser.add_argument('--ablate',action='store_true')
    args=parser.parse_args()
    try:
        rounds_grid=[None if r=='ideal' else int(r) for r in args.rounds]
    except ValueError:
        parser.error('update counts must be integers or ideal')
    if (any(r is not None and not 2<=r<=32 for r in rounds_grid)
        or any(not 1<=q<=2048 for q in args.groups)
        or any(not 38<=u<=256 for u in args.supports)):
        parser.error('invalid rounds, group count, or support')
    data=attack_cache.build(args.tilts,directory=args.directory)
    caps=screen.baseline.authenticated_caps()
    cdf=screen.baseline.integer_cdf(screen.baseline.weighted_cdf_upper(caps,1<<128,full_weight=screen.Q(10,9)))
    shells=screen.baseline.weighted_union_shells(caps,full_weight=screen.Q(10,9))
    counts={u:min(screen.Q(cdf[u]),shells[u]) for u in args.supports}
    for tilt,(base,_,details) in data.items():
        for rounds in rounds_grid:
            ops=screen.as_array(retarget(base,details['spectrum'],tilt,2,rounds))
            region=screen.float_placement(ops,max(args.groups))
            for q in args.groups:
                for u,count in counts.items():
                    value,p=screen.score(region,q,u,count,tilt,cutoff=193986)
                    print('MIXING q/tilt/rounds',q,tilt,rounds,'u',u,'log2',value,'p',p,flush=True)
        if not args.ablate:continue
        original=screen.as_array(base)
        variants={}
        for start in (1,4,7,9):
            v=original.copy();v[start:,Z,Z]=0;variants[f'delete-ZZ-from-{start}']=v
        for name,source,target in (('CZ',C,Z),('CC',C,C),('FZ',F,Z),('FC',F,C)):
            v=original.copy();v[1:,source,target]=0;variants['delete-'+name]=v
        v=original.copy();v[1:,C,:]=0;variants['delete-C-all']=v
        for name,ops in variants.items():
            region=screen.float_placement(ops,max(args.groups))
            for q in args.groups:
                for u,count in counts.items():
                    value,p=screen.score(region,q,u,count,tilt,cutoff=193986)
                    print('INVALID ABLATION',name,'q/tilt',q,tilt,'u',u,'log2',value,'p',p,flush=True)


if __name__=='__main__':main()

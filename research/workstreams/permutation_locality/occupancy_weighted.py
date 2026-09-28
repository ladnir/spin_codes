"""Selected-point screen retaining total input weight across all active groups.

Exact integer outer shell/total caps; positive tilted operators. Floating
optimization and selected supports do not constitute a distance certificate.
"""
import argparse
from itertools import product
from math import comb,log

import numpy as np
from scipy.special import logsumexp
from flint import ctx,fmpz_poly

from bch_joint_support import authenticated_caps,support_caps
from basis_lattice import improve_caps
from shortened_bound import dimension_caps
from occupancy_count_gap import oa_refinement
from occupancy_memory import prepare,epoch_operators,rounded,TERMINAL
from occupancy_model import local_data,placement
from occupancy_screen import optimize
from joint_support import span


def weight_caps(spectrum,u,total_cap):
    """All ordered four-tuples with union support exactly u, by total weight."""
    polynomial=fmpz_poly([1]+spectrum[1:u+1])**4
    return [min(total_cap,int(polynomial[w])) if u<=w<=4*u else 0
            for w in range(4*u+1)]


def allocate(shells,total_cap,penalty):
    """Maximize sum A_w penalty^-w given upper shells and total count.

    The integer greedy allocation is exact for these relaxed constraints.
    Only evaluation of the exponential weights below uses binary64.
    """
    assert penalty>0 and total_cap>=0 and all(c>=0 for c in shells)
    order=range(len(shells)) if penalty>=1 else range(len(shells)-1,-1,-1)
    remaining=total_cap
    result=[]
    for w in order:
        count=min(remaining,shells[w])
        if count:
            result.append((w,count))
            remaining-=count
        if not remaining:
            break
    assert sum(c for _,c in result)==min(total_cap,sum(shells))
    return result


def log_weighted_cap(shells,total_cap,penalty):
    terms=[log(c)-w*log(penalty) for w,c in allocate(shells,total_cap,penalty)]
    return float(logsumexp(terms)) if terms else -float('inf')


def self_test():
    from fractions import Fraction
    checks=0
    for rows in ([1,2,4],[0b1111000,0b1100110,0b1010101]):
        words=span(rows)
        n=max(words).bit_length()
        spectrum=[sum(w.bit_count()==i for w in words) for i in range(n+1)]
        actual=[[0]*(4*n+1) for _ in range(n+1)]
        for tup in product(words,repeat=4):
            u=(tup[0]|tup[1]|tup[2]|tup[3]).bit_count()
            actual[u][sum(w.bit_count() for w in tup)]+=1
        for u in range(1,n+1):
            total=sum(actual[u])
            shells=weight_caps(spectrum,u,total)
            assert all(a<=b for a,b in zip(actual[u],shells))
            for rho in (Fraction(7,8),Fraction(1),Fraction(9,8)):
                upper=sum(c*rho**(-w) for w,c in allocate(shells,total,rho))
                observed=sum(c*rho**(-w) for w,c in enumerate(actual[u]))
                assert upper>=observed
                checks+=1
    print('Exact support/weight caps and weighted count allocation:',checks,'small-code checks',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=16,choices=range(4,33))
    parser.add_argument('--tilts',nargs='+',default=['.004','.0064','.008','.01','.0128'])
    parser.add_argument('--penalties',nargs='+',default=['.99','.995','1','1.0025','1.005','1.01','1.02'])
    parser.add_argument('--supports',nargs='+',type=int,default=[76,80,84,96,128,176,216])
    args=parser.parse_args()
    self_test()
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    oa,_=oa_refinement(caps)
    totals={u:sum(row[u] for row in oa) for u in args.supports}
    shells={u:weight_caps(spectrum,u,totals[u]) for u in args.supports}
    count_terms={(u,p):log_weighted_cap(shells[u],totals[u],float(p)) for u in args.supports for p in args.penalties}
    prepared=prepare(local_data(4))
    ctx.prec=192
    best={u:(float('inf'),None,None,None) for u in args.supports}
    for tilt in args.tilts:
        for penalty in args.penalties:
            ops=epoch_operators(prepared,tilt,maximum_groups=args.groups,pair_conditioned=True,
                                input_penalty=penalty if penalty!='1' else 1)
            regions=placement(ops,rounding=rounded)
            arrays=[np.array([[float(t[i,j]) for j in range(9)] for i in range(9)]) for t in regions]
            for u in args.supports:
                value,_=optimize(arrays,(u,)*args.groups,float(tilt),TERMINAL)
                score=(value+args.groups*count_terms[u,penalty]+log(comb(2048,args.groups)))/log(2)
                if score<best[u][0]:
                    best[u]=score,tilt,penalty,value/log(2)
            print('BINARY64 input-weight screen',tilt,penalty,'best',[(u,round(best[u][0],3)) for u in args.supports],flush=True)
    for u,row in best.items():
        print('RESULT u, log2 contribution, tilt, penalty, log2 weighted moment',u,*row,flush=True)
    print('Selected equal-support points only; OA counts require dual distance >=30; no full-code certificate.')


if __name__=='__main__':
    main()

"""Diagnostic refinement of fresh-state output moments, not a certificate.

Preserves the existing nine-coordinate invariant and all density/cancellation
bounds. Only fresh-to-mass and fresh-to-refresh entries are tightened.
"""
import argparse
from collections import Counter
from math import comb,log,log2

import numpy as np
from flint import arb,arb_mat,ctx

from group_moment import maps
from group_rank_one_verify import up
from occupancy_model import local_data,placement
from occupancy_memory import prepare,epoch_operators,rounded,shape_classes,Z,F,M,C,U,TERMINAL
from occupancy_screen import optimize
from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from occupancy_count_gap import oa_refinement


def fresh_census(prepared):
    images,columns,spectrum=maps()
    atoms=prepared[0][0][3]
    levels={a:Counter() for a in range(1,5)}
    moments={}
    inputs={b:[] for b in range(1,5)}
    for start in range(0,128,4):
        for mask in range(1,16):
            inputs[mask.bit_count()].append(mask<<start)
    for a in range(1,5):
        distribution=atoms[(a,)]
        assert sum(distribution.values())==32*comb(4,a) and not distribution[0]
        for state,count in distribution.items():
            levels[a][images[state].bit_count()]+=count
        for b,words in inputs.items():
            row=Counter()
            for state,count in distribution.items():
                for word in words:
                    row[(images[state]^word).bit_count()]+=count
            assert sum(row.values())==32*comb(4,a)*32*comb(4,b)
            moments[a,b]=row
    return spectrum,levels,moments


def refine(prepared,census,tilt,groups,full_penalty=1,input_penalty=1,odd_penalty=1):
    spectrum,levels,moments=census
    result=epoch_operators(prepared,tilt,maximum_groups=groups,pair_conditioned=True,full_penalty=full_penalty,input_penalty=input_penalty,odd_penalty=odd_penalty)
    rho=arb(full_penalty)
    powers=[up((-arb(tilt)*w).exp()) for w in range(129)]
    def average(hist,weight=0):
        return up(sum((count*powers[max(0,w-weight)] for w,count in hist.items()),arb(0))/sum(hist.values()))
    for j,t in enumerate(result):
        if j==1:
            moment=max(up(average(hist)*rho**int(b==4)*arb(input_penalty)**b*arb(odd_penalty)**(b%2)) for (a,b),hist in moments.items())
        elif j==0:
            moment=max(average(hist) for hist in levels.values())
        else:
            # Keep the all-one penalty coupled to the actual weight shape.
            by_weight={w:max(average(hist,w) for hist in levels.values())
                       for w in range(j,4*j+1)}
            moment=max(up(by_weight[w]*rho**full*arb(input_penalty)**w*arb(odd_penalty)**odd) for w,_,full,odd in shape_classes(j,True,True))
        if j:
            t[F,M]=min(t[F,M],up(moment/2))
        # Empty lazy transitions stay in F and retain their pointwise
        # multiplier. Replacing that entry by an average would be invalid.
        for k,w in enumerate(sorted(spectrum)):
            t[F,U+k]=min(t[F,U+k],up(moment*spectrum[w]/(2*((1<<19)-1))))
    return result


def self_test(prepared,census):
    spectrum,levels,moments=census
    checks=0
    def reference(hist,tilt,shift=0):
        # Reference balls need more precision than the tested outward
        # coefficient, or an exact upper endpoint can overlap their radius.
        saved=ctx.prec
        try:
            ctx.prec=saved+128
            return sum((count*(-arb(tilt)*max(0,w-shift)).exp() for w,count in hist.items()),arb(0))/sum(hist.values())
        finally:
            ctx.prec=saved
    for tilt in ('0','.0016','.01'):
        for penalty in ('1','.5'):
            operators=refine(prepared,census,tilt,4,penalty)
            original=epoch_operators(prepared,tilt,maximum_groups=4,pair_conditioned=True,full_penalty=penalty)
            assert operators[0][F,F]==original[0][F,F]
            cases=[(0,hist,arb(1)) for hist in levels.values()]
            cases += [(1,hist,arb(penalty)**int(b==4)) for (a,b),hist in moments.items()]
            for j,hist,scale in cases:
                actual=scale*reference(hist,tilt)
                if j:
                    assert actual/2<=operators[j][F,M]
                for k,w in enumerate(sorted(spectrum)):
                    assert actual*spectrum[w]/(2*((1<<19)-1))<=operators[j][F,U+k]
                checks+=1
            for j in range(2,5):
                for weight,_,full in shape_classes(j,True):
                    for hist in levels.values():
                        actual=arb(penalty)**full*reference(hist,tilt,weight)
                        assert actual/2<=operators[j][F,M],(tilt,penalty,j,weight,full,actual/2,operators[j][F,M])
                        checks+=1
            assert all(operators[j][a,b]<=original[j][a,b] for j in range(5) for a in range(9) for b in range(9))
    print('Fresh averaged-mass/refresh inequalities and unchanged lazy density:',checks,'checks',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=16,choices=range(4,33))
    parser.add_argument('--tilts',nargs='+',default=['.004','.005','.0064','.008','.01','.0128','.016'])
    parser.add_argument('--remove-mature-return',action='store_true',
                        help='INVALID AS A BOUND: suppress C-to-zero only to diagnose sensitivity')
    args=parser.parse_args()
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    shells,_=oa_refinement(caps)
    prepared=prepare(local_data(4))
    census=fresh_census(prepared)
    print('Exact fresh E-weight histograms:',census[1],flush=True)
    ctx.prec=192
    self_test(prepared,census)
    points=(72,76,80,83,84,96,128,176,216,256)
    best={u:(0.,None) for u in points}
    for tilt in args.tilts:
        refined=refine(prepared,census,tilt,args.groups)
        original=epoch_operators(prepared,tilt,maximum_groups=args.groups,pair_conditioned=True)
        assert all(refined[j][a,b]<=original[j][a,b] for j in range(args.groups+1) for a in range(9) for b in range(9))
        if args.remove_mature_return:
            for t in refined:
                t[C,Z]=arb(0)
        region=placement(refined,rounding=rounded)
        arrays=[np.array([[float(t[i,j]) for j in range(9)] for i in range(9)]) for t in region]
        for u in points:
            value,_=optimize(arrays,(u,)*args.groups,float(tilt),TERMINAL)
            if value<best[u][0]:
                best[u]=value,tilt
        print('BINARY64 fresh-moment tilt',tilt,flush=True)
    print('u original_counts OA_shell_counts inner_log2 tilt',flush=True)
    for u,(value,tilt) in best.items():
        constant=value/log(2)+log2(comb(2048,args.groups))
        print(u,constant+args.groups*log2(sum(row[u] for row in caps)),
              constant+args.groups*log2(sum(row[u] for row in shells)),value/log(2),tilt,flush=True)
    if args.remove_mature_return:
        print('COUNTERFACTUAL ONLY: mature lazy returns removed. These numbers are NOT probability bounds.')
    else:
        print('Selected binary64 points only; OA counts additionally require dual distance >=30.')


if __name__=='__main__':
    main()

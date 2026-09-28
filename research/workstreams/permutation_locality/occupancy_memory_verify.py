"""Complete support-box replay for a fixed occupancy of one-column groups.

Ideal uniform routing at K=2^20, IMT(128,19), bad output weight <=209715.
This checks only the requested number of active groups, not all messages.
"""
import argparse
from collections import Counter
from itertools import combinations_with_replacement
from math import comb,factorial,log

import numpy as np
from flint import arb,arb_mat,ctx

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from group_rank_one_verify import up
from occupancy_model import local_data
from occupancy_memory import prepare,build,self_test,C,TERMINAL
from occupancy_screen import optimize,matrix_for_probabilities
from two_group_screen import log_power_moment,log_binomial_mass

TILTS=('.00064','.001','.00125','.001375','.0016','.002','.0025')
DENOMINATOR=10**9


def buckets(caps,step):
    counts=[sum(row[u] for row in caps) for u in range(257)]
    starts=sorted(set(list(range(38,257,step))+[38,57,67,72,257]))
    result=[(lo,end-1,counts[end-1]) for lo,end in zip(starts,starts[1:])]
    assert [u for lo,hi,_ in result for u in range(lo,hi+1)]==list(range(38,257))
    return result


def endpoint_lower(probability,lo,hi):
    return min(arb((arb(comb(256,u))*probability**u*(1-probability)**(256-u)).lower()) for u in (lo,hi))


def replay(region,tilt,ps,box):
    probabilities=[arb(p)/DENOMINATOR for p in ps]
    masses=[arb(1)]
    for p in probabilities:
        updated=[arb(0)]*(len(masses)+1)
        for k,m in enumerate(masses):
            updated[k]+=m*(1-p)
            updated[k+1]+=m*p
        masses=updated
    matrix=sum((mass*r for mass,r in zip(masses,region)),arb_mat(9,9))
    powered=matrix**256
    value=sum((powered[0,j] for j in range(9) if j!=C),arb(0))*(arb(tilt)*209715).exp()
    for p,(lo,hi,_) in zip(probabilities,box):
        lower=endpoint_lower(p,lo,hi)
        assert lower>0
        value=up(value/lower)
    return min(arb(1),up(value))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,choices=(3,4),default=3)
    parser.add_argument('--step',type=int,default=8)
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--screen-only',action='store_true')
    args=parser.parse_args()
    assert args.precision>=128 and 1<=args.step<=32
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    intervals=buckets(caps,args.step)
    prepared=prepare(local_data(args.groups))
    ctx.prec=args.precision
    self_test(prepared)
    operators={}
    for tilt in TILTS:
        exact,approximate=build(prepared,tilt)
        probabilities=[]
        for lo,hi,_ in intervals:
            _,ps=optimize(approximate,((lo+hi)/2,)*args.groups,float(tilt),TERMINAL)
            probabilities.append(DENOMINATOR if lo==256 else max(1,min(DENOMINATOR-1,round(ps[0]*DENOMINATOR))))
        operators[tilt]=exact,approximate,probabilities
        print('Built',args.precision,'bit operators and fixed witnesses for tilt',tilt,flush=True)
    total=arb(0)
    float_terms=[]
    worst=[]
    covered=0
    location_count=comb(2048,args.groups)
    number=comb(len(intervals)+args.groups-1,args.groups)
    progress_step=max(500,(number//10//500)*500)
    for index,indices in enumerate(combinations_with_replacement(range(len(intervals)),args.groups),1):
        box=[intervals[i] for i in indices]
        multiplicity=factorial(args.groups)
        for repeat in Counter(indices).values():
            multiplicity//=factorial(repeat)
        covered+=multiplicity
        witness=None
        for tilt,(_,approximate,ps) in operators.items():
            nums=[ps[i] for i in indices]
            probabilities=[p/DENOMINATOR for p in nums]
            value=log_power_moment(matrix_for_probabilities(approximate,probabilities),256,TERMINAL)+float(tilt)*209715
            value-=sum(min(log_binomial_mass(256,u,p) for u in (lo,hi)) for p,(lo,hi,_) in zip(probabilities,box))
            if witness is None or value<witness[0]:
                witness=value,tilt,nums
        value,tilt,ps=witness
        count=location_count*multiplicity
        for _,_,cap in box:
            count*=cap
        score=min(0.,value)+log(count)
        float_terms.append(score)
        worst.append((score/log(2),[(lo,hi) for lo,hi,_ in box],tilt))
        if not args.screen_only:
            total=up(total+replay(operators[tilt][0],tilt,ps,box)*count)
        if index%progress_step==0:
            print('Support boxes',index,'/',number,flush=True)
    assert covered==len(intervals)**args.groups
    from scipy.special import logsumexp
    print('BINARY64 full-support box union log2',logsumexp(float_terms)/log(2),flush=True)
    print('Largest boxes:',*sorted(worst,reverse=True)[:10],sep='\n',flush=True)
    if not args.screen_only:
        print('OUTWARD',args.precision,'bits, occupancy',args.groups,'total upper',total,
              'margin',-total.log()/arb(2).log(),flush=True)
        assert 0<total<arb(2)**-45
        print('VERIFIED for this occupancy: all supports and all group locations, <2^-45',flush=True)
    print('Full-code certificate: NO. Other higher occupancies remain outside this replay.')


if __name__=='__main__':
    main()

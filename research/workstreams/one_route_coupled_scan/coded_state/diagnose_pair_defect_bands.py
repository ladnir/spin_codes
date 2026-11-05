"""Charge intervals for one defective episode with regional sharing.

Retains actual RS zeros. Bounds all charges in a band, not sampled charges.
Fixed construction, restricted geometry, no generated numerical files.
"""
import argparse
import itertools
import math
import time
from functools import lru_cache

from diagnose_pair_activation_route import M,B,Q,NEG,choose,logceil_ratio,episode_table
from diagnose_pair_defect_envelope import initial_masks,history_cap
from diagnose_pair_eligible_rows import eligible_majorant
from diagnose_pair_region_rank import count_cap
from diagnose_pair_shared_regions import sharing_price, case as shared_clean_case


@lru_cache(maxsize=128)
def band_table(active,length,low,high):
    assert 1<=low<=high<=B*length
    best=[NEG]*(2*active+1)
    births=sum(choose(B*length,n)*length**n for n in range(high+1))
    # A fixed initial support has only length-1 intermediate support masks.
    universal=1 << (B*(length-1))
    coefficients=[0]+[initial_masks(B,w,high)*min(length**w*births,universal)
                      * Q**(active-w) for w in range(1,active+1)]
    for left in range(1,length+1):
        right=length+1-left
        for nleft in range(1,active+1):
            for nright in range(active+1):
                numerator=sum(coefficients[w]
                              *choose((left-1)*min(B,w+high),nleft-w)
                              *choose(right*min(B,w+high),nright)
                              for w in range(1,nleft+1))
                if numerator:
                    denominator=choose(M,nleft)*choose(M,nright)*Q**(active+low)
                    deficit=2*active-nleft-nright
                    best[deficit]=max(best[deficit],logceil_ratio(numerator,denominator))
    return tuple(best)


def profile(low,high,shared=1):
    began=time.monotonic()
    lengths=(3,)*5
    dimension=3
    cap=count_cap(eligible_majorant(lengths,dimension+shared),256)
    terms=[]
    for active in range(8,cap+1):
        budget=(dimension+shared-1)*active
        dp=[0]+[NEG]*budget
        for index in range(5):
            table=band_table(active,3,low,high) if index==0 else episode_table(active,3)
            nxt=[NEG]*(budget+1)
            for used,previous in enumerate(dp):
                if previous==NEG:continue
                for extra,value in enumerate(table[:budget-used+1]):
                    if value!=NEG:nxt[used+extra]=max(nxt[used+extra],previous+value)
            dp=nxt
        envelope=tuple(itertools.accumulate(dp,max))
        zeros=[]
        for z in range(2*active+1):
            value=envelope[z+shared*active]
            if value!=NEG:
                zeros.append(((choose(131*active,z)-1).bit_length()
                              +32*(3*active-z)+value,z))
        if zeros:
            best=max(zeros)
            terms.append(((choose(M,active)-1).bit_length()+best[0]
                          +(len(zeros)-1).bit_length()+sharing_price(active,shared),active,best[1]))
        if time.monotonic()-began>35:
            raise TimeoutError('35-second defect-band cap; incomplete')
    worst=max(terms)
    # Density, source counts, endpoints, and choice of exceptional episode.
    upper=worst[0]+1+(len(terms)-1).bit_length()+156+3
    result=dict(charge_band=(low,high),shared=shared,lengths=lengths,regions=131,
                eligible_cap=cap,source_counts=len(terms),worst=worst,
                good_route_upper=upper,route_exception_upper=-100,
                family_upper=max(upper,-100)+1,whole_code_certificate=False,
                seconds=time.monotonic()-began)
    print(result,flush=True)
    return result


def test():
    # Charge-band envelopes dominate each exact-charge numerator factor.
    checks=0
    for width,length in ((4,3),(8,2)):
        for low,high in ((1,1),(1,2),(2,4)):
            universal=1 << (width*(length-1))
            for w in range(1,width+1):
                count=min(history_cap(width,length,w,high),universal)
                masks=initial_masks(width,w,high)
                for exact in range(low,high+1):
                    assert count>=min(history_cap(width,length,w,exact),universal)
                    assert masks>=initial_masks(width,w,exact)
                    assert min(width,w+high)>=min(width,w+exact)
                    assert low<=exact
                    checks+=1
    print(dict(charge_band_component_checks=checks),flush=True)


def replay():
    began=time.monotonic()
    bands=((1,2),(3,8),(9,32),(33,128),(129,1536))
    assert bands[0][0]==1 and bands[-1][1]==3*B
    assert all(right+1==nxt for (_,right),(nxt,_) in zip(bands,bands[1:]))
    clean=shared_clean_case((3,)*5,3,1)
    bounds=[clean['good_route_upper']]
    for low,high in bands:
        result=profile(low,high,1)
        assert result['good_route_upper']<=-100
        bounds.append(result['good_route_upper'])
    good=max(bounds)+(len(bounds)-1).bit_length()
    print(dict(kind='one_unrestricted_closed_episode',charge_cover=(0,3*B),
               regions=131,shared=1,lengths=(3,)*5,
               good_route_union_upper=good,route_exception_once=-100,
               family_upper=max(good,-100)+1,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--low',type=int,default=1)
    p.add_argument('--high',type=int,default=2)
    p.add_argument('--shared',type=int,default=1)
    p.add_argument('--test',action='store_true')
    p.add_argument('--replay',action='store_true')
    args=p.parse_args()
    assert 1<=args.low<=args.high<=1536 and 0<=args.shared<=1
    if args.test:test()
    elif args.replay:replay()
    else:profile(args.low,args.high,args.shared)

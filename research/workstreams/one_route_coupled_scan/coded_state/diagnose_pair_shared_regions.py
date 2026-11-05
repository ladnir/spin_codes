"""Shared-region correction to the clean activation-route count.

Split overlapping regional fragments only for the bound, never construction.
Charge exact finite-population dependence and allocation choices explicitly.
"""
import argparse
import itertools
import math
import time
from fractions import Fraction as F

from diagnose_pair_activation_route import (
    M, NEG, choose, logceil_ratio, episode_table)
from diagnose_pair_eligible_rows import eligible_majorant
from diagnose_pair_region_rank import count_cap


def sharing_price(active, shared):
    # Each shared region has at most two interval fragments. At fixed
    # allocation n1+n2=n<=s, the true/virtual denominator ratio is at most
    # 2^s*(M/(M-s+1))^s. There are at most s+1 allocations, which must be
    # summed (they are random), not merely maximized.
    return shared*(active+logceil_ratio(M**active,(M-active+1)**active)
                   + active.bit_length())


def deficit_envelope(active,lengths,budget):
    dp=[0]+[NEG]*budget
    for length in lengths:
        table=episode_table(active,length)
        nxt=[NEG]*(budget+1)
        for used,previous in enumerate(dp):
            if previous==NEG:continue
            for extra,value in enumerate(table[:budget-used+1]):
                if value!=NEG:
                    nxt[used+extra]=max(nxt[used+extra],previous+value)
        dp=nxt
    return tuple(itertools.accumulate(dp,max))


def case(lengths, dimension, shared):
    began=time.monotonic()
    virtual_dimension=dimension+shared
    placement=26*(len(lengths)+1)
    cap=count_cap(eligible_majorant(lengths,virtual_dimension),placement+100)
    minimum=1<<max(lengths)
    terms=[]
    for active in range(minimum,cap+1):
        if time.monotonic()-began>35:
            raise TimeoutError('35-second shared profile cap; incomplete')
        envelope=deficit_envelope(active,lengths,(virtual_dimension-1)*active)
        zero_terms=[]
        for zeros in range((dimension-1)*active+1):
            joint=envelope[shared*active+zeros]
            if joint!=NEG:
                # Fix the actual zero coordinates before routing. Each
                # specified zero costs one RS dimension. Splitting g regions
                # adds exactly g*s virtual zeros, but no message dimensions.
                zero_terms.append(((choose((128+dimension)*active,zeros)-1).bit_length()
                                   +32*(dimension*active-zeros)+joint,zeros))
        if zero_terms:
            best=max(zero_terms)
            terms.append(((choose(M,active)-1).bit_length()+best[0]
                          +(len(zero_terms)-1).bit_length()
                          +sharing_price(active,shared),active,best[1]))
    if terms:
        worst=max(terms)
        upper=1+worst[0]+(len(terms)-1).bit_length()+placement
    else:
        worst=None
        upper=NEG
    return dict(lengths=lengths,regions=128+dimension,shared=shared,
                virtual_regions=128+virtual_dimension,eligible_cap=cap,
                source_counts=len(terms),worst=worst,good_route_upper=upper,
                route_exception_upper=-100,family_upper=max(upper,-100)+1,
                whole_code_certificate=False)


def run(dimension=3,shared=1,cover=False,maximum=8):
    began=time.monotonic()
    if not cover:
        result=case((3,3,3,3,3),dimension,shared)
        print(dict(result,seconds=time.monotonic()-began),flush=True)
        return
    worst=NEG; witness=None; cases=excluded=0
    for lengths in itertools.combinations_with_replacement(range(1,maximum+1),5):
        result=case(lengths,dimension,shared)
        cases+=1
        excluded+=result['source_counts']==0
        if result['good_route_upper']>worst:
            worst,witness=result['good_route_upper'],result
        if time.monotonic()-began>40:
            print(dict(incomplete=True,processed=cases,worst=witness),flush=True)
            raise TimeoutError('40-second shared-region cover; incomplete')
    print(dict(regions=128+dimension,shared=shared,maximum_length=maximum,
               cases=cases,excluded=excluded,
               good_route_upper=worst,route_exception_upper=-100,
               family_upper=max(worst,-100)+1,worst=witness,
               whole_code_certificate=False,seconds=time.monotonic()-began),flush=True)


def test():
    began=time.monotonic()
    ratios=images=zeros=0
    for population in range(2,17):
        for active in range(1,population):
            bound=2**active*F(population,population-active+1)**active
            for left in range(active+1):
                for right in range(active-left+1):
                    ratio=F(math.comb(population,left)*math.comb(population,right),
                            math.comb(population,left+right))
                    assert ratio<=bound
                    ratios+=1
    # One real sample in disjoint fragments versus two virtual samples.
    # Sum over allocations before comparing; the population has eight slots.
    for a in range(5):
        for b in range(5):
            if a+b>8:continue
            for count in range(9):
                direct=sum(all(x<a or a<=x<a+b for x in sample)
                           for sample in itertools.combinations(range(8),count))
                convolution=sum(choose(a,left)*choose(b,count-left)
                                for left in range(count+1))
                assert direct==convolution==choose(a+b,count)
                images+=1
    # Exact GF4 RS[4,2] has 3 nonzero constants (zero zeros) and 12
    # nonconstant affine words (one zero). Convolve over active rows.
    for active in range(1,7):
        for total in range(active+1):
            exact=math.comb(active,total)*12**total*3**(active-total)
            bound=math.comb(4*active,total)*4**(2*active-total)
            assert exact<=bound
            zeros+=1
    # Splitting a region into two virtual copies preserves total nonzero
    # count and adds exactly s virtual zero coordinates.
    for active in range(1,17):
        for left in range(active+1):
            for right in range(active-left+1):
                assert 2*active-left-right == active+(active-left-right)
                zeros+=1
    print(dict(finite_population_ratio_checks=ratios,allocation_image_checks=images,
               zero_class_and_split_checks=zeros,seconds=time.monotonic()-began),flush=True)
    assert time.monotonic()-began<10


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--test',action='store_true')
    parser.add_argument('--cover',action='store_true')
    parser.add_argument('--shared',type=int,default=1)
    parser.add_argument('--dimension',type=int,default=3)
    parser.add_argument('--maximum-length',type=int,default=8)
    args=parser.parse_args()
    assert 1<=args.shared<=5 and 1<=args.dimension<=3
    assert 1<=args.maximum_length<=8
    if args.test:test()
    else:run(args.dimension,args.shared,args.cover,args.maximum_length)

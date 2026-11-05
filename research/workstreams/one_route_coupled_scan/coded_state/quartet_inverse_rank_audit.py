"""Exact checks for the all-row-count inverse rank bound. No full cert.

Exhaustive GF4 outputs for two RS[4,2] rows and two width4 updates.
Raw active-quartet point law is modeled by summing compatible outer words,
not by replacing the sampled code in a claimed distance certificate.
File-free, bounded180s. No native timing or construction changes.
"""
from collections import defaultdict
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time


def mul(a,b):
    p=0
    for _ in range(2):
        if b&1:p^=a
        b>>=1;a<<=1
        if a&4:a^=7
    return p


def inverse_test():
    # Route: each coordinate-region contains one coordinate of each RS row.
    code=[tuple(a^mul(b,j) for j in range(4)) for a,b in itertools.product(range(4),repeat=2)]
    outer=[tuple(v for j in range(4) for v in (left[j],right[j])) for left,right in itertools.product(code,repeat=2)]
    projections={mask:defaultdict(int) for mask in (0,240)}
    for mask,table in projections.items():
        for word in outer:table[tuple(word[i] for i in range(8) if not (mask>>i&1))]+=1
    totals=defaultdict(int);meta={};rejected=0
    for y in itertools.product(range(4),repeat=8):
        if not any(y):continue
        state=y[:4]+tuple(y[4+i]^mul(2,y[i]) for i in range(4))
        s=sum((v!=0)<<i for i,v in enumerate(state));j=sum((v!=0)<<i for i,v in enumerate(y))
        live=240 if any(y[:4]) else 0;births=j&~live;covered=live|births
        deaths=((s&15)&~(s>>4)).bit_count();weight=j.bit_count()
        rank=sum(min(2,sum(bool(covered>>(2*c+r)&1) for c in range(4))) for r in range(2))
        count=projections[live].get(tuple(y[i] for i in range(8) if not (live>>i&1)),0)
        impossible=any(0<sum(bool(births>>(2*c+r)&1) for c in range(4))
                       and sum(bool(covered>>(2*c+r)&1) for c in range(4))<=2 for r in range(2))
        if impossible:assert count==0;rejected+=1
        # Fixed common denominator q^8 for the sum of raw probabilities.
        totals[s,j]+=count*4**(8-live.bit_count())
        meta[s,j]=(weight,deaths,rank)
    for key,total in totals.items():
        weight,deaths,rank=meta[key]
        assert total<=4**(weight-deaths-rank+8)
    print(dict(exhaustive_nonzero_GF4_outputs=4**8-1,mask_classes=len(totals),
               exact_joint_rank_comparisons=len(totals),short_row_rejections=rejected),flush=True)


def support_test():
    checked=0
    for length in range(1,8):
        paths={(0,0,0):1};patterns=defaultdict(int)
        for _ in range(length):
            after=defaultdict(int)
            for (old,j,d),ways in paths.items():
                for new in (0,1):
                    outputs=(0,1) if old and new else (int(old!=new),)
                    for out in outputs:after[new,j+out,d+int(old and not new)]+=ways
            paths=after
        for (end,j,d),ways in paths.items():patterns[j,d]+=ways
        for (j,d),ways in patterns.items():
            assert ways<=math.comb(length,j)*math.comb(j,d);checked+=1
        # Multiple independent state coordinates, not a terminal reset.
        joint={(0,0):1}
        for width in range(1,5):
            after=defaultdict(int)
            for (j,d),ways in joint.items():
                for (jj,dd),value in patterns.items():after[j+jj,d+dd]+=ways*value
            joint=dict(after)
            for (j,d),ways in joint.items():
                assert ways<=math.comb(length*width,j)*math.comb(j,d);checked+=1
    print(dict(exact_output_death_support_counts=checked),flush=True)


def certificate_checks():
    q=1<<32;N=1<<21;K=1<<20;W=943718
    assert W+65538==1009256 and K-(W+65538)==39320
    # (1+1/(2q))^N <=1/(1-N/(2q))<2 via binomial domination.
    correction=1/(1-F(N,2*q))
    assert 2*F(q)**(-65538)*2**N*correction<F(2)**-62
    # Half-death envelope: (1+2^-17)^N<e^16<(11/4)^16<2^24.
    assert N==(1<<17)*16 and F(11,4)**16<2**24
    assert 2*F(q)**(-65539)*2**N*F(11,4)**16<F(2)**-71
    # Every (w,d) profile receives <=2^-107; sum <2^-63.
    assert (N+1)**2<2**44
    assert F(2)**-107*(N+1)**2<F(2)**-63
    print(dict(high_rank_threshold=W+65538,high_rank_failure_less_than='2^-62',
               half_death_threshold=W+65539,half_death_failure_less_than='2^-71',
               all_output_death_profiles_failure_less_than='2^-63'),flush=True)


def census_test():
    M=4;regions=4;cap=2;checks=0
    for sizes in ((0,0,4,4),(2,2,2,2),(1,1,3,3),(0,2,3,4),(0,0,0,0),(4,4,4,4)):
        subsets=[list(itertools.combinations(range(M),n)) for n in sizes]
        counts=defaultdict(int)
        for sets in itertools.product(*subsets):
            occupancy=[sum(i in group for group in sets) for i in range(M)]
            counts[sum(max(0,cap-v) for v in occupancy)]+=1
        denominator=math.prod(map(len,subsets));total=sum(sizes)
        distribution=[F(1)];conditioning=F(1)
        for n in sizes:
            p=F(n,M);new=[F(0)]*(len(distribution)+1)
            for k,v in enumerate(distribution):new[k]+=v*(1-p);new[k+1]+=v*p
            distribution=new;conditioning*=math.comb(M,n)*p**n*(1-p)**(M-n)
        assert conditioning>=F(1,M+1)**regions
        mean=F(total,M*regions)
        for t,rho in itertools.product((F(2),F(4)),(F(1,2),F(1),F(2))):
            actual=sum(F(v,denominator)*t**p for p,v in counts.items())
            marginal=sum(v*t**max(0,cap-k)*rho**k for k,v in enumerate(distribution))
            bound=marginal**M/rho**total/conditioning
            equalized=sum(math.comb(regions,k)*mean**k*(1-mean)**(regions-k)
                          *t**max(0,cap-k)*rho**k for k in range(regions+1))
            assert actual<=bound and marginal<=equalized;checks+=1
        if sizes==(0,0,4,4):assert counts=={0:1}
    print(dict(exact_conditioned_row_census_checks=checks,
               concentrated_half_region_control=True),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        inverse_test();support_test();certificate_checks();census_test()
        print(dict(seconds=time.monotonic()-start,full_distance_certificate=False),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

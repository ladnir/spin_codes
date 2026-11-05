"""Bounded diagnostic: does a global quota still lose closed suffix witnesses?

Common-marker or exact-global-count relaxations only; not the actual regional
ensemble. Source-only serial Peach run, <=180s, no files. FLOAT proposes
suffix witnesses; rational logarithm endpoints verify two fixed witnesses.
"""
import math
import os
from pathlib import Path
import shutil
import threading
import time
import itertools
from fractions import Fraction as F
from quartet_two_update import log_upper
from quartet_closed_paid_sets import canonical

N=1<<21;R=537600;CAP=255463
LAM=707/2560;BETA=6/5;D=31669.6


def lc(n,k):
    if not 0<=k<=n:return -math.inf
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1))/math.log(2)


def best_suffix(ell):
    best=(-math.inf,None,None)
    for length in range(512,N-R+CAP+1,512):
        lo=max(0,R-N+length);hi=min(length,R,CAP)
        if lo>hi:continue
        # Discrete successive-term ratio decreases with c.
        while lo<hi:
            c=(lo+hi)//2
            ratio=(length-c)*(R-c)/((c+1)*(N-length-R+c+1))
            if math.log2(ratio)>BETA*ell:lo=c+1
            else:hi=c
        c=lo
        value=lc(length,c)+lc(N-length,R-c)-255*lc(8192,2100)+ell*(LAM*length-BETA*c-D)
        if value>best[0]:best=value,length,c
    return best


def tests():
    full=(1<<16)-1;checked=0
    for length,c in ((4,1),(8,1),(8,2),(12,2),(12,3)):
        count=0
        for pp in itertools.combinations(range(16-length,16),c):
            paid=sum(1<<i for i in pp)
            for uu in itertools.combinations(range(16-length),4-c):
                unpaid=sum(1<<i for i in uu);inactive=paid|unpaid
                support,cover,score=canonical(full^unpaid)
                assert inactive&cover==paid and score>=2*length
                count+=1;checked+=1
        assert count==math.comb(length,c)*math.comb(16-length,4-c)
    print(dict(exact_closed_suffix_pairs=checked),flush=True)


def directed():
    lam=F(707,2560);beta=F(6,5);threshold=F(158348,5)
    for length,c,ell in ((910336,121999,F(3,2)),(1077760,94501,F(5,2))):
        assert length%512==0 and c<=CAP and length//8192*2100>c
        slope=lam*length-beta*c-threshold;assert slope>0
        lower=(-log_upper(-1,math.comb(length,c))-log_upper(-1,math.comb(N-length,R-c))
               -log_upper(255,math.comb(8192,2100))+ell*slope)
        assert lower>0
        compact=F(lower.numerator*(1<<32)//lower.denominator,1<<32)
        print(dict(suffix_length=length,paid=c,log2_b_at_least=str(ell),
                   directed_global_quota_log2_lower=str(compact),decimal_lower=float(compact),
                   exact_closed_witnesses=True,actual_regional_event_impossible=True,
                   minimum_actual_payments=length//8192*2100),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests()
        for ell in (0.,.125,.25,.5,1.,1.5,2.,2.5,3.):
            value,length,c=best_suffix(ell)
            print(dict(log2_b=ell,float_global_quota_lower=value,suffix_length=length,paid=c,
                       actual_regions_require_at_least=length//8192*2100,seconds=time.monotonic()-start),flush=True)
        directed()
    finally:timer.cancel()


if __name__=='__main__':main()

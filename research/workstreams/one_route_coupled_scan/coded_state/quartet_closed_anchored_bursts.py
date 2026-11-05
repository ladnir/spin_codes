"""Actual-region lower diagnostic for disjoint closed crossing-burst witnesses.

Vary burst lengths and payment counts; anchor each paid set in its last
512-position slab so different lengths are disjoint. No code-failure claim.
Serial180s, small arrays, no files. FLOAT locates candidates only.
"""
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np


def lc(n,k):
    if not 0<=k<=n:return -math.inf
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1))/math.log(2)


def logsum(values):
    peak=max(values)
    return peak+math.log2(sum(2**(v-peak) for v in values))


def anchor_log(v,c):
    whole=lc(v,c);old=lc(v-512,c)
    if c==0:return -math.inf
    if old==-math.inf:return whole
    return whole+math.log2(-math.expm1((old-whole)*math.log(2)))


def probe(ell):
    lam=707/2560;tau=381/512;d=31669.6;normal=lc(8192,2100);rows=[]
    # Each joined pair has at least four full updates, hence at least
    # three selected updates. Its positive run belongs to the canonical
    # maximizer, which is necessary for every payment to be closed.
    for half in range(2,16):
        v=512*half;terms=[]
        # 256*997 is below the surviving global payment cap255463.
        for c in range(max(1,2100-(8192-v)),min(997,v,2100)+1):
            terms.append(anchor_log(v,c)+lc(8192-v,2100-c)-normal+ell*(lam*v-6*c/5))
        if terms:rows.append((half,logsum(terms)))
    per_region=logsum([value for half,value in rows])
    lower=256*per_region+normal-ell*(d+512*128*(lam+tau))
    best=max(rows,key=lambda p:p[1])
    return lower,best,per_region


def tests():
    import itertools
    checks=0
    # One original-inactive set, different anchored prefixes produce
    # distinct paid sets. Sum of exact counts is therefore legitimate.
    for width,slab,r in ((8,2,3),(12,3,4)):
        for inactive in itertools.combinations(range(width),r):
            seen=set()
            for length in range(slab,width,slab):
                paid=tuple(i for i in inactive if i<length)
                if not any(i>=length-slab for i in paid):continue
                assert paid not in seen;seen.add(paid);checks+=1
        for length in range(slab,width,slab):
            for c in range(r+1):
                exact=sum(sum(i<length for i in inactive)==c and any(length-slab<=i<length for i in inactive)
                          for inactive in itertools.combinations(range(width),r))
                expected=(math.comb(length,c)-math.comb(length-slab,c))*math.comb(width-length,r-c) if c<=length and r-c<=width-length else 0
                assert exact==expected;checks+=1
    print(dict(exact_anchor_count_and_disjointness_checks=checks),flush=True)
    def canonical_cover(allowed):
        masks=[(allowed>>(4*t))&15 for t in range(8)]
        free=[masks[t] if masks[t+1]==15 else 0 for t in range(7)]+[masks[-1]];selected=[0]*8
        for bit in range(4):
            t=0
            while t<8:
                if not (free[t]>>bit)&1:t+=1;continue
                end=t
                while end<8 and (free[end]>>bit)&1:end+=1
                if end==8 or end-t>=3:
                    for j in range(t,end):selected[j]|=1<<bit
                t=end
        cover=[selected[t]|(15 if t and selected[t-1] else 0) for t in range(8)]
        score=2*sum(s.bit_count() for s in selected)-5*sum((selected[t]&~selected[t+1]).bit_count() for t in range(7))
        return sum(v<<(4*t) for t,v in enumerate(cover)),score
    masks=[sum(1<<i for i in p) for p in itertools.combinations(range(16),2)];checks=0;full=(1<<32)-1
    for left,right in itertools.product((2,3),repeat=2):
        burst=((1<<(4*(left+right)))-1)<<(16-4*left)
        for low,high in itertools.product(masks,repeat=2):
            inactive=low|(high<<16);paid=inactive&burst
            cover,score=canonical_cover(full^(inactive^paid))
            assert inactive&cover==paid and score>=8*(left+right-1)-20;checks+=1
    print(dict(exact_variable_length_crossing_burst_closure_checks=checks),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests()
        for ell in (1.,1.5,2.,2.25,2.5,2.75,3.):
            value,best,region=probe(ell)
            print(dict(log2_b=ell,float_actual_region_moment_lower=value,best_half=best[0],
                       sum_all_anchored_lengths=True,paid_per_region_at_most=997,
                       covered_counts_summed=True,negative_is_not_certificate=True,
                       seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

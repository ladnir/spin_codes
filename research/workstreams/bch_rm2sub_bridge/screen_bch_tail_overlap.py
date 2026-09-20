"""Bivariate OA29 Christoffel discovery for T intersect (T+d).

For a fixed difference d of weight c, split coordinates into its support
and complement. Total-degree-14 product Krawtchouk polynomials are mutually
orthogonal under uniform C because their products have degree <=28.
Point caps then sum over both tail-weight conditions. Floating SCREEN ONLY.
"""
import math
from pathlib import Path
import numpy as np
import bridge as base
from christoffel_caps import kraw
import screen_exponential_modes as cap_source


def basis(n):
    return np.array([[kraw(n,a,i)**2/math.comb(n,a) if a<=n else 0 for i in range(n+1)] for a in range(15)])


def run():
    caps,sources=cap_source.latest_caps();tail_lower=base.read(base.HERE/'generated/joint_tail_lower_80/lower.json')['lower']
    tail_upper=sum(caps[w] for w in range(38,81,2));rows=[]
    for c in range(38,161,2):
        first=basis(c);second=basis(256-c)
        kernel=sum(first[a,:,None]*np.sum(second[:15-a],axis=0)[None,:] for a in range(15))
        total=0.0
        for i in range(c+1):
            for j in range(257-c):
                a=i+j;b=c-i+j
                if not (38<=a<=80 and 38<=b<=80 and a%2==0 and b%2==0):continue
                total+=min(2.0**128/kernel[i,j],math.comb(c,i)*math.comb(256-c,j),caps[a],caps[b])
        total=min(total,float(tail_upper))
        pair_upper=min(float(tail_upper)**2,caps[c]*total)
        row=dict(difference_weight=c,intersection_upper_screen=total,
            intersection_bits=math.log2(total) if total else None,
            pair_fraction_upper_bits=math.log2(pair_upper)-2*math.log2(tail_lower) if pair_upper else None)
        rows.append(row)
    print('Tail overlap diagnostics',[(r['difference_weight'],round(r['intersection_bits'],3),round(r['pair_fraction_upper_bits'],3))
        for r in rows if r['difference_weight'] in [38,48,60,70,80,90,100,110,120,128,140,150,160]],flush=True)
    base.write_new(base.HERE/'generated/bch_tail_overlap_screen.json',dict(status='BIVARIATE_OA29_OVERLAP_SCREEN_ONLY',
        tail_lower=tail_lower,rows=rows,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),base.HERE/'christoffel_caps.py',
            base.HERE/'generated/joint_tail_lower_80/lower.json']+sources}))


if __name__=='__main__':run()

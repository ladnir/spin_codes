"""Exact toy matching moments and rational Chernoff-threshold checks."""
from collections import Counter
from fractions import Fraction as F
from math import comb
import time
from layered_matching_chernoff import matching_envelope,threshold,TILTS
from layered_route_gate import graph_envelope


def matchings(vertices):
    if not vertices:
        yield ();return
    a=vertices[0]
    for i in range(1,len(vertices)):
        b=vertices[i]
        for rest in matchings(vertices[1:i]+vertices[i+1:]):yield ((a,b),)+rest


def run():
    began=time.monotonic();moments=0;factorials=0;tails=0
    for M in (2,4,6,8,10):
        pairs=list(matchings(tuple(range(M))))
        for s in range(2,M//2+1):
            counts=Counter(sum(a<s and b<s for a,b in mat) for mat in pairs)
            n=s//2;p=F(s-1 if s%2==0 else s,M-1)
            for j in range(n+1):
                actual=sum(F(c*comb(z,j),len(pairs)) for z,c in counts.items() if z>=j)
                assert actual<=comb(n,j)*p**j;factorials+=1
            for t in TILTS:
                actual=sum(F(c,len(pairs))*t**z for z,c in counts.items())
                assert actual<=(1+p*(t-1))**n;moments+=1
    for n in (1,4,16,64,128):
        for p in (F(1,8),F(3,10),F(1,2),F(7,8)):
            for bits in (1,7,31):
                for t in TILTS:
                    k=threshold(n,p,bits,t)
                    if k<n:assert (1+p*(t-1))**n/t**k<=F(1,1<<bits)
                    # At k=n the violated event Z>k is empty.
                    assert 0<=k<=n;tails+=1
    old=graph_envelope(8192,193);new=matching_envelope(8192,193)
    assert all(0<=a<=b for a,b in zip(new,old))
    assert all(new[s]==193*(s-4096)+new[8192-s] for s in range(4097,8193))
    print(dict(exact_matching_moments=moments,factorial_moments=factorials,
        rational_thresholds=tails,improved_cardinalities=sum(a<b for a,b in zip(new,old)),
        peak_edge_old=old[2411],peak_edge_new=new[2411],elapsed_seconds=time.monotonic()-began))
    assert time.monotonic()-began<60


if __name__=='__main__':run()

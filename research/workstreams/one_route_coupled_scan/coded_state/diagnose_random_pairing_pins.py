"""Joint random-matching/pin transfer. Exact local checks, no certificate.

Source only; serial; <=30s tests/report. q is the comparison field size.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import math
import time
from diagnose_random_pairing import crossing_law,all_matchings,matchings,log2_fraction


def polynomial_mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def linear_power(a,b,n):
    return [math.comb(n,j)*a**(n-j)*b**j for j in range(n+1)]


def kernel(n,w,q,y,pin,beta):
    """Coefficient v: relaxed amplitude sum into next weight v.

    y marks output weight; pin prices a birth from an empty pair; beta
    prices an empty-pair coordinate remaining zero. Field point loss q^-D
    is included; the accepted-matrix density is NOT included here.
    """
    out=[F(0)]*(n+1)
    for k,prob in crossing_law(n,w).items():
        values=linear_power(y,1+q*y,w)
        values=polynomial_mul(values,linear_power(1,q*y,k))
        values=polynomial_mul(values,linear_power(beta,q*y*pin,n-w-k))
        price=prob/F(q**(w+k))
        out=[a+price*b for a,b in zip(out,values)]
    return out


def row_value(n,w,q,y,pin,beta,potential):
    """Evaluate kernel at a geometric state potential without expansion."""
    return sum((prob/F(q**(w+k))
                *(y+(1+q*y)*potential)**w
                *(1+q*y*potential)**k
                *(beta+q*y*pin*potential)**(n-w-k)
                for k,prob in crossing_law(n,w).items()),F(0))


def literal(n,w,q,y,pin,beta):
    out=[F(0)]*(n+1)
    active=(1<<w)-1
    for pairs in all_matchings(tuple(range(n))):
        live=0
        for a,b in pairs:
            if active&((1<<a)|(1<<b)):live|=(1<<a)|(1<<b)
        price=F(1,q**live.bit_count())/matchings(n)
        for nxt in range(1<<n):
            common=active&nxt; births=nxt&~active
            pins=(births&~live).bit_count()
            fixed_zero=n-(live|nxt).bit_count()
            one_output=(active^nxt).bit_count()
            # Each common active coordinate either has zero output, gaining
            # one equality, or nonzero output, retaining a q-valued variable.
            factor=F(q**births.bit_count())*y**one_output
            factor*=(1+q*y)**common.bit_count()
            factor*=pin**pins*beta**fixed_zero
            out[nxt.bit_count()]+=price*factor
    return out


def pin_identity_test():
    """Count temporal components literally for arbitrary changing matchings."""
    n=4; cases=0
    match=list(all_matchings(tuple(range(n))))
    # States Q1,Q2,Q3; Q0=0. Equality on common-active Q1,Q2 at a
    # deterministic selected subset; test all subsets to retain zeros jointly.
    for supports in product(range(16),repeat=3):
        common=[(t,i) for t in range(2) for i in range(n)
                if supports[t]&supports[t+1]&(1<<i)]
        # Bound exhaustive work: all equality subsets only for <=3 common
        # edges, and zero/full equality fixtures for larger edge sets.
        masks=range(1<<len(common)) if len(common)<=3 else (0,(1<<len(common))-1)
        for mask in masks:
            parent={(t,i):(t,i) for t in range(3) for i in range(n) if supports[t]>>i&1}
            def root(v):
                while parent[v]!=v:v=parent[v]
                return v
            z=0
            for j,(t,i) in enumerate(common):
                if mask>>j&1:parent[root((t+1,i))]=root((t,i));z+=1
            for m1,m2 in product(match,repeat=2):
                pinned=[root((0,i)) for i in range(n) if supports[0]>>i&1]
                for t,pairs in enumerate((m1,m2)):
                    live=0
                    for a,b in pairs:
                        if supports[t]&((1<<a)|(1<<b)):live|=(1<<a)|(1<<b)
                    pinned.extend(root((t+1,i)) for i in range(n)
                                  if supports[t+1]>>i&1 and not live>>i&1)
                components={root(v) for v in parent}
                assert len(components)==len(parent)-z
                assert len(pinned)==len(set(pinned))
                assert len(components)-len(pinned)==len(parent)-z-len(pinned)
                cases+=1
    return cases


def test():
    began=time.monotonic();checks=0
    for n in (2,4,6):
        for w in range(n+1):
            for q,y,pin,beta in ((4,F(1,2),F(1,8),F(3,4)),
                                  (8,F(2,3),F(1,3),F(4,5))):
                predicted=kernel(n,w,q,y,pin,beta)
                assert predicted==literal(n,w,q,y,pin,beta)
                h=F(3,2)
                assert sum(v*h**i for i,v in enumerate(predicted))==row_value(n,w,q,y,pin,beta,h)
                checks+=1
    pins=pin_identity_test()
    assert time.monotonic()-began<30
    print(dict(exact_kernel_rows=checks,arbitrary_matching_pin_histories=pins,
               seconds=time.monotonic()-began),flush=True)


def report():
    began=time.monotonic();n=512;w=256;q=1<<32;p=F(1,2)
    pin=p/q;beta=1-p+p/q;law=crossing_law(n,w)
    # Factor out the common q^-w inverse price. No output-amplitude sum
    # in this diagnostic: it isolates how the next-support pins tilt matching.
    for births in (0,256):
        weights={}
        for k,prob in law.items():
            j=k if births else 0
            weights[k]=prob/F(q**k)*pin**(births-j)*beta**(n-w-k-births+j)
        total=sum(weights.values());mean=sum(k*v for k,v in weights.items())/total
        print(dict(current_weight=w,next_inactive_births=births,
                   weighted_mean_crossing=float(mean),no_cross_fraction=float(weights[0]/total),
                   transform_log2=log2_fraction(total),comparison_only=True),flush=True)
    print(dict(seconds=time.monotonic()-began,whole_code_certificate=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--test',action='store_true')
    args=parser.parse_args()
    if args.test:test()
    else:report()

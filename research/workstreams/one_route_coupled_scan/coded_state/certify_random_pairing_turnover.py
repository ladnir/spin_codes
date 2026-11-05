"""Directed replay for <=256 swaps, constant weight256, two suffixes.

No whole-code claim. Exact rational moments and integer logarithm bounds.
Serial <=30s, memory bounded by short polynomial/integer powers, no files.
"""
from fractions import Fraction as F
import math
import time
from diagnose_random_pairing import crossing_law

SCALE=256


def logs(x):
    a=x.numerator**SCALE;b=x.denominator**SCALE
    d=a.bit_length()-b.bit_length()
    return d-1,d+1


def ceildiv(a,b):return -(-a//b)


def replay(R,cutoff,r,t):
    n=512;M=8192;w=256;T=16*R-1;limit=256
    x=F(1,1<<20)
    polynomial=sum((math.comb(w,b)*math.comb(n-w,b)*x**b for b in range(w+1)),F(0))
    hscaled=logs(F(math.comb(n,w)))[1]+T*logs(polynomial)[1]-limit*logs(x)[0]
    history=ceildiv(hscaled,SCALE)
    # Pad each history's union to this area, independently of routing.
    area=(T+1)*w+limit;p=F(area,R*M);center=area//M
    moment=sum((math.comb(R,b)*p**b*(1-p)**(R-b)
                *r**max(0,b-cutoff)*t**(center-b) for b in range(R+1)),F(0))
    budget=120+history
    # Independent Bernoulli surrogate, exact column-size conditioning.
    # At the mean integer column size, its mass is >=1/(M+1).
    numerator=(M*logs(moment)[1]+(area-center*M)*logs(t)[1]
               +R*logs(F(M+1))[1]+budget*SCALE)
    rlow=logs(r)[0];assert rlow>0
    excess=ceildiv(numerator,rlow)
    assert numerator-excess*rlow<=0
    epsilon=logs(1+F(2)**(R-4096+31*cutoff))[1]
    Vscaled=M*epsilon+31*excess*SCALE
    mw=max(prob/F(math.comb(n-w,k)) for k,prob in crossing_law(n,w).items())
    W=460*(1<<21)//1000;zeros=max(0,(T+1)*w-W)
    outchoices=(math.comb(T*w,zeros)-1).bit_length()
    # Use only baseline required zeros and drop all death penalties; this
    # is uniform over B<=limit. Keep the worst matching factor 2^limit.
    good=((1+history+limit+outchoices-32*zeros)*SCALE
          +T*logs(mw)[1]+Vscaled)
    good_bits=ceildiv(good,SCALE)
    return dict(R=R,w=w,total_swaps_at_most=limit,history_upper_bits=history,
                cutoff=cutoff,r=str(r),t=str(t),excess_cap=excess,
                good_route_upper_bits=good_bits,route_union_upper_bits=-120,
                passes=good_bits<=-120,whole_code_certificate=False)


def main():
    began=time.monotonic()
    for args in ((192,120,F(3,2),F(1)),(256,124,F(3,2),F(4,3))):
        result=replay(*args);print(result,flush=True)
        assert result['passes']
        if time.monotonic()-began>25:raise TimeoutError('turnover replay budget')
    print(dict(seconds=time.monotonic()-began,joint_restricted_upper_bits=-118))


if __name__=='__main__':main()

"""Exact-integer conditional matching envelope and bounded critical-class test."""
from math import comb,isqrt
from fractions import Fraction
from itertools import permutations
import time
import numpy as np
from layered_mixed_gate import M,b,N,envelopes,mixed_dimension,no_empty_point


def ceil_root(num,den):
    root=isqrt(num//den)
    return root+(root*root*den<num)


def joint_edge_cap(R,s,h,bits):
    """Conditional mean plus a sampling-without-replacement range bound."""
    n=M//2
    if s==0:return 0
    if s==M:return R*n
    # Sum h_j(s-h_j)/n <= (s*h-h^2/R)/n.
    mean=(R*s*h-h*h+R*n-1)//(R*n)
    # Reveal heads or tails, choosing the smaller TOTAL draw count d.
    # For d_j draws, squared ranges <=d_j*(n+1-d_j)/(n-1), including d_j=0.
    # Sum and Jensen: <=((n+1)*R*d-d^2)/(R*(n-1)).
    d=min(h,R*s-h)
    numerator=(n+1)*R*d-d*d
    range_num=max(0,numerator)
    range_den=R*(n-1)
    # Also retain the elementary d-range bound.
    slack=min(ceil_root(7*d*(bits+77),20),
              ceil_root(7*range_num*(bits+77),20*range_den))
    return min(h,R*s-h,mean+slack)


def dimension(R,P,r,f,e,env):
    H,E,sizes=env;D=P+r+f-2*e;c=r+f-e
    best=0;witness=None;choose=1
    for s in range(M+1):
        if s:choose=choose*(M-s+1)//s
        # Below Rs/2 both the conditional mean and the head bound increase.
        # Above it the range bound used below is replaced by a constant maximum
        # over that interval, ensuring the decreasing branch is monotone.
        lo=(R*s+1)//2;hi=int(H[s])
        if hi<lo:continue
        # Constant slack at the symmetric draw count is safe for all h>=Rs/2
        # only when s<=M/2. For larger s, use the old unconditional bound.
        if s>M//2:
            degree=min(D,P-e+int(E[s]),hi+c,R*s)
            value=degree-128*s
        else:
            bits=choose.bit_length();n=M//2;dmax=R*s//2
            range_num=max(0,(n+1)*R*dmax-dmax*dmax)
            slack=min(ceil_root(7*dmax*(bits+77),20),
                      ceil_root(7*range_num*(bits+77),20*R*(n-1))) if dmax else 0
            def decreasing(h):
                mu=(R*s*h-h*h+R*n-1)//(R*n)
                return P-e+min(int(E[s]),mu+slack)
            left,right=lo,hi
            while left<right:
                mid=(left+right)//2
                if mid+c>=decreasing(mid):right=mid
                else:left=mid+1
            degree=max(min(D,h+c,decreasing(h),R*s)
                       for h in {lo,hi,left,max(lo,left-1)})
            value=degree-128*s
        if value>best:best=value;witness=s
    return min(D//2,max(0,R*M-N//2),max(0,best)),witness


def tests():
    cases=0
    for n in range(2,9):
        for draws in range(n+1):
            exact=sum((Fraction(n-draws,n-i))**2 for i in range(1,draws+1)) if draws<n else 0
            upper=Fraction(draws*(n+1-draws),n-1)
            assert exact<=upper and exact<=draws
            cases+=1
    laws=0
    for n in range(1,7):
        for heads in range(n+1):
            for tails in range(n+1):
                counts=[0]*(n+1)
                for pi in permutations(range(n)):
                    counts[sum(pi[i]<tails for i in range(heads))]+=1
                denom=sum(counts)
                for edges,count in enumerate(counts):
                    numerator=comb(tails,edges)*comb(n-tails,heads-edges) if 0<=edges<=tails and 0<=heads-edges<=n-tails else 0
                    assert Fraction(count,denom)==Fraction(numerator,comb(n,heads))
                    laws+=1
    return dict(range_checks=cases,conditional_hypergeometric_checks=laws)


def main():
    started=time.monotonic();print(tests())
    R=192;P=(R*M-b)//2;r=0;f=589632;e=368520;B=49136
    env=envelopes(R)
    old=mixed_dimension(R,P,0,r,f,e,env)
    new,s=dimension(R,P,r,f,e,env)
    print(dict(old_dimension=old,joint_dimension=new,maximizing_subset_size=s,
               point=no_empty_point(R,f,r,e,env,births=B,dimension=new),
               whole_code_certificate=False))
    assert time.monotonic()-started<20


if __name__=='__main__':main()

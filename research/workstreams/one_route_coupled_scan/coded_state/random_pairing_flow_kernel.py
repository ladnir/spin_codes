"""Exact small checks of the finite-population causal-flow kernel.

The formula applies to512 coordinates; enumeration is deliberately capped
at8. Actual regional target counts are retained, not treated as IID.
Standard library only, serial <=30s, no output files.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import math
import shutil
import tempfile
import time
from random_pairing_demand_flow import causal_pair,self_test

TYPES=tuple((a,b) for a in range(5) for b in range(a,5))
INDEX={pair:i for i,pair in enumerate(TYPES)}


@lru_cache(maxsize=1000)
def matching_profiles(histogram):
    n=sum(histogram);assert n<=8 and n%2==0
    if not n:return {(0,)*15:F(1)}
    a=next(i for i,x in enumerate(histogram) if x)
    rest=list(histogram);rest[a]-=1;out={}
    for b,choices in enumerate(rest):
        if not choices:continue
        rem=rest.copy();rem[b]-=1;index=INDEX[min(a,b),max(a,b)]
        for profile,weight in matching_profiles(tuple(rem)).items():
            key=list(profile);key[index]+=1;key=tuple(key)
            out[key]=out.get(key,F(0))+F(choices,n-1)*weight
    return out


def profile_probability(histogram,profile):
    n=sum(histogram)
    denominator=math.prod(range(1,n,2))
    for (a,b),e in zip(TYPES,profile):
        denominator*=math.factorial(e)*(2**e if a==b else 1)
    return F(math.prod(math.factorial(h) for h in histogram),denominator)


def matchings(vertices):
    if not vertices:
        yield ();return
    a=vertices[0]
    for i in range(1,len(vertices)):
        for tail in matchings(vertices[1:i]+vertices[i+1:]):yield ((a,vertices[i]),)+tail


def tag_probability(n,j,regional,remaining=None):
    if regional:
        m,r=remaining if remaining is not None else (n+4,n//2+1)
        assert m>=n and 0<=r<=m
        return F(math.comb(m-n,r-j),math.comb(m,r)) if 0<=r-j<=m-n else F(0)
    p=F(1,4)
    return p**j*(1-p)**(n-j)


def kernel(histogram,regional,remaining=None):
    n=sum(histogram);out={}
    for profile,weight in matching_profiles(histogram).items():
        pairs=[pair for pair,e in zip(TYPES,profile) for _ in range(e)]
        for tags in product(range(3),repeat=n//2):
            new=[0]*5;loss=0;ways=1;j=sum(tags)
            for (a,b),targets in zip(pairs,tags):
                na,nb,_,ell=causal_pair(a-2,b-2,targets)
                new[na+2]+=1;new[nb+2]+=1;loss+=ell
                ways*=math.comb(2,targets)
            mass=weight*ways*tag_probability(n,j,regional,remaining)
            key=(tuple(new),loss,j);out[key]=out.get(key,F(0))+mass
    return {k:v for k,v in out.items() if v}


def brute(histogram,regional,remaining=None):
    states=[i-2 for i,h in enumerate(histogram) for _ in range(h)];n=len(states)
    pairs=list(matchings(tuple(range(n))));out={}
    for matching in pairs:
        for tags in range(1<<n):
            new=[0]*5;loss=0;j=tags.bit_count()
            for a,b in matching:
                targets=((tags>>a)&1)+((tags>>b)&1)
                na,nb,_,ell=causal_pair(states[a],states[b],targets)
                new[na+2]+=1;new[nb+2]+=1;loss+=ell
            mass=tag_probability(n,j,regional,remaining)/len(pairs)
            key=(tuple(new),loss,j);out[key]=out.get(key,F(0))+mass
    return {k:v for k,v in out.items() if v}


def verify_delayed_pairing_gain(histogram):
    """Exact IID three-step loss mean; not an exponential tail bound."""
    p=F(1,4);p0=(1-p)**2;p1=2*p*(1-p);n=sum(histogram)
    law=kernel(histogram,False)
    a,b,c=histogram[:3]
    predicted=p0**2*(p0*(a+b+c)+p1*F(a*(a+b-1),n-1))
    # Given the first transition, loss after two further updates has
    # mean p0^2*(h[-2]+h[-1]); only the first matching is relevant here.
    direct=sum(mass*p0**2*(new[0]+new[1]) for (new,_,_),mass in law.items())
    assert direct==predicted
    if all(h%2==0 for h in histogram):
        fixed=p0**2*(p0*(a+b+c)+p1*a)
        gain=p0**2*p1*F(a*(n-a-b),n-1)
        assert fixed-predicted==gain
    return predicted


def falling(x,k):
    return math.prod(range(x-k+1,x+1)) if x>=k else 0


def verify_regional_delayed_gain(histogram,m=16,r=5):
    """Three-step means under actual fixed-count sampling within one region."""
    n=sum(histogram);assert m>=3*n
    a,b,c=histogram[:3]
    predicted=(F(falling(m-r,2),falling(m,2))*a,
               F(falling(m-r,4),falling(m,4))*(a+b),
               F(falling(m-r,6),falling(m,6))*(a+b+c)
               +F(2*r*falling(m-r,5),falling(m,6))*F(a*(a+b-1),n-1))
    states={(histogram,r):F(1)};remaining_m=m
    for step in range(3):
        moment=sum(mass*F(falling(remaining_m-tags,2),falling(remaining_m,2))*h[0]
                   for (h,tags),mass in states.items())
        assert moment==predicted[step]
        if step==2:break
        new={}
        for (h,tags),mass in states.items():
            for (out,_,used),chance in kernel(h,True,(remaining_m,tags)).items():
                key=(out,tags-used);new[key]=new.get(key,F(0))+mass*chance
        assert sum(new.values())==1
        states=new;remaining_m-=n
    return predicted


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();profiles=0;kernels=0
    print(self_test(),flush=True)
    for histogram in ((2,0,0,0,2),(0,1,2,1,0),(1,1,2,1,1),(2,1,2,1,2)):
        law=matching_profiles(histogram)
        assert sum(law.values())==1
        for e,probability in law.items():
            assert probability==profile_probability(histogram,e);profiles+=1
        for regional in (False,True):
            actual=kernel(histogram,regional)
            assert actual==brute(histogram,regional) and sum(actual.values())==1
            kernels+=1
        verify_delayed_pairing_gain(histogram)
        assert time.monotonic()-began<25
    verify_delayed_pairing_gain((2,0,0,0,2))
    for h in ((2,0,0,0,2),(1,1,1,1,0),(4,0,0,0,0)):
        verify_regional_delayed_gain(h)
    print(dict(exact_matching_profiles=profiles,exact_histogram_loss_tag_kernels=kernels,
               exact_delayed_pairing_gain_checks=5,
               exact_three_step_regional_mean_checks=9,
               seconds=time.monotonic()-began,whole_code_certificate=False))


if __name__=='__main__':main()

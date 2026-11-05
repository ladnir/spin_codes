"""Bounded fixed-construction audit of the 224/32 quartet family.

Positive integer coefficient extraction; no FFT or floating-point value is
trusted for the directed bound. No files are written. Working arrays <16 MiB;
fixed local enumeration, 128 polynomial factors, and a 60-second watchdog.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
from functools import lru_cache

import numpy as np
from diagnose_pairing_sufficiency import choices

S,W,D,L,M,R = 512,224,32,4096,8192,1974
K,N = 1<<20,1<<21
P=F(387,577)
ZERO=F(70947,593351)
HOLD,NEXT,DEATH=F(907217,98518),F(62443,278617),F(37217,974128)
H,X,Y=HOLD*NEXT,NEXT,DEATH/NEXT
THETA=F(3,16)
ROOTS={}
FROZEN_MARKERS=('583109/878204','69994/585015','8934667/981625','101103/460037','32227/860918')
FROZEN_DUAL=('3356808999274797/357705217375993856','8993/326442','-2345/999543','-21527/582869')
FROZEN_UPPER=F(-4578692708060433,43980465111040)


def root_upper(x):
    if x not in ROOTS:
        cube=x**THETA.numerator;num=cube.numerator<<(48*THETA.denominator);den=cube.denominator
        lo,hi=0,1<<52
        while hi**THETA.denominator*den<num:hi*=2
        while hi-lo>1:
            mid=(lo+hi)//2
            if mid**THETA.denominator*den>=num:hi=mid
            else:lo=mid
        ROOTS[x]=F(hi,1<<48)
        assert ROOTS[x]**THETA.denominator>=cube and F(hi-1,1<<48)**THETA.denominator<cube
    return ROOTS[x]


def natural_interval(x):
    z=(x-1)/(x+1);assert 0<=z<=F(1,3)+F(1,1<<79)
    value=2*sum(z**(2*j+1)/F(2*j+1) for j in range(32))
    return value,value+2*z**65/(65*(1-z*z))


@lru_cache(maxsize=64)
def log_interval(x):
    x=F(x);assert x>0
    if x.numerator&(x.numerator-1)==0 and x.denominator&(x.denominator-1)==0:
        exact=x.numerator.bit_length()-x.denominator.bit_length();return F(exact),F(exact)
    exponent=x.numerator.bit_length()-x.denominator.bit_length()
    mantissa=x/F(2)**exponent
    if mantissa<1:exponent-=1;mantissa*=2
    assert 1<=mantissa<2
    floor=mantissa.numerator*(1<<80)//mantissa.denominator
    lo=natural_interval(F(floor,1<<80))[0]
    hi=natural_interval(F(floor+1,1<<80))[1]
    ln2lo,ln2hi=natural_interval(F(2))
    lower=exponent+lo/ln2hi;upper=exponent+hi/ln2lo
    lo_int=lower.numerator*(1<<40)//lower.denominator
    hi_int=-((-upper.numerator*(1<<40))//upper.denominator)
    return F(lo_int,1<<40),F(hi_int,1<<40)


def log_upper(c,x):
    lo,hi=log_interval(F(x));return c*(hi if c>=0 else lo)


def local():
    coefficients=defaultdict(F);records=[]
    for a in range(5):
        before=(1<<a)-1
        for after,target in itertools.product(range(16),repeat=2):
            counts=choices(4,before,after,target,True)
            if not any(counts):continue
            d=(before&~after).bit_count();b=(after&~before).bit_count()
            t=target.bit_count();covered=t if a else (target&after).bit_count()
            active_covered=4-t if a else after.bit_count()
            poly=sum(c*ZERO**j for j,c in enumerate(counts))
            value=(math.comb(4,a)*P**t*(1-P)**(4-t)
                   *root_upper(poly/F(1<<(32*covered)))*H**a*X**b*Y**d)
            coefficients[a,b,d]+=value
            records.append((value,a,b,d,t,active_covered))
    total=sum(coefficients.values())
    return {key:value/total for key,value in coefficients.items()},total,records


def coefficient_upper(probabilities,power=128,target=(224,32,32),precision=31):
    """Outward fixed-point DP; truncation is exact for nonnegative degrees.

    Round each coefficient and each stage upward on a 2^-precision grid.
    Scalar max bounds prove every uint64 accumulation is below overflow.
    """
    scale=1<<precision
    integers={k:(v.numerator*scale+v.denominator-1)//v.denominator for k,v in probabilities.items()}
    shape=tuple(x+1 for x in target)
    assert math.prod(shape)*8*3 < 16*(1<<20)
    current=np.zeros(shape,dtype=np.uint64);current[0,0,0]=scale
    nxt=np.zeros_like(current)
    maximum=scale;coefficient_sum=sum(integers.values())
    for step in range(power):
        assert maximum*coefficient_sum+scale-1 < 1<<64
        nxt.fill(0)
        # Loop bounds are construction constants. The slices avoid allocating
        # arrays outside the required coefficient box.
        for (a,b,d),value in integers.items():
            if any(k>n for k,n in zip((a,b,d),target)):continue
            extent=tuple(min(n-k,4*step)+1 for k,n in zip((a,b,d),target))
            source=tuple(slice(0,e) for e in extent)
            dest=tuple(slice(k,k+e) for k,e in zip((a,b,d),extent))
            nxt[dest]+=current[source]*np.uint64(value)
        nxt+=np.uint64(scale-1);nxt>>=np.uint64(precision)
        maximum=(maximum*coefficient_sum+scale-1)//scale
        assert int(nxt.max())<=maximum
        current,nxt=nxt,current
    return F(int(current[target]),scale)


def atom_cap(service):
    """Fourier atom bound for total adjacent-overlap service >0."""
    assert service>0
    # On |t|<=1, sin(|t|/2)>=(23/24)|t|/2. Integrating the
    # resulting Gaussian and using pi>3 bounds the central integral
    # by 24/(23 sqrt(6s)). The normalized tail is at most
    # exp(-529s/1152) <= 2^-floor(529s/1152).
    square=F(96,529)/service
    if square>=1:return F(1)
    lo,hi=0,1<<40
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid*mid*square.denominator >= square.numerator*(1<<80):hi=mid
        else:lo=mid
    central=F(hi,1<<40);assert central*central>=square
    exponent=529*service/1152
    tail=F(1,1<<(exponent.numerator//exponent.denominator))
    return min(F(1),central+tail)


def tests():
    p={(0,0,0):F(1,2),(1,0,0):F(1,4),(0,1,0):F(1,4)}
    exact={(0,0,0):F(1)}
    for power in range(1,5):
        nxt=defaultdict(F)
        for a,x in exact.items():
            for b,y in p.items():nxt[tuple(u+v for u,v in zip(a,b))]+=x*y
        exact=nxt
        for target in ((1,1,0),(2,0,0),(0,0,0)):
            upper=coefficient_upper(p,power,target,precision=20)
            assert upper>=exact.get(target,F(0))
            assert upper-exact.get(target,F(0))<F(1,10000)
    # Exact two-step support composition on n=4, weight=2, one death.
    # The transition weight depends only on the support intersection.
    supports=[s for s in range(16) if s.bit_count()==2]
    kernel=lambda u,v:F(3,7) if (u&~v).bit_count()==1 else F(0)
    row=sum(kernel(supports[0],v) for v in supports)
    for u in supports:
        assert sum(kernel(u,v) for v in supports)==row
        assert sum(kernel(u,v)*kernel(v,w) for v in supports for w in supports)==row**2
    thirds={(0,0,0):F(1,3),(1,0,0):F(2,3)}
    for power in range(1,5):
        for k in range(power+1):
            exact=F(math.comb(power,k)*2**k,3**power)
            assert coefficient_upper(thirds,power,(k,0,0),precision=12)>=exact
    # Exact characteristic-function check at (3+4i)/5, including laws with
    # nonadjacent support and a degenerate law (zero adjacency overlap).
    for law in ((F(1,2),F(1,2)),(F(1,4),F(1,2),F(1,4)),(F(1,2),F(0),F(1,2)),(F(1),)):
        real=imag=F(0);cr,ci=F(1),F(0)
        for probability in law:
            real+=probability*cr;imag+=probability*ci
            cr,ci=F(3,5)*cr-F(4,5)*ci,F(4,5)*cr+F(3,5)*ci
        adjacency=sum(law[j]*law[j+1] for j in range(len(law)-1))
        assert real*real+imag*imag<=1-F(4,5)*adjacency
    for exponent in range(-12,13):assert log_interval(F(2)**exponent)==(exponent,exponent)
    lo,hi=log_interval(F(3,2));lo2,hi2=log_interval(F(9,4))
    assert lo2<=2*hi and hi2>=2*lo
    for n in (2,4,16,32,80,128):
        assert F(math.comb(n,n//2),1<<n)<=atom_cap(F(n,4))
    # Independent exact check of the THREE-variable partition identity:
    # n=8, two quartets, incoming/next weight 3, one birth and one death.
    weight=lambda a,b,d:F(a+b+1,d+2)
    polynomial={}
    for a in range(5):
        for b in range(5-a):
            for d in range(a+1):
                polynomial[a,b,d]=math.comb(4,a)*math.comb(4-a,b)*math.comb(a,d)*weight(a,b,d)
    coefficient=sum(value*polynomial.get((3-a,1-b,1-d),F(0))
                    for (a,b,d),value in polynomial.items())
    incoming={0,1,2};direct=F(0)
    for block in itertools.combinations(range(8),4):
        block=set(block);other=set(range(8))-block
        for removed in incoming:
            for born in set(range(8))-incoming:
                after=(incoming-{removed})|{born};term=F(1)
                for group in (block,other):
                    a=len(incoming&group);b=len((after-incoming)&group);d=len((incoming-after)&group)
                    term*=weight(a,b,d)
                direct+=term
    assert direct/math.comb(8,4)==coefficient/math.comb(8,3)
    print(dict(integer_dp_tests=26,two_update_row_sum_checks=6,
               characteristic_function_checks=4,log_interval_checks=26,atom_cap_checks=6,
               exact_partition_transitions=1050),flush=True)


def bound(total,probability):
    A=L*W;z=A+(L-1)*D-46*N//100
    pi=math.comb(M,R)*P**R*(1-P)**(M-R)
    fixed=(log_upper(1,math.comb(M,R))+log_upper(1,math.comb(S,W))
           +log_upper(W,1-P)+log_upper(-256,pi)+THETA
           +32*THETA*(A-K+128*R-z)+log_upper(-THETA*z,ZERO))
    kernel=total**128*probability/(math.comb(S,W)*H**W*X**D*Y**D)
    return fixed+log_upper(L-1,kernel)


def quota_bound(records,frozen=False):
    """Uniform atom bound after conditioning each step's (a,b,d) totals.

    For a local T-count law v, s=sum_j v_j v_(j+1) gives
    |phi(t)| <= exp(-2s sin(t/2)^2). Certify an LP dual lower bound on
    sum s over the 128 local groups. Only the rational dual is trusted.
    """
    buckets=defaultdict(lambda:[F(0)]*5)
    for value,a,b,d,t,active in records:buckets[a,b,d][t]+=value
    overlaps={}
    for key,values in buckets.items():
        norm=sum(values);prob=[v/norm for v in values]
        overlaps[key]=sum(prob[j]*prob[j+1] for j in range(4))
    keys=sorted(overlaps)
    if frozen:dual=list(map(F,FROZEN_DUAL))
    else:
        from scipy.optimize import linprog
        rows=np.array([[1,*key] for key in keys],dtype=float)
        fit=linprog(-np.array([128,224,32,32]),A_ub=rows,
                    b_ub=[float(overlaps[k]) for k in keys],bounds=[(None,None)]*4,method='highs')
        assert fit.success
        dual=[F(float(x)).limit_denominator(1<<20) for x in fit.x]
        violation=max(F(0),max(sum(c*x for c,x in zip([1,*key],dual))-overlaps[key] for key in keys))
        # Round the required intercept correction upward on a compact dyadic grid.
        shift=(violation.numerator*(1<<40)+violation.denominator-1)//violation.denominator
        dual[0]-=F(shift+1,1<<40)
    assert all(sum(c*x for c,x in zip([1,*key],dual))<=overlaps[key] for key in keys)
    service=sum(c*x for c,x in zip([128,224,32,32],dual))
    assert service>0
    first,ordinary=atom_cap(15*service),atom_cap(16*service)
    gain=log_upper(1,first)+log_upper(255,ordinary)
    return dict(dual=list(map(str,dual)),service_lower=str(service),service_float=float(service),
                first_cap=str(first),ordinary_cap=str(ordinary),
                quota_log2_upper=str(gain),quota_log2_float=float(gain)),gain


def main():
    global P,ZERO,HOLD,NEXT,DEATH,H,X,Y,THETA
    parser=argparse.ArgumentParser();parser.add_argument('--tests-only',action='store_true')
    parser.add_argument('--tuned',action='store_true')
    parser.add_argument('--frozen',action='store_true')
    args=parser.parse_args()
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        tests()
        if args.tests_only:return
        if args.frozen:
            THETA=F(7,40);P,ZERO,HOLD,NEXT,DEATH=map(F,FROZEN_MARKERS)
            H,X,Y=HOLD*NEXT,NEXT,DEATH/NEXT
        elif args.tuned:
            from scipy.optimize import minimize
            from quartet_trajectories import profile,objective
            THETA=F(7,40)
            f=objective(profile('224',[224],[32]),free_theta=True)
            initial=[math.log2(float(P/(1-P))),math.log2(float(ZERO)),math.log2(float(HOLD)),
                     math.log2(float(NEXT)),math.log2(float(DEATH))]
            fit=minimize(lambda v:f([*v,float(32*THETA)])/K,initial,method='L-BFGS-B',
                         bounds=[(-10,10),(-12,12),(-16,16),(-16,16),(-16,16)],
                         options={'maxiter':120,'ftol':1e-12,'gtol':1e-8})
            P=F(1/(1+2**-float(fit.x[0]))).limit_denominator(1000000)
            ZERO,HOLD,NEXT,DEATH=[F(2**float(v)).limit_denominator(1000000) for v in fit.x[1:]]
            H,X,Y=HOLD*NEXT,NEXT,DEATH/NEXT
            print(dict(theta=str(THETA),p=str(P),zero=str(ZERO),h_old=str(HOLD),
                       next=str(NEXT),death=str(DEATH)),flush=True)
        probs,total,records=local()
        means=[float(sum(value*record[j] for record in records for value in [record[0]])/total)
               for j in range(1,6)]
        print(dict(local_terms=len(probs),tilted_mean_a_b_d_inactive_activecovered=means,
                   diagnostic_only=True),flush=True)
        mass=coefficient_upper(probs)
        upper=bound(total,mass);chernoff=bound(total,F(1))
        print(dict(coefficient_mass_upper=str(mass),directed_log2_upper=str(upper),
                   upper_decimal=float(upper),old_evaluation_upper=float(chernoff),
                   coefficient_gain_bits=float(chernoff-upper),whole_code=False,
                   exact_profile='224 active, 32 deaths per transition, 1974 inactive rows, 46%'),flush=True)
        quota,gain=quota_bound(records,args.frozen)
        if args.frozen:
            assert mass==F(287823,1<<31)
            assert upper+gain==FROZEN_UPPER and FROZEN_UPPER < -104
        print(dict(**quota,combined_log2_upper=str(upper+gain),combined_float=float(upper+gain),
                   directed=True,whole_code=False),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

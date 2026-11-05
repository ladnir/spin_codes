"""Group state histories inside a fixed footprint before the moment.

Floating diagnostics, NOT a certificate. Load uniform_route_moment from
the global helper. Root-owned, file-free, bounded serial Peach runs only.
"""
from math import comb,log,log2,lgamma
from fractions import Fraction as F
from itertools import combinations
import os,resource,shutil,sys,threading,time
import numpy as np

N=1<<21;W=943718;T=4096;LN2=log(2)
BINOM=np.array([(lgamma(513)-lgamma(u+1)-lgamma(513-u))/LN2 for u in range(513)])


def ladd(x,y):
    hi=max(x,y);return hi+log2(1+2**(-abs(x-y)))


def envelope(theta,logz,logx,logh):
    if not (-64<=logz<=0 and -64<=logx<=32 and -64<=logh<=32):return None
    # Sum the two-state chains first using potential(1,1+1/q):
    # their allowed-slot factor is lambda=1+(q+1)*x, x=z^(1/theta).
    logb=theta*ladd(0.,32+log2(1+2**-32)+logz/theta)
    t=logx+logh
    partial=np.logaddexp.reduce(np.array([0.,2+t,log2(6)+2*t,2+3*t])*LN2)/LN2
    local=np.array([ladd(partial,4*t+u*logb)+log2((1,4,6,4,1)[u]) for u in range(5)])*LN2
    poly=np.array([0.])
    for _ in range(128):
        after=np.full(len(poly)+4,-np.inf)
        for u in range(5):after[u:u+len(poly)]=np.logaddexp(after[u:u+len(poly)],poly+local[u])
        poly=after
    ratios=poly/LN2-BINOM-logh*np.arange(513)
    return T*float(ratios.max())+max(0.,512*(logb-logh)),int(np.argmax(ratios)),logb


def scan():
    for fraction in (.125,.25,.4,.45,.5,.625,.75,.8):
        start=time.monotonic();covered=round(fraction*N);best=(float('inf'),None)
        for theta in (1/16,1/8,.2,1/4,1/2,1):
            route=uniform_route_moment(covered,theta)
            def objective(point):
                result=envelope(theta,*point)
                if result is None:return float('inf')
                return theta-32*theta-W*point[0]-covered*point[1]+route+result[0]
            for initial in ([-32*theta,-16*theta,0.],[-32*theta,-32*theta,32*theta]):
                point=initial;value=objective(point)
                for step in (4.,2.,1.,.5,.25,.125):
                    for _ in range(3):
                        changed=False
                        for axis in range(3):
                            for sign in (-1,1):
                                candidate=point.copy();candidate[axis]+=sign*step
                                trial=objective(candidate)
                                if trial<value:point=candidate;value=trial;changed=True
                        if not changed:break
                if value<best[0]:best=(value,(theta,point,envelope(theta,*point)[1:]))
        print(dict(footprint_fraction=fraction,log2_diagnostic=best[0],parameters=best[1],
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def forced_envelope(theta,logz,logx,loghf,loghg):
    if not (-64<=logz<=0 and -64<=logx<=32 and -64<=loghf<=32 and -64<=loghg<=32):return None
    logm=32*theta+theta*log2(1+2**-32)+logz
    logb=theta*ladd(0.,logm/theta)
    t=logx+loghf
    partial=np.logaddexp.reduce(np.array([0.,2+t,log2(6)+2*t,2+3*t])*LN2)/LN2
    rates=[]
    for forced in range(5):
        for optional in range(5-forced):
            full=4*(logx+loghg)+forced*logm+optional*logb
            row=full if forced else ladd(full,partial)
            rates.append((row-forced*loghf-optional*loghg,(forced,optional)))
    maximum=max(rates)
    initial=512*ladd(0.,logx+loghf)
    terminal=512*max(0.,logm-loghf,logb-loghg)
    return initial+128*(T-1)*maximum[0]+terminal,maximum[1]


def forced_scan():
    # Bound the conditional product by its largest quartet row. This is
    # rigorous algebraically, but loses averaging over the input labels.
    for fraction in (.125,.25,.4,.45,.5,.625,.75,.8):
        covered=round(fraction*N);best=(float('inf'),None)
        for theta in (1/16,1/8,.2,1/4,1/2,1):
            route=uniform_route_moment(covered,theta)
            def objective(point):
                result=forced_envelope(theta,*point)
                if result is None:return float('inf')
                return theta-32*theta-W*point[0]-covered*point[1]+route+result[0]
            for initial in ([-32*theta,-16*theta,0.,0.],[0.,-32*theta,16*theta,16*theta]):
                point=initial;value=objective(point)
                for step in (4.,2.,1.,.5,.25,.125,.0625):
                    for _ in range(5):
                        changed=False
                        for axis in range(4):
                            for sign in (-1,1):
                                candidate=point.copy();candidate[axis]+=sign*step
                                trial=objective(candidate)
                                if trial<value:point=candidate;value=trial;changed=True
                        if not changed:break
                if value<best[0]:best=(value,(theta,point,forced_envelope(theta,*point)[1]))
        print(dict(forced_footprint_fraction=fraction,log2_diagnostic=best[0],
                   parameters=best[1],certificate=False),flush=True)


def exact_checks():
    checks=0
    q=16
    for x in (F(1,32),F(1,16),F(1,2),F(1)):
        lam=1+(q+1)*x;v=(F(1),F(q+1,q))
        assert v[0]+q*x*v[1]<=lam*v[0]
        assert x*v[0]+(1+q*x)*v[1]<=lam*v[1]
        for steps in range(1,7):
            for allowed in range(1<<steps):
                row=[F(1),F(0)]
                for t in range(steps):
                    row=[row[0]+x*row[1],
                         (q*x*row[0]+(1+q*x)*row[1]) if (allowed>>t)&1 else F(0)]
                assert sum(row)<=lam**allowed.bit_count()
                checks+=1
    forced_checks=0
    for x in (F(1,32),F(1,16),F(1,2),F(1)):
        c=(q+1)*x;lam=1+c;totals={}
        for first in range(16):
            for second in range(16):
                cover1=first;cover2=15 if first else second
                birth=(second&~first).bit_count();death=(first&~second).bit_count()
                survived=(first&second).bit_count()
                weight=(q*x)**(first.bit_count()+birth)*x**death*(1+q*x)**survived
                key=(cover1,cover2);totals[key]=totals.get(key,F(0))+weight
        for cover1 in range(16):
            for cover2 in range(16):
                allowed1=cover1 if cover2==15 else 0
                forced1=cover1;forced2=cover2 if cover2!=15 else 0
                if forced1&~allowed1:
                    assert totals.get((cover1,cover2),0)==0
                else:
                    allowed=allowed1.bit_count()+cover2.bit_count()
                    forced=forced1.bit_count()+forced2.bit_count()
                    assert totals.get((cover1,cover2),0)<=c**forced*lam**(allowed-forced)
                forced_checks+=1
    kernel_checks=0;b=F(3,2)
    for width in (4,8):
        polynomial={(0,0):F(1)}
        for _ in range(width//4):
            after={}
            for (u,v),value in polynomial.items():
                for i in range(5):
                    for j in range(5):
                        term=comb(4,i)*comb(4,j)*(b**i if j==4 else 1)
                        key=(u+i,v+j);after[key]=after.get(key,F(0))+value*term
            polynomial=after
        partitions=[(15,)] if width==4 else []
        if width==8:
            for others in combinations(range(1,8),3):
                first=1+sum(1<<i for i in others);partitions.append((first,255^first))
        for u in range(width+1):
            incoming=(1<<u)-1;truth=[F(0)]*(width+1)
            for parts in partitions:
                for covered in range(1<<width):
                    allowed=sum((incoming&part).bit_count() for part in parts if part&covered==part)
                    truth[covered.bit_count()]+=b**allowed/len(partitions)
            for v in range(width+1):
                assert polynomial[u,v]/comb(width,u)==truth[v]
                kernel_checks+=1
    print(dict(exact_allowed_mask_chain_checks=checks,exact_footprint_kernel_checks=kernel_checks,
               exact_forced_birth_grouping_checks=forced_checks,
               full_code_certificate=False),flush=True)


def main():
    if shutil.disk_usage('/tmp').free<20*(1<<30):raise RuntimeError('20GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(2*(1<<30),2*(1<<30)))
    timer=threading.Timer(175,lambda:os._exit(124));timer.daemon=True;timer.start()
    if '--checks' in sys.argv:exact_checks()
    elif '--forced' in sys.argv:forced_scan()
    else:scan()
    timer.cancel()


if __name__=='__main__':main()

"""Bounded regional quiet-interval diagnostic, unchanged quartet construction.

FLOAT screens are proposals, not certificates. No files are generated; runs
are serial, guarded by 20 GiB free space and a 60-second watchdog.
"""
import argparse
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np

import quartet_two_update as q
import quartet_history_band as band

ROW_CACHE={}
NEXT_WEIGHT_ROWS={}
P_GRID=(0.24,0.30,0.36,0.42,0.48,0.54,0.60,112627/167924)


def log_choose(n, k):
    if not 0 <= k <= n:
        return -math.inf
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1))/math.log(2)


def normalized_power(values, power):
    """FLOAT normalized polynomial coefficients via log-domain recurrence."""
    degree=len(values)-1
    previous=np.array([0.0])
    for step in range(power):
        current=np.full(degree*(step+1)+1,-np.inf)
        old=np.array([log_choose(degree*step,u) for u in range(len(previous))])
        for a,value in enumerate(values):
            if value<=0:continue
            terms=previous+old+log_choose(degree,a)+math.log2(value)
            current[a:a+len(previous)]=np.logaddexp2(current[a:a+len(previous)],terms)
        current-=np.array([log_choose(degree*(step+1),u) for u in range(len(current))])
        previous=current
    return previous


def local_rows(p, tau=0):
    """All next weights/turnover, fixed incoming weight; full local laws."""
    theta=7/40;ell=float(band.BRIDGE_MARKERS[-1].split('/')[0])/float(band.BRIDGE_MARKERS[-1].split('/')[1])
    X=97905/432928
    sigma=((2**32)*ell)**theta
    death=X**-1*sigma**(tau-1)
    result=[]
    for a in range(5):
        before=(1<<a)-1;value=0.0
        for after in range(16):
            for target in range(16):
                coeff=q.choices(4,before,after,target,True)
                if not any(coeff):continue
                t=target.bit_count();c=t if a else (after&target).bit_count()
                b=(after&~before).bit_count();d=(before&~after).bit_count()
                poly=sum(n*ell**j for j,n in enumerate(coeff))
                value+=p**t*(1-p)**(4-t)*(poly*2**(-32*c))**theta*X**b*death**d
        result.append(value*ell**(-theta*a)*sigma**(tau*a))
    return result


def bounded_death_rows(p):
    """FLOAT 2-D coefficient count; retain deaths, relax next weight only."""
    if p in ROW_CACHE:return ROW_CACHE[p]
    theta=7/40;ell=13592/664659;X=97905/432928
    sigma=((2**32)*ell)**theta
    terms={}
    for a in range(5):
        before=(1<<a)-1
        for after in range(16):
            b=(after&~before).bit_count();d=(before&~after).bit_count()
            for target in range(16):
                coeff=q.choices(4,before,after,target,True)
                if not any(coeff):continue
                t=target.bit_count();c=t if a else (after&target).bit_count()
                poly=sum(n*ell**j for j,n in enumerate(coeff))
                term=p**t*(1-p)**(4-t)*(poly*2**(-32*c))**theta
                term*=X**(b-d)*sigma**-d*ell**(-theta*a)
                terms[a,d]=terms.get((a,d),0.0)+term
    current=np.full((233,41),-np.inf);current[0,0]=0
    for step in range(128):
        out=np.full_like(current,-np.inf)
        for (a,d),value in terms.items():
            nu=min(232-a,4*step)+1;nd=min(40-d,4*step)+1
            out[a:a+nu,d:d+nd]=np.logaddexp2(out[a:a+nu,d:d+nd],
                current[:nu,:nd]+math.log2(value)+log_choose(4,a))
        current=out
    current-=np.array([log_choose(512,u) for u in range(233)])[:,None]
    ROW_CACHE[p]=current
    return current


def retained_next_rows(p):
    """Directed coefficients, FLOAT logs: diagnostic of outgoing relaxation."""
    if p in NEXT_WEIGHT_ROWS:return NEXT_WEIGHT_ROWS[p]
    band.configure()
    q.P,q.H,q.X,q.Y,q.ZERO=map(F,band.BRIDGE_MARKERS)
    q.P=F(p).limit_denominator(1000000);q.ROOTS.clear()
    probability,total,_=q.local()
    box=band.coefficient_box(probability,(232,56,40))
    rows=np.full((17,17),-np.inf)
    for u in range(216,233):
        for d in range(24,41):
            mass=sum(int(box[u,v-u+d,d]) for v in range(216,233))/(2**31)
            rows[u-216,d-24]=(128*math.log2(float(total))+math.log2(mass)-log_choose(512,u)
                -u*math.log2(float(q.H))-(7/40)*u*math.log2(float(q.ZERO))
                -d*math.log2(float(q.X*q.Y))-(7/40)*d*(32+math.log2(float(q.ZERO))))
    NEXT_WEIGHT_ROWS[p]=rows
    return rows


def quota_cap(p, r, nonquiet):
    """Exact identity in log form; FLOAT evaluation only.

    Quiet positions contribute an independent binomial count. Maximize its
    mass over every feasible number of inactive labels in nonquiet positions.
    """
    quiet=8192-512*nonquiet
    lo=max(0,r-512*nonquiet);hi=min(r,quiet)
    mode=max(lo,min(hi,math.floor((quiet+1)*p)))
    pmf=lambda n,k:log_choose(n,k)+k*math.log2(p)+(n-k)*math.log2(1-p)
    return pmf(quiet,mode)-pmf(8192,r)


def block_coefficients(g,a,b,length=16):
    """16-step two-state operator, graded by nonquiet updates.

    State 0 is zero, state 1 is the live band. g,a,b are log2 row bounds.
    Even if g includes resets, the excess is a harmless positive overcount.
    """
    matrix=np.full((2,2,length+1),-np.inf)
    matrix[0,0,0]=matrix[1,1,0]=0
    edges=((0,0,0,0.0),(0,1,1,a),(1,0,1,b),(1,1,1,g))
    for step in range(length):
        out=np.full_like(matrix,-np.inf)
        for i,j,h,value in edges:
            out[:,j,h:h+step+1]=np.logaddexp2(out[:,j,h:h+step+1],matrix[:,i,:step+1]+value)
        matrix=out
    return matrix


def region_power(block, count=256):
    """Keep total nonquiet count; sum all intermediate zero/live states."""
    current=np.full((2,1),-np.inf);current[0,0]=0
    for step in range(count):
        nxt=np.full((2,16*(step+1)+1),-np.inf)
        for i in range(2):
            for j in range(2):
                for h in range(17):
                    if not math.isfinite(block[i,j,h]):continue
                    nxt[j,h:h+current.shape[1]]=np.logaddexp2(
                        nxt[j,h:h+current.shape[1]],current[i]+block[i,j,h])
        current=nxt
    return current


def tests():
    # Exact conditioning check with a nonconstant nonnegative integrand on
    # nonquiet positions. No independence is assumed after conditioning.
    checks=0;p=F(1,3);n=8
    for h in range(n+1):
        for r in range(n+1):
            pi=math.comb(n,r)*p**r*(1-p)**(n-r)
            cap=max((math.comb(n-h,r-s)*p**(r-s)*(1-p)**(n-h-r+s)
                     for s in range(h+1) if 0<=r-s<=n-h),default=F(0))/pi
            unconditioned=F(0);conditioned=F(0)
            for mask in range(1<<h):
                s=mask.bit_count();value=F(2**s+mask%3)
                unconditioned+=value*p**s*(1-p)**(h-s)
                if 0<=r-s<=n-h:
                    conditioned+=value*F(math.comb(n-h,r-s),math.comb(n,r))
            assert conditioned<=cap*unconditioned;checks+=1
    # Exhaustive short paths test nonquiet grading and both boundary states.
    edges=((F(1),F(1,8)),(F(1,32),F(1,4)));path_checks=0
    for length in range(1,6):
        predicted=block_coefficients(-2,-3,-5,length)
        exact={}
        for initial in range(2):
            for path in itertools.product(range(2),repeat=length):
                i=initial;value=F(1);h=0
                for j in path:
                    value*=edges[i][j];h+=bool(i or j);i=j
                key=initial,i,h;exact[key]=exact.get(key,F(0))+value
        for key,value in exact.items():
            assert abs(predicted[key]-math.log2(float(value)))<1e-12;path_checks+=1
    # Interpolation of the two upper bounds, checked without real powers.
    interpolation_checks=0
    for exponent in range(-10,11):
        ratio=F(2)**exponent
        assert min(F(1),ratio)**50<=ratio**3
        interpolation_checks+=1
    print(dict(exact_quota_checks=checks,short_path_coefficients=path_checks,
               exact_interpolation_checks=interpolation_checks),flush=True)


def screen(tau=0, retained_next=False):
    ell=13592/664659;X=97905/432928;theta=7/40
    candidates=[]
    for p in P_GRID:
        # Two expensive coefficient boxes only; other p values retain the
        # cheaper next-weight relaxation. This is a bounded diagnostic.
        retained=(retained_next_rows(p) if retained_next and p>=.6
                  else bounded_death_rows(p)[216:233,24:41])
        score=np.arange(216,233)[:,None]+np.arange(24,41)[None,:]
        g=float(max(np.logaddexp2.reduce(retained+tau*theta*(32+math.log2(ell))*score,axis=1)))
        reset=[1.0]+[sum(math.comb(4,t)*p**t*(1-p)**(4-t)*2**(-32*theta*t)
                         for t in range(a)) for a in range(1,5)]
        reset_rows=normalized_power(reset,128)
        b=max(reset_rows[w]-32*theta*w-w*math.log2(X)-2*theta*w*math.log2(ell)
              +2*tau*theta*(32+math.log2(ell))*w
              for w in range(216,233))
        a=float(np.logaddexp2.reduce([log_choose(512,w)+w*math.log2((1-p)*X)
                                     for w in range(216,233)]))
        candidates.append((p,block_coefficients(g,a,b)))
        if tau==0:print(dict(p=p,full_next_band_row=g,activation=a,reset=b),flush=True)
    curves={}
    for r in (1850,1974,2100):
        by_h=np.full((2,2,17),np.inf);selected={}
        for p,coeff in candidates:
            value=coeff+np.array([quota_cap(p,r,h) for h in range(17)])
            for h in range(17):
                if value[1,1,h]<by_h[1,1,h]:selected[h]=p
            by_h=np.minimum(by_h,value)
        block=np.logaddexp2.reduce(by_h,axis=2)
        # Fixed f on the band, f(0)=1; terminal bound uses its band minimum.
        terminal=np.array([0.0,-232*(math.log2(X)+theta*math.log2(ell)-tau*theta*(32+math.log2(ell)))])
        for _ in range(256):terminal=np.logaddexp2.reduce(block+terminal[None,:],axis=1)
        base=(log_choose(8192,r)+theta+32*theta*(943718-1048576+128*r)
              +theta*943718*math.log2(ell)-tau*theta*943718*(32+math.log2(ell)))
        print(dict(inactive_rows=r,tau=tau,block_log_entries=block.tolist(),log2_union_screen=base+terminal[0],
                   directed=False,whole_code=False),flush=True)
        counted=region_power(by_h)
        counted[1]-=232*(math.log2(X)+theta*math.log2(ell)-tau*theta*(32+math.log2(ell)))
        by_count=np.logaddexp2.reduce(counted,axis=0)
        # For the actual retained histories (not their enlarged envelopes),
        # A+D <= 272 h. The zero-selection-free count is also available.
        correction=np.minimum(0,(1-tau)*theta*(32+math.log2(ell))*(272*np.arange(4097)-943718))
        refined=base+by_count+correction
        curves[r]=refined
        where=int(np.argmax(refined))
        print(dict(inactive_rows=r,tau=tau,zero_free_combined_screen=float(np.logaddexp2.reduce(refined)),
                   worst_nonquiet=where,worst_count_screen=float(refined[where]),
                   directed=False,whole_code=False),flush=True)
    return curves


def main():
    global P_GRID
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser();parser.add_argument('--retain-next',action='store_true')
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--extended-p',action='store_true')
        args=parser.parse_args()
        tests()
        if args.tests_only:return
        if args.extended_p:P_GRID=P_GRID+(.70,.73,.76)
        best={}
        for tau in ((0,.03,.06,.1,.15) if args.extended_p else
                    (0,.03,.06,.1) if args.retain_next else (0,.01,.03,.06,.1,.15,.2)):
            curves=screen(tau,args.retain_next)
            for r,curve in curves.items():
                best[r]=np.minimum(best.get(r,np.full(4097,np.inf)),curve)
        for r,curve in best.items():
            print(dict(inactive_rows=r,combined_tau_per_count=float(np.logaddexp2.reduce(curve)),
                       worst_nonquiet=int(np.argmax(curve)),directed=False,whole_code=False),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

"""FLOAT diagnostic: keep inactive occupancy through a complete region.

Unchanged quartet construction. A coefficientwise live-state envelope sums
all next weights/turnover; hence it bounds, but does not exactly count, the
old restricted band. No numerical certificate is claimed. Serial, <16 MiB
arrays, 60-second watchdog, no output files, and 20 GiB disk guard.
"""
import argparse
from fractions import Fraction as F
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np

import quartet_two_update as q
import quartet_quiet_regions as quiet


def log_convolve(left,right,cap):
    """Positive polynomial product in log space; no FFT or tail cutoff."""
    if len(left)<len(right):left,right=right,left
    out=np.full(min(cap+1,len(left)+len(right)-1),-np.inf)
    for j,value in enumerate(right):
        if not math.isfinite(value) or j>=len(out):continue
        count=min(len(left),len(out)-j)
        out[j:j+count]=np.logaddexp2(out[j:j+count],left[:count]+value)
    return out


def occupancy_power(terms,power=128,umax=232):
    """Coefficient table in incoming weight and inactive-label count."""
    current=np.full((min(umax,4*power)+1,4*power+1),-np.inf)
    current[0,0]=0
    for step in range(power):
        nxt=np.full_like(current,-np.inf)
        for (a,t),value in terms.items():
            nu=min(umax-a,4*step)+1;nt=4*step+1
            if nu<=0:continue
            nxt[a:a+nu,t:t+nt]=np.logaddexp2(nxt[a:a+nu,t:t+nt],
                current[:nu,:nt]+math.log2(value)+quiet.log_choose(4,a))
        current=nxt
    current-=np.array([quiet.log_choose(4*power,u) for u in range(current.shape[0])])[:,None]
    return current


def local_terms(p,tau,reset=False,X=97905/432928,theta=7/40,ell=13592/664659,death_marker=1):
    sigma=((2**32)*ell)**theta
    terms={}
    for a in range(5):
        before=(1<<a)-1
        for after in ((0,) if reset else range(16)):
            b=(after&~before).bit_count();d=(before&~after).bit_count()
            for target in range(16):
                coeff=q.choices(4,before,after,target,True)
                if not any(coeff):continue
                t=target.bit_count();c=t if a else (after&target).bit_count()
                poly=sum(n*ell**j for j,n in enumerate(coeff))
                value=p**t*(1-p)**(4-t)*(poly*2**(-32*c))**theta
                value*=X**(b-d)*sigma**((tau-1)*d+tau*a)*ell**(-theta*a)*death_marker**d
                terms[a,t]=terms.get((a,t),0.0)+value
    return terms


def region_matrix(edges,length=16,cap=2100):
    current=[[np.array([0.0]),np.array([-np.inf])],
             [np.array([-np.inf]),np.array([0.0])]]
    for _ in range(length):
        nxt=[[None,None],[None,None]]
        for i in range(2):
            for j in range(2):
                terms=[log_convolve(current[i][k],edges[k][j],cap) for k in range(2)]
                size=max(map(len,terms));value=np.full(size,-np.inf)
                for term in terms:value[:len(term)]=np.logaddexp2(value[:len(term)],term)
                nxt[i][j]=value
        current=nxt
    return current


def tests():
    checks=0
    for left in ((F(1),F(2)),(F(1,8),F(0),F(3,7))):
        for right in ((F(2),F(1,3)),(F(0),F(1))):
            for cap in range(5):
                log=lambda seq:np.array([math.log2(float(x)) if x else -np.inf for x in seq])
                result=log_convolve(log(left),log(right),cap)
                for i,v in enumerate(result):
                    exact=sum((left[j]*right[i-j] for j in range(len(left)) if 0<=i-j<len(right)),F(0))
                    assert (not exact and v==-np.inf) or abs(v-math.log2(float(exact)))<1e-12
                    checks+=1
    # The scalar graded helper and the polynomial matrix code must agree.
    edges=[[np.array([0.,-np.inf]),np.array([-np.inf,-3.])],
           [np.array([-np.inf,-5.]),np.array([-np.inf,-2.])]]
    for length in range(1,5):
        expected=quiet.block_coefficients(-2,-3,-5,length)
        actual=region_matrix(edges,length,length)
        for i in range(2):
            for j in range(2):
                assert np.allclose(actual[i][j],expected[i,j],atol=1e-12)
    terms={(0,0):F(1,2),(1,2):F(1,3),(3,1):F(1,5),(4,4):F(1,7)}
    normalized=occupancy_power(terms,power=2,umax=8)
    exact={}
    for (a,t),x in terms.items():
        for (b,s),y in terms.items():
            key=a+b,t+s
            exact[key]=exact.get(key,F(0))+math.comb(4,a)*math.comb(4,b)*x*y
    for (a,t),value in exact.items():
        assert abs(normalized[a,t]-math.log2(float(value/F(math.comb(8,a)))))<1e-12
    print(dict(polynomial_coefficients=checks,matrix_comparisons=16,
               normalized_occupancy_coefficients=len(exact)),flush=True)


def run(tau,p,X,theta=7/40,ell=13592/664659,next_tilts=False):
    stamp=time.monotonic()
    live=occupancy_power(local_terms(p,tau,X=X,theta=theta,ell=ell))
    if next_tilts:
        for z,y in ((.75,1),(1.25,1),(1.5,1),(2,1),(1,.75),(1.25,.75),(1.5,.75),(2,.75)):
            alternative=occupancy_power(local_terms(p,tau,X=X*z,theta=theta,ell=ell,death_marker=y))
            # Bound v in [216,232], d in [24,40] coefficientwise in t.
            # X*z inserts z^(v-u); an independent marker inserts y^d.
            v_edge=216 if z>=1 else 232;d_edge=24 if y>=1 else 40
            alternative+=(np.arange(233)-v_edge)[:,None]*math.log2(z)-d_edge*math.log2(y)
            live=np.minimum(live,alternative)
            print(dict(phase='outgoing tilt',z=z,y=y,seconds=time.monotonic()-stamp),flush=True)
    g=live[216:233].max(axis=0)
    print(dict(phase='live coefficients',seconds=time.monotonic()-stamp),flush=True)
    reset=occupancy_power(local_terms(p,tau,reset=True,X=X,theta=theta,ell=ell))
    b=reset[216:233].max(axis=0)
    print(dict(phase='reset coefficients',seconds=time.monotonic()-stamp),flush=True)
    c=np.array([quiet.log_choose(512,t)+t*math.log2(p)+(512-t)*math.log2(1-p) for t in range(513)])
    a=np.full(513,-np.inf)
    for t in range(297):
        a[t]=c[t]+np.logaddexp2.reduce([quiet.log_choose(512-t,w)+w*math.log2(X)
                                       for w in range(216,min(232,512-t)+1)])
    region=region_matrix([[c,a],[b,g]])
    print(dict(phase='joint region',seconds=time.monotonic()-stamp),flush=True)
    result={}
    for r in (1850,1974,2100):
        pi=quiet.log_choose(8192,r)+r*math.log2(p)+(8192-r)*math.log2(1-p)
        matrix=np.array([[region[i][j][r]-pi for j in range(2)] for i in range(2)])
        terminal=np.array([0.,-232*(math.log2(X)+theta*math.log2(ell)-tau*theta*(32+math.log2(ell)))])
        for _ in range(256):terminal=np.logaddexp2.reduce(matrix+terminal[None,:],axis=1)
        base=(quiet.log_choose(8192,r)+theta+32*theta*(943718-1048576+128*r)
              +theta*943718*math.log2(ell)-tau*theta*943718*(32+math.log2(ell)))
        print(dict(inactive_rows=r,tau=tau,p=p,X=X,theta=theta,ell=ell,region_entries=matrix.tolist(),
                   log2_screen=float(base+terminal[0]),directed=False,whole_code=False,
                   live_next_weights_and_deaths_relaxed=not next_tilts,
                   outgoing_constraints_chernoff=next_tilts),flush=True)
        result[r]=float(base+terminal[0])
    return result


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:
        parser=argparse.ArgumentParser();parser.add_argument('--tau',type=float,default=.1)
        parser.add_argument('--p',type=float,default=112627/167924)
        parser.add_argument('--x',type=float,default=97905/432928)
        parser.add_argument('--theta',type=float,default=7/40)
        parser.add_argument('--ell',type=float,default=13592/664659)
        parser.add_argument('--next-tilts',action='store_true')
        parser.add_argument('--compare-p',action='store_true')
        parser.add_argument('--tests-only',action='store_true');args=parser.parse_args()
        assert 0<=args.tau<=1 and 0<args.p<1 and args.x>0 and 0<args.theta<=1
        assert args.ell>0 and (1<<32)*args.ell>=1
        assert not (args.compare_p and args.next_tilts)  # Combined run must stay bounded.
        tests()
        if not args.tests_only:
            result=run(args.tau,args.p,args.x,args.theta,args.ell,args.next_tilts)
            if args.compare_p:
                alternate=run(args.tau,.5,args.x,args.theta,args.ell)
                error=max(abs(result[r]-alternate[r]) for r in result)
                assert error<1e-6
                print(dict(artificial_tilt_invariance_max_log2_error=error),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

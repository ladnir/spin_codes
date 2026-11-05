"""FLOAT state-weight flow diagnostic for the unchanged quartet code.

Retain every weight 216..232 plus zero, and exact permitted deaths 24..40.
Compare with the scalar collapse using the SAME transition coefficients.
Not a directed certificate, not whole-code coverage. Serial, <16 MiB main
coefficient arrays, 60-second watchdog, no files, >=20 GiB free.
"""
import argparse
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np

import quartet_two_update as q
import quartet_history_band as band
import quartet_quiet_regions as quiet

THETA=7/40
P0,H,X,Y,ELL=(float(q.F(x)) for x in band.BRIDGE_MARKERS)
LOGSIGMA=THETA*(32+math.log2(ELL))
WEIGHTS=np.arange(216,233)


def floating_box(p):
    terms={};active_moment={};inactive={}
    for a in range(5):
        before=(1<<a)-1
        for after in range(16):
            b=(after&~before).bit_count();d=(before&~after).bit_count()
            for target in range(16):
                coefficients=q.choices(4,before,after,target,True)
                if not any(coefficients):continue
                t=target.bit_count();covered=t if a else (target&after).bit_count()
                poly=sum(n*ELL**j for j,n in enumerate(coefficients))
                value=math.comb(4,a)*p**t*(1-p)**(4-t)*(poly*2**(-32*covered))**THETA
                value*=H**a*X**b*Y**d
                key=a,b,d;terms[key]=terms.get(key,0)+value
                active_cover=4-t if a else after.bit_count()
                active_moment[key]=active_moment.get(key,0)+active_cover*value
                inactive[key]=inactive.get(key,0)+t*value
    total=sum(terms.values());terms={key:value/total for key,value in terms.items()}
    active_moment={key:value/total for key,value in active_moment.items()}
    inactive={key:value/total for key,value in inactive.items()}
    shape=(233,57,41);assert math.prod(shape)*8*3<16*(1<<20)
    current=np.zeros(shape);current[0,0,0]=1.;nxt=np.zeros_like(current)
    for step in range(128):
        nxt.fill(0)
        for (a,b,d),value in terms.items():
            extent=tuple(min(n-k,4*step)+1 for k,n in zip((a,b,d),(232,56,40)))
            source=tuple(slice(0,e) for e in extent)
            dest=tuple(slice(k,k+e) for k,e in zip((a,b,d),extent))
            nxt[dest]+=current[source]*value
        current,nxt=nxt,current
    # nxt still contains the 127-factor polynomial. Differentiate the
    # positive 128th power by marking one quartet; no second large DP.
    for u,v,d in ((216,216,31),(224,224,32),(232,232,40)):
        target=(u,v-u+d,d);moments=[]
        for marked in (active_moment,inactive):
            numerator=0.
            for degrees,value in marked.items():
                index=tuple(n-k for n,k in zip(target,degrees))
                if min(index)>=0:numerator+=value*nxt[index]
            moments.append(128*numerator/current[target])
        print(dict(p=p,transition=(u,v,d),tilted_active_coverage=moments[0],
                   tilted_inactive_labels=moments[1],conditional_on_transition=True,
                   regional_quota_not_conditioned=True),flush=True)
    logs=np.empty((17,17,17))
    for i,u in enumerate(WEIGHTS):
        for j,v in enumerate(WEIGHTS):
            for k,d in enumerate(range(24,41)):
                mass=current[u,v-u+d,d]
                assert mass>0 and math.isfinite(mass)  # Never certify an underflow.
                logs[i,j,k]=(128*math.log2(total)+math.log2(mass)-quiet.log_choose(512,int(u))
                    -u*math.log2(H)-THETA*u*math.log2(ELL)-d*math.log2(X*Y)-LOGSIGMA*d)
    reset=[1.]+[sum(math.comb(4,t)*p**t*(1-p)**(4-t)*2**(-32*THETA*t)
                    for t in range(a)) for a in range(1,5)]
    r=quiet.normalized_power(reset,128)
    resets=np.array([r[w]-32*THETA*w-w*math.log2(X)-2*THETA*w*math.log2(ELL) for w in WEIGHTS])
    activation=np.array([quiet.log_choose(512,int(w))+w*math.log2((1-p)*X) for w in WEIGHTS])
    return logs,resets,activation


def step_matrix(data,tau):
    logs,resets,activation=data
    result=np.full((18,18),-np.inf);result[0,0]=0
    scores=WEIGHTS[:,None,None]+np.arange(24,41)[None,None,:]
    result[1:,1:]=np.logaddexp2.reduce(logs+tau*LOGSIGMA*scores,axis=2)
    result[0,1:]=activation;result[1:,0]=resets+2*tau*LOGSIGMA*WEIGHTS
    return result


def region_coefficients(matrix,length=16):
    """Matrix polynomial power; only the zero-to-zero edge is quiet."""
    n=len(matrix);nonquiet=matrix.copy();nonquiet[0,0]=-np.inf
    current=np.full((1,n,n),-np.inf);np.fill_diagonal(current[0],0)
    for step in range(length):
        nxt=np.full((step+2,n,n),-np.inf)
        # Shapes: degree, initial, middle, next; sum the shared middle state.
        nxt[1:]=np.logaddexp2.reduce(current[:,:,:,None]+nonquiet[None,None,:,:],axis=2)
        nxt[:-1,:,0]=np.logaddexp2(nxt[:-1,:,0],current[:,:,0]+matrix[0,0])
        current=nxt
    return current


def complete(regions,r,tau,weights):
    matrix=np.logaddexp2.reduce(regions,axis=0)
    terminal=np.array([0.]+[-w*(math.log2(X)+THETA*math.log2(ELL)-tau*LOGSIGMA) for w in weights])
    for _ in range(256):terminal=np.logaddexp2.reduce(matrix+terminal[None,:],axis=1)
    base=(quiet.log_choose(8192,r)+THETA+32*THETA*(943718-1048576+128*r)
          +THETA*943718*math.log2(ELL)-tau*LOGSIGMA*943718)
    return float(base+terminal[0])


def bounds(data,tau,with_scalar=False):
    full={r:np.full((17,18,18),np.inf) for r in (1850,1974,2100)}
    scalar={r:np.full((17,2,2),np.inf) for r in full}
    for p,record in data.items():
        matrix=step_matrix(record,tau);coeff=region_coefficients(matrix)
        if with_scalar:
            small=np.array([[0.,np.logaddexp2.reduce(matrix[0,1:])],
                            [max(matrix[1:,0]),max(np.logaddexp2.reduce(matrix[1:,1:],axis=1))]])
            scoeff=region_coefficients(small)
        for r in full:
            cap=np.array([quiet.quota_cap(p,r,h) for h in range(17)])[:,None,None]
            full[r]=np.minimum(full[r],coeff+cap)
            if with_scalar:scalar[r]=np.minimum(scalar[r],scoeff+cap)
    result={}
    for r in full:
        f=complete(full[r],r,tau,WEIGHTS)
        s=complete(scalar[r],r,tau,[232]) if with_scalar else None
        if with_scalar:assert f<=s+1e-7
        result[r]=(f,s)
    return result


def tests():
    matrix=np.array([[0.,-3.],[-5.,-2.]])
    for length in range(1,6):
        actual=region_coefficients(matrix,length)
        expected=quiet.block_coefficients(-2,-3,-5,length).transpose(2,0,1)
        assert np.allclose(actual,expected,atol=1e-12)
    # Endpoint aggregation on an arbitrary positive multi-state kernel must
    # dominate every iterate when the band terminal vector is bounded.
    raw=np.array([[1.,.03,.02],[.01,.1,.02],[.02,.01,.08]])
    scalar=np.array([[raw[0,0],sum(raw[0,1:])],[max(raw[1:,0]),max(raw[1:,1:].sum(axis=1))]])
    vector=np.ones(3);envelope=np.ones(2)
    for _ in range(12):
        vector=raw@vector;envelope=scalar@envelope
        assert vector[0]<=envelope[0]+1e-15 and max(vector[1:])<=envelope[1]+1e-15
    print(dict(graded_matrix_tests=5,scalar_domination_steps=12),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser();parser.add_argument('--coverage-only',action='store_true')
        parser.add_argument('--optimize-tau',action='store_true')
        args=parser.parse_args();tests();data={}
        if args.coverage_only:
            floating_box(P0)
            print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
            return
        # Matched compact grid for the first causal comparison; the full
        # eight-point collection exceeded the per-run watchdog.
        for p in (.42,.54,.60,P0):
            data[p]=floating_box(p)
            print(dict(phase='transition coefficients',p=p,seconds=time.monotonic()-start),flush=True)
        for tau in (0,.03,.06,.1):
            for r,(f,s) in bounds(data,tau,True).items():
                print(dict(inactive_rows=r,tau=tau,state_flow=f,scalar_same_coefficients=s,
                           gain_bits=s-f,directed=False,whole_code=False),flush=True)
        if args.optimize_tau:
            from scipy.optimize import minimize_scalar
            result=minimize_scalar(lambda t:max(v[0] for v in bounds(data,t).values()),
                                   bounds=(.06,.12),method='bounded',
                                   options=dict(maxiter=20,xatol=1e-6))
            chosen=round(float(result.x),6)
            print(dict(proposed_tau=chosen,proposal_success=bool(result.success),
                       evaluations=result.nfev,finite_certificate=False),flush=True)
            for r,(f,s) in bounds(data,chosen,True).items():
                print(dict(inactive_rows=r,proposed_tau=chosen,state_flow=f,
                           scalar_same_coefficients=s,gain_bits=s-f,
                           directed=False,whole_code=False),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

"""FLOAT all-population screen with independent seven-type potentials.

No construction changes; all destination types are summed. This tests whether
balancing the three selected-age ratios caused the global envelope to fail.
No generated files. Run serially; 180-second deadline, <256 MiB arrays.
"""
from fractions import Fraction as F
import math
import os
from pathlib import Path
import shutil
import sys
import threading
import time
import numpy as np
from quartet_two_cohort import compositions
from quartet_closed_population import multinomial
from quartet_closed_factored import quartet_value

LOCAL = np.array(list(compositions(4, 7)), dtype=np.int16)
BALANCED = (73,)*6+(74,)
MEAN_ATOM = F(multinomial(BALANCED))*math.prod(F(n,512)**n for n in BALANCED)
ATOM_COST = -math.log2(MEAN_ATOM)
ROW_COST = math.log2(math.comb(8192,2100))


def coefficients(v):
    h = np.exp2(np.r_[0., v[:6]])
    x = 2**v[6]; a = 2**(707/1024); d = 2**(-1905/1024)
    paid_marker = 2**v[7] if len(v)>7 else 1.
    t = 1+x*paid_marker/8; base = h[1]+a*t*h[4]
    nonfull = sum(math.comb(4,j)*x**j*base**(4-j) for j in range(1,5))
    live = np.array([h[1]+a*h[4], h[2], h[3], 1., a*h[5], a*h[6], d*h[3]+a*h[6]])/h
    paid = np.array([base,h[2]/h[1],h[3]/h[2]])
    selected = LOCAL[:,4:].sum(axis=1)>0; blocked = LOCAL[:,3]>0
    values = np.zeros(len(LOCAL)); keep = selected & ~blocked
    values[keep] = t**4*np.prod(live**LOCAL[keep],axis=1)
    keep = ~selected
    values[keep] = nonfull/np.prod(h**LOCAL[keep],axis=1)
    keep = ~selected & ~blocked
    values[keep] += np.prod(paid**LOCAL[keep,:3],axis=1)
    return values


def elevation(degree):
    pops = np.array(list(compositions(degree,7)),dtype=np.int16)
    weights = np.empty((len(pops),len(LOCAL)))
    assert weights.nbytes < 128*(1<<20)
    choose = np.zeros((degree+1,5))
    for n in range(degree+1):
        for k in range(min(n,4)+1): choose[n,k]=math.comb(n,k)
    for j,m in enumerate(LOCAL):
        weights[:,j]=np.prod([choose[pops[:,i],m[i]] for i in range(7)],axis=0)/math.comb(degree,4)
    assert max(abs(weights.sum(axis=1)-1))<1e-12
    return pops,weights


def objective(v,weights,temp=0.):
    vals=weights@coefficients(v); logs=np.log2(vals[vals>0]); peak=max(logs)
    bound=peak if not temp else peak+temp*math.log2(sum(np.exp2((logs-peak)/temp)))
    cap = -255463*v[7]/(1<<21) if len(v)>7 else 0.
    return .25*bound-2100/8192*v[6]-79174/(1<<21)+ATOM_COST/512-255*ROW_COST/(1<<21)+max(0.,-min(v[3:6]))/4096+cap


def gradient(v,w,t):
    grad=np.zeros(len(v)); eps=2e-5
    for i in range(len(v)):
        a=v.copy(); b=v.copy(); a[i]+=eps; b[i]-=eps
        grad[i]=(objective(a,w,t)-objective(b,w,t))/(2*eps)
    return objective(v,w,t),grad


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124)); timer.daemon=True; timer.start()
    start=time.monotonic()
    try:
        v=np.array([-.044,-.031,-.055,-2.,-1.5,-1.,-.123])
        cap = '--paid-cap' in sys.argv
        if cap:v=np.r_[v,-1.]
        for trial in (v,np.array([.1,-.2,.3,-.4,.5,-.6,.7])):
            h=np.exp2(np.r_[0.,trial[:6]])
            y=(2**trial[7] if len(trial)>7 else 1.)/8
            reference=np.array([quartet_value(m,h,2**(707/1024),2**(-1905/1024),2**trial[6],y)/math.prod(h**m) for m in LOCAL])
            assert max(abs(coefficients(trial)-reference))<1e-9
        print(dict(factored_coefficient_checks=2*len(LOCAL),atom_cost=ATOM_COST),flush=True)
        for degree,temp in ((6,.03),(8,.005),(8,.0005)):
            pops,weights=elevation(degree); inv=np.eye(len(v)); value,grad=gradient(v,weights,temp)
            for iteration in range(120):
                direction=-inv@grad
                if grad@direction>=0: direction=-grad;inv=np.eye(len(v))
                direction/=max(1.,max(abs(direction))/2);rate=1.
                for _ in range(25):
                    candidate=np.clip(v+rate*direction,-16,16);step=candidate-v
                    if cap:candidate[7]=min(0.,candidate[7]);step=candidate-v
                    if objective(candidate,weights,temp)<=value+1e-4*grad@step:break
                    rate*=.5
                if max(abs(step))<1e-9:break
                nv,ng=gradient(candidate,weights,temp);delta=ng-grad;curvature=delta@step
                if curvature>1e-12:
                    correction=np.eye(len(v))-np.outer(step,delta)/curvature
                    inv=correction@inv@correction.T+np.outer(step,step)/curvature
                else: inv=np.eye(len(v))
                v,value,grad=candidate,nv,ng
            vals=weights@coefficients(v);idx=int(np.argmax(vals))
            print(dict(degree=degree,temperature=temp,parameters=v.tolist(),worst=pops[idx].tolist(),
                       float_global_log2_upper=float(objective(v,weights)*(1<<21)),seconds=time.monotonic()-start),flush=True)
        pops,weights=elevation(16);vals=weights@coefficients(v);idx=int(np.argmax(vals))
        print(dict(final_degree=16,worst=pops[idx].tolist(),float_global_log2_upper=float(objective(v,weights)*(1<<21)),
                   without_mean_atom_penalty=float(objective(v,weights)*(1<<21)-4096*ATOM_COST),
                   directed=False,seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

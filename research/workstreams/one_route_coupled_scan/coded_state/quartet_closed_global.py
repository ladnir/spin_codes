"""All-population product-potential bound using the exact factorization.

The five source cohorts retain Z,V0,V1,B and total selected count. Selected
age weights are balanced; no histories or destinations are omitted.
FLOAT optimizer proposes a potential, never certifies it. Serial Peach180s,
<=256MiB, no generated files.
"""
from fractions import Fraction as F
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
from quartet_two_cohort import compositions
from quartet_closed_population import multinomial


LOCAL=list(compositions(4,5))
P=F(2100,8192)
BALANCED=(102,102,102,103,103)
MEAN_ATOM=F(multinomial(BALANCED))*math.prod(F(n,512)**n for n in BALANCED)
ATOM_COST=-math.log2(MEAN_ATOM)
ROW_COST=math.log2(math.comb(8192,2100))


def parameters(v):
    u,w,blocked,gap,x=np.exp2(v)
    activation=2**(707/1024);death=2**(-1905/1024);rho=activation+gap
    l2=death*blocked/gap;l1=activation*l2/rho;l0=activation*l1/rho
    return (1.,u,w,blocked,l0,l1,l2),float(x),float(rho)


def coefficients(v):
    h,x,rho=parameters(v);a=2**(707/1024);t=1+x/8
    base=h[1]+a*t*h[4]
    nonfull=sum(math.comb(4,j)*x**j*base**(4-j) for j in range(1,5))
    ordinary=(h[1]+a*h[4],h[2]/h[1],h[3]/h[2])
    paid=(base,h[2]/h[1],h[3]/h[2]);values=[]
    for z,v0,v1,b,s in LOCAL:
        if s and b:values.append(0.);continue
        if s:
            value=t**4*ordinary[0]**z*ordinary[1]**v0*ordinary[2]**v1*rho**s
        else:
            value=nonfull/(h[1]**v0*h[2]**v1*h[3]**b)
            if not b:value+=paid[0]**z*paid[1]**v0*paid[2]**v1
        values.append(value)
    return np.array(values)


def elevation(degree):
    pops=np.array(list(compositions(degree,5)),dtype=np.int16)
    choose=np.zeros((degree+1,5))
    for n in range(degree+1):
        for k in range(min(4,n)+1):choose[n,k]=math.comb(n,k)
    denominator=math.comb(degree,4)
    weights=np.empty((len(pops),len(LOCAL)))
    for j,m in enumerate(LOCAL):weights[:,j]=np.prod([choose[pops[:,i],m[i]] for i in range(5)],axis=0)/denominator
    assert max(abs(weights.sum(axis=1)-1))<1e-12
    return pops,weights


def objective(v,weights,temperature=0.):
    values=weights@coefficients(v);logs=np.log2(values);peak=max(logs)
    bound=peak if not temperature else peak+temperature*math.log2(sum(np.exp2((logs-peak)/temperature)))
    h,x,rho=parameters(v);tail=max(0.,-math.log2(min(h[4:])))
    return .25*bound-float(P)*v[4]-79174/(1<<21)+ATOM_COST/512-255*ROW_COST/(1<<21)+tail/4096


def gradient(v,weights,temperature):
    value=objective(v,weights,temperature);grad=np.zeros(5);epsilon=2e-5
    for i in range(5):
        plus=v.copy();minus=v.copy();plus[i]+=epsilon;minus[i]-=epsilon
        grad[i]=(objective(plus,weights,temperature)-objective(minus,weights,temperature))/(2*epsilon)
    return value,grad


def optimize():
    start=time.monotonic();v=np.array([0.,0.,0.,-1.4,0.])
    lo=np.array([-16.,-16.,-16.,-16.,-8.]);hi=np.array([16.,16.,16.,10.,8.])
    for degree,temp in ((8,.05),(16,.01),(16,.002)):
        pops,weights=elevation(degree);inverse=np.eye(5);value,grad=gradient(v,weights,temp)
        for iteration in range(100):
            direction=-inverse@grad
            if grad@direction>=0:direction=-grad;inverse=np.eye(5)
            direction/=max(1.,max(abs(direction))/2);rate=1.
            for _ in range(25):
                candidate=np.clip(v+rate*direction,lo,hi);step=candidate-v
                nxt=objective(candidate,weights,temp)
                if nxt<=value+1e-4*(grad@step):break
                rate*=.5
            if max(abs(step))<1e-9:break
            nv,ng=gradient(candidate,weights,temp);delta=ng-grad;curvature=delta@step
            if curvature>1e-12:
                correction=np.eye(5)-np.outer(step,delta)/curvature
                inverse=correction@inverse@correction.T+np.outer(step,step)/curvature
            else:inverse=np.eye(5)
            v,value,grad=candidate,nv,ng
        print(dict(degree=degree,temperature=temp,parameters=list(map(float,v)),
                   float_global_log2_upper=float(objective(v,weights)*(1<<21)),
                   seconds=time.monotonic()-start),flush=True)
    pops,weights=elevation(32);k=weights@coefficients(v);index=int(np.argmax(k))
    print(dict(final_parameters=list(map(float,v)),source_degree32_max=pops[index].tolist(),
               log2_quartic_upper=float(math.log2(k[index])),atom_cost=ATOM_COST,
               float_global_log2_upper=float(objective(v,weights)*(1<<21)),
               all_source_populations=True,all_destinations=True,
               all_original_quota_allocations=True,directed=False,
               seconds=time.monotonic()-start),flush=True)
    return v


def tests():
    # The multinomial point mass at its mean is minimized by balanced
    # counts. Exhaustive small checks supplement the convexity argument.
    checked=0
    for width in range(1,13):
        for categories in range(2,6):
            low=width//categories;balanced=(low,)*(categories-width%categories)+(low+1,)*(width%categories)
            expected=F(multinomial(balanced))*math.prod(F(n,width)**n for n in balanced if n)
            for pop in compositions(width,categories):
                atom=F(multinomial(pop))*math.prod(F(n,width)**n for n in pop if n)
                assert atom>=expected;checked+=1
    # Exact integer degree-elevation weights sum to binom(D,4).
    elevations=0
    for degree in (4,5,8):
        for pop in compositions(degree,5):
            count=sum(math.prod(math.comb(n,m) if n>=m else 0 for n,m in zip(pop,local)) for local in LOCAL)
            assert count==math.comb(degree,4);elevations+=1
    print(dict(exact_mean_atom_comparisons=checked,exact_elevation_rows=elevations),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:tests();optimize()
    finally:timer.cancel()


if __name__=='__main__':main()

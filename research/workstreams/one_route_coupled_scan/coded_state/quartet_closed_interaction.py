"""FLOAT diagnostic for a blocked/selected population interaction potential.

Phi(n)=prod h_i^n_i * 2^(-gamma*n_B*n_S/512). No histories are deleted.
Pure-source screens guide proposals; mixed-source checks sum all outgoing
populations but are NOT global coverage. Serial, 180s, <256MiB, no files.
"""
import math
import os
from pathlib import Path
import shutil
import threading
import time
from collections import defaultdict
import numpy as np
from quartet_closed_quotient import kernel
from quartet_closed_population import multinomial
from quartet_closed_factored import source_fit
from quartet_two_cohort import compositions

WIDTH=512
COUNTS=np.arange(WIDTH+1)
CHOOSE=np.array([math.log2(math.comb(WIDTH,int(j))) for j in COUNTS])
ROW_COST=math.log2(math.comb(8192,2100))


def logsum(values):
    peak=max(values);return peak+math.log2(sum(np.exp2(values-peak)))


def pure_rows(v,gamma):
    h=np.exp2(np.r_[0.,v[:6]]);x=2**v[6];a=2**(707/1024);d=2**(-1905/1024);t=1+x/8
    base=h[1]+a*t*h[4];n=sum(math.comb(4,j)*x**j*base**(4-j) for j in range(1,5))
    rows=[math.log2(x+base),math.log2(n+h[2]**4)/4-v[0],
          math.log2(n+h[3]**4)/4-v[1],math.log2(n)/4-v[2],
          math.log2(a*t)+v[4]-v[3],math.log2(a*t)+v[5]-v[4]]
    terms=CHOOSE+COUNTS*math.log2(d*h[3])+(WIDTH-COUNTS)*math.log2(a*h[6])-gamma*COUNTS*(WIDTH-COUNTS)/WIDTH
    rows.append(math.log2(t)-v[5]+logsum(terms)/WIDTH)
    # V1 ignores its outgoing interaction penalty, so that particular
    # row is an upper envelope. Other homogeneous rows are exact.
    return np.array(rows)


def objective(v,gamma,temp=0.):
    rows=pure_rows(v,gamma);peak=max(rows)
    value=peak if not temp else peak+temp*math.log2(sum(np.exp2((rows-peak)/temp)))
    return value-2100/8192*v[6]-79174/(1<<21)-255*ROW_COST/(1<<21)+max(0.,-min(v[3:6]))/4096


def gradient(v,gamma,temp):
    g=np.zeros(7);eps=2e-5
    for i in range(7):
        a=v.copy();b=v.copy();a[i]+=eps;b[i]-=eps
        g[i]=(objective(a,gamma,temp)-objective(b,gamma,temp))/(2*eps)
    return objective(v,gamma,temp),g


def optimize(v,gamma):
    for temp in (.005,.0002):
        inverse=np.eye(7);value,grad=gradient(v,gamma,temp)
        for _ in range(100):
            direction=-inverse@grad
            if grad@direction>=0:direction=-grad;inverse=np.eye(7)
            direction/=max(1.,max(abs(direction))/2);rate=1.
            for _ in range(25):
                nxt=np.clip(v+rate*direction,-16,16);step=nxt-v
                if objective(nxt,gamma,temp)<=value+1e-4*grad@step:break
                rate*=.5
            if max(abs(step))<1e-9:break
            nv,ng=gradient(nxt,gamma,temp);delta=ng-grad;curve=delta@step
            if curve>1e-12:
                correction=np.eye(7)-np.outer(step,delta)/curve
                inverse=correction@inverse@correction.T+np.outer(step,step)/curve
            else:inverse=np.eye(7)
            v,value,grad=nxt,nv,ng
    return v


def mixed_row(pop,v,gamma,width=512):
    h=np.exp2(np.r_[0.,v[:6]]);x=2**v[6];a=2**(707/1024);d=2**(-1905/1024);t=1+x/8
    z,logq,cauchy,residual=source_fit(pop,h,a,d,x,1/8)
    poly=defaultdict(float)
    for m in compositions(4,7):
        if any(c and not z[i] for i,c in enumerate(m)):continue
        factor=multinomial(m)*math.prod((z/h)**m)
        hist=tuple(i for i,count in enumerate(m) for _ in range(count))
        for (dest,u,covered,activity,deaths),ways in kernel(hist):
            b=dest.count(3);s=sum(i>=4 for i in dest)
            poly[b,s]+=factor*ways*x**u*t**covered*a**activity*d**deaths*math.prod(h[i] for i in dest)
    total=sum(poly.values());assert abs(math.log2(total)-logq)<1e-7
    poly={k:v/total for k,v in poly.items() if v}
    mass=interaction_mass(poly,width//4,gamma,width)
    return cauchy+math.log2(mass)+gamma*pop[3]*sum(pop[4:])/width,residual


def interaction_mass(poly,groups,gamma,width):
    current=np.ones((1,1))
    for _ in range(groups):
        old=current.shape[0];after=np.zeros((old+4,old+4))
        for (b,s),value in poly.items():after[b:b+old,s:s+old]+=value*current
        current=after
    assert abs(current.sum()-1)<1e-9
    counts=np.arange(width+1)
    penalty=np.exp2(-gamma*counts[:,None]*counts[None,:]/width)
    mass=float(np.sum(current*penalty));assert mass>0
    return mass


def tests():
    from fractions import Fraction as F
    checked=0
    for poly in ({(0,0):F(1,3),(1,1):F(2,3)},
                 {(4,0):F(1,4),(0,4):F(1,4),(1,2):F(1,2)}):
        for groups in range(1,5):
            current={(0,0):F(1)};width=4*groups
            for _ in range(groups):
                after=defaultdict(F)
                for (b,s),p in current.items():
                    for (db,ds),q in poly.items():after[b+db,s+ds]+=p*q
                current=after
            exact=sum(p*F(2)**(-b*s) for (b,s),p in current.items())
            got=interaction_mass({k:float(p) for k,p in poly.items()},groups,width,width)
            assert abs(got-float(exact))<1e-12;checked+=1
    print(dict(exact_interaction_polynomial_checks=checked),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests()
        v=np.array([-.044758,-.042774,-.055881,-1.832615,-1.421249,-1.009882,-.120615])
        proposals=[]
        for gamma in (1.,2.,3.,4.,6.):
            v=optimize(v,gamma);proposals.append((objective(v,gamma),gamma,v.copy()))
            print(dict(gamma=gamma,parameters=v.tolist(),homogeneous_only_log2_envelope=objective(v,gamma)*(1<<21),
                       seconds=time.monotonic()-start),flush=True)
        _,gamma,v=min(proposals)
        profiles=((104,176,72,32,16,16,96),(104,120,64,32,16,16,160),
                  (104,56,32,64,32,32,192),(0,0,0,64,0,0,448),(128,0,0,0,0,0,384))
        for pop in profiles:
            row,residual=mixed_row(pop,v,gamma)
            result=4096*row-537600*v[6]-79174-255*ROW_COST+512*max(0.,-min(v[3:6]))
            print(dict(population=pop,gamma=gamma,all_outgoing=True,float_log2_row_upper=row,
                       repeated_row_diagnostic=result,source_residual=residual,seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

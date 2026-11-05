"""Exact two-update interface for homogeneous seven-type source populations.

The intermediate blocked and selected coordinates share their actual second
partition. Full-size evaluations are FLOAT, not directed certificates. Mixed
B/L2 sources have a separate positive upper relaxation, explicitly labeled.
Serial Peach180s, arrays <64MiB, no files.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import math
import os
from pathlib import Path
import shutil
import sys
import threading
import time
import numpy as np
from quartet_closed_quotient import kernel, NAMES
from quartet_closed_factored import quartet_value
from quartet_closed_population import multinomial
from quartet_two_cohort import compositions


def first_choices(source,a,d,x,y):
    t=1+x*y
    if source in (0,1,2,3):
        result=[(0,x),(1,1),(4,a*t)]
        if source in (1,2):result.append((source+1,1))
        return result
    if source in (4,5):return [(source+1,a*t)]
    return [(3,d*t),(6,a*t)]


def polynomial(source,h,a,d,x,y):
    """Second quartet sum; markers count first-step full targets and Zs."""
    choices=first_choices(source,a,d,x,y);result=defaultdict(lambda:0)
    for counts in compositions(4,len(choices)):
        mid=[0]*7;weight=multinomial(counts)
        for n,(kind,value) in zip(counts,choices):mid[kind]+=n;weight*=value**n
        target=mid[source+1] if source in (1,2) else 0
        result[target,mid[0]]+=weight*quartet_value(mid,h,a,d,x,y)
    return dict(result)


@lru_cache(maxsize=4)
def occupancy(width):
    """Positive coefficients of ((1+z)^4-1)^g, for every needed g."""
    groups=width//4;table=np.zeros((groups+1,width+1));table[0,0]=1.
    for g in range(1,groups+1):
        for j in range(1,5):table[g,j:]+=math.comb(4,j)*table[g-1,:width+1-j]
    return table


def first_partition_factor(source,width):
    groups=width//4;occ=occupancy(width);answer=np.zeros((width+1,width+1))
    targets=range(0,width+1,4) if source in (1,2) else (0,)
    for target in targets:
        f=target//4;left=groups-f
        cluster=math.comb(groups,f)/math.comb(width,target)
        for zero in range(left,4*left+1):
            answer[target,zero]=cluster*occ[left,zero]/math.comb(width-target,zero)
    assert answer.max()<=1+1e-12
    return answer


def positive_power(poly,groups):
    total=sum(poly.values());assert total>0
    has_target=any(k[0] and value for k,value in poly.items())
    current=np.ones((1,1))
    for _ in range(groups):
        rows,cols=current.shape
        after=np.zeros((rows+4 if has_target else 1,cols+4))
        for (t,z),value in poly.items():
            if value:after[t:t+rows,z:z+cols]+=(float(value)/total)*current
        current=after
    assert abs(current.sum()-1)<1e-9
    return current,groups*math.log2(total)


def row(source,h,a,d,x,y,width=512):
    poly=polynomial(source,h,a,d,x,y)
    if source in (0,4,5,6):return width/4*math.log2(sum(poly.values()))-width*math.log2(h[source])
    current,scale=positive_power(poly,width//4)
    factor=first_partition_factor(source,width)[:current.shape[0],:]
    mass=float(np.sum(current*factor));assert mass>0
    return scale+math.log2(mass)-width*math.log2(h[source])


def row_gradient(source,v,width=512,need_gradient=True):
    """One positive convolution, then a coefficient-adjoint derivative.

    Local derivatives propose FLOAT potentials only. They never determine
    whether a final directed certificate accepts.
    """
    a=2**(707/1024);d=2**(-1905/1024);y=1/8;groups=width//4
    def local(v):return polynomial(source,np.exp2(np.r_[0.,v[:6]]),a,d,2**v[6],y)
    poly=local(v);total=sum(poly.values());keys=list(poly)
    if source in (0,4,5,6):adjoint={key:1. for key in keys}
    else:
        current,_=positive_power(poly,groups-1)
        factor=first_partition_factor(source,width);rows,cols=current.shape
        adjoint={key:float(np.sum(current*factor[key[0]:key[0]+rows,key[1]:key[1]+cols])) for key in keys}
    mass=sum(adjoint[key]*poly[key] for key in keys)
    value=groups*math.log2(total)+math.log2(mass/total)-(width*v[source-1] if source else 0.)
    if not need_gradient:return value,None
    grad=np.zeros(7);eps=2e-5
    for i in range(7):
        plus=v.copy();minus=v.copy();plus[i]+=eps;minus[i]-=eps
        lp=local(plus);lm=local(minus)
        deriv=sum(adjoint[key]*(lp[key]-lm[key])/(2*eps) for key in keys)
        grad[i]=groups*deriv/(math.log(2)*mass)
    if source:grad[source-1]-=width
    return value,grad


def optimize_homogeneous(v):
    """Restricted homogeneous sources, NOT an all-population certificate."""
    start=time.monotonic()
    constant=-79174/(1<<21)-255*math.log2(math.comb(8192,2100))/(1<<21)
    def evaluate(v,temp,need_gradient=True):
        results=[row_gradient(i,v,need_gradient=need_gradient) for i in range(7)]
        rows=np.array([r[0]/1024 for r in results]);peak=max(rows)
        weights=np.exp2((rows-peak)/temp);weights/=sum(weights)
        value=peak+temp*math.log2(sum(np.exp2((rows-peak)/temp)))-2100/8192*v[6]+constant
        tail=int(np.argmin(v[3:6]))+3
        value+=max(0.,-v[tail])/4096
        if not need_gradient:return value,None,rows
        grad=weights@np.array([r[1]/1024 for r in results]);grad[6]-=2100/8192
        if v[tail]<0:grad[tail]-=1/4096
        return value,grad,rows
    for temp in (.002,.0002):
        inverse=np.eye(7);value,grad,rows=evaluate(v,temp)
        for iteration in range(40):
            direction=-inverse@grad
            if grad@direction>=0:direction=-grad;inverse=np.eye(7)
            direction/=max(1.,max(abs(direction))/2);rate=1.
            for _ in range(15):
                nxt=np.clip(v+rate*direction,-16,16);step=nxt-v
                nv,_,_=evaluate(nxt,temp,False)
                if nv<=value+1e-4*grad@step:break
                rate*=.5
            if max(abs(step))<1e-8:break
            nv,ng,rows=evaluate(nxt,temp);delta=ng-grad;curve=delta@step
            if curve>1e-12:
                correction=np.eye(7)-np.outer(step,delta)/curve
                inverse=correction@inverse@correction.T+np.outer(step,step)/curve
            else:inverse=np.eye(7)
            v,value,grad=nxt,nv,ng
            if iteration%10==9:print(dict(iteration=iteration+1,temp=temp,smooth_expression=value*(1<<21),seconds=time.monotonic()-start),flush=True)
        _,_,rows=evaluate(v,temp,False)
        sharp=max(rows)-2100/8192*v[6]+constant+max(0.,-min(v[3:6]))/4096
        print(dict(homogeneous_only=True,temp=temp,parameters=v.tolist(),float_log2_expression=sharp*(1<<21),
                   row_bits=(rows*1024).tolist(),seconds=time.monotonic()-start),flush=True)
    return v


def mixed_blocked_selected(blocked,h,a,d,x,y,width=512):
    """Retain the first B/L2 partition; relax first non-full B quartets."""
    if blocked%4:return -math.inf
    t=1+x*y;ordinary=((0,x),(1,1),(4,a*t));live=((3,d*t),(6,a*t))
    coeff=np.zeros(5)
    for b in range(5):
        for left in compositions(b,3):
            for right in compositions(4-b,2):
                mid=[0]*7;weight=math.comb(4,b)*multinomial(left)*multinomial(right)
                for n,(kind,value) in zip(left+right,ordinary+live):mid[kind]+=n;weight*=value**n
                coeff[b]+=weight*quartet_value(mid,h,a,d,x,y)
    total=sum(coeff);p=coeff/total;current=np.array([1.])
    for _ in range(width//4):current=np.convolve(current,p)
    assert abs(sum(current)-1)<1e-9
    if not current[blocked]:return -math.inf
    partition=math.comb(width//4,blocked//4)/math.comb(width,blocked)
    return (math.log2(partition)+width/4*math.log2(total)+math.log2(current[blocked])
            -math.log2(math.comb(width,blocked))-blocked*math.log2(h[3])-(width-blocked)*math.log2(h[6]))


def exact_local(hist,a,d,x,y):
    out=defaultdict(F)
    for (dest,u,v,activity,deaths),ways in kernel(hist):
        pop=tuple(dest.count(i) for i in range(7))
        out[pop]+=ways*x**u*(1+x*y)**v*a**activity*d**deaths
    return dict(out)


def exact_transition(pop,a,d,x,y):
    width=sum(pop);assert width in (4,8)
    if width==4:return exact_local(tuple(i for i,n in enumerate(pop) for _ in range(n)),a,d,x,y)
    out=defaultdict(F)
    active=[i for i,n in enumerate(pop) if n]
    for counts in compositions(4,len(active)):
        left=[0]*7
        for i,n in zip(active,counts):left[i]=n
        if any(n>p for n,p in zip(left,pop)):continue
        right=tuple(p-n for n,p in zip(left,pop))
        prob=F(math.prod(math.comb(p,n) for n,p in zip(left,pop)),math.comb(width,4))
        lh=tuple(i for i,n in enumerate(left) for _ in range(n));rh=tuple(i for i,n in enumerate(right) for _ in range(n))
        for l,lv in exact_local(lh,a,d,x,y).items():
            for r,rv in exact_local(rh,a,d,x,y).items():out[tuple(u+v for u,v in zip(l,r))]+=prob*lv*rv
    return dict(out)


def exact_one(pop,h,a,d,x,y):
    width=sum(pop)
    if width==4:return quartet_value(pop,h,a,d,x,y)
    active=[i for i,n in enumerate(pop) if n];value=F(0)
    for counts in compositions(4,len(active)):
        left=[0]*7
        for i,n in zip(active,counts):left[i]=n
        if any(n>p for n,p in zip(left,pop)):continue
        right=tuple(p-n for n,p in zip(left,pop))
        prob=F(math.prod(math.comb(p,n) for n,p in zip(left,pop)),math.comb(width,4))
        value+=prob*quartet_value(left,h,a,d,x,y)*quartet_value(right,h,a,d,x,y)
    return value


def exact_formula(source,h,a,d,x,y,width):
    poly=polynomial(source,h,a,d,x,y);current={(0,0):F(1)};groups=width//4
    for _ in range(groups):
        after=defaultdict(F)
        for (t,z),value in current.items():
            for (dt,dz),weight in poly.items():after[t+dt,z+dz]+=value*weight
        current=dict(after)
    occ=[[1]]
    for g in range(1,groups+1):
        after=[0]*(4*g+1)
        for i,value in enumerate(occ[-1]):
            for j in range(1,5):after[i+j]+=value*math.comb(4,j)
        occ.append(after)
    answer=F(0)
    for (target,z),value in current.items():
        if source in (0,4,5,6):factor=1
        elif target%4:continue
        else:
            f=target//4;left=groups-f
            if z>4*left:continue
            factor=F(math.comb(groups,f)*occ[left][z],math.comb(width,target)*math.comb(width-target,z))
        answer+=value*factor
    return answer/h[source]**width


def tests():
    h=tuple(map(F,(1,2,3,4,5,6,7)));a=F(3,2);d=F(1,4);x=F(2,3);y=F(1,8)
    checks=0
    for width in (4,8):
        for source in range(7):
            pop=tuple(width*int(i==source) for i in range(7))
            exact=sum(value*exact_one(mid,h,a,d,x,y) for mid,value in exact_transition(pop,a,d,x,y).items())/h[source]**width
            assert exact==exact_formula(source,h,a,d,x,y,width)
            got=row(source,tuple(map(float,h)),float(a),float(d),float(x),float(y),width)
            assert abs(math.log2(exact)-got)<1e-10;checks+=1
        for blocked in range(0,width+1,2):
            pop=(0,0,0,blocked,0,0,width-blocked)
            exact=sum(value*exact_one(mid,h,a,d,x,y) for mid,value in exact_transition(pop,a,d,x,y).items())/h[3]**blocked/h[6]**(width-blocked)
            got=mixed_blocked_selected(blocked,tuple(map(float,h)),float(a),float(d),float(x),float(y),width)
            assert (exact==0 and got==-math.inf) or (exact>0 and math.log2(exact)<=got+1e-10);checks+=1
    print(dict(exact_two_update_checks=checks),flush=True)
    v=np.array([.1,-.2,.3,-.4,.5,-.6,.7]);derivatives=0
    for source in range(7):
        value,grad=row_gradient(source,v,width=8)
        for i in range(7):
            plus=v.copy();minus=v.copy();plus[i]+=1e-5;minus[i]-=1e-5
            numeric=(row_gradient(source,plus,width=8,need_gradient=False)[0]-row_gradient(source,minus,width=8,need_gradient=False)[0])/2e-5
            assert abs(numeric-grad[i])<1e-7;derivatives+=1
    print(dict(coefficient_adjoint_checks=derivatives),flush=True)


def quota_tests():
    """Degree+1 rational evaluations prove the small quota-polynomial identities."""
    h=tuple(map(F,(1,2,3,4,5,6,7)));a=F(3,2);d=F(1,4);y=F(1,8);checks=0
    for width in (4,8):
        for source in range(7):
            pop=tuple(width*int(i==source) for i in range(7))
            for node in range(2*width+1):
                x=F(node,4)
                direct=sum(value*exact_one(mid,h,a,d,x,y) for mid,value in exact_transition(pop,a,d,x,y).items())/h[source]**width
                assert direct==exact_formula(source,h,a,d,x,y,width);checks+=1
    print(dict(exact_quota_polynomial_evaluations=checks,homogeneous_profiles=14,
               degree_bounds=[8,16],all_small_original_count_coefficients=True),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests()
        if '--quota-tests' in sys.argv:quota_tests();return
        logs=np.array([0.,-.04475834141004524,-.042774336422858605,-.05588058545549474,
                       -1.8326154305711793,-1.4212486763304475,-1.0098819495674625])
        logx=-.12061508248717494;h=np.exp2(logs);a=2**(707/1024);d=2**(-1905/1024);x=2**logx;y=1/8
        if '--optimize' in sys.argv:
            v=optimize_homogeneous(np.r_[logs[1:],logx]);logs=np.r_[0.,v[:6]];logx=v[6];h=np.exp2(logs);x=2**logx
        constant=-537600*logx-79174-255*math.log2(math.comb(8192,2100))+512*max(0.,-min(logs[4:]))
        for source in range(7):
            value=row(source,h,a,d,x,y)
            old=128*math.log2(quartet_value(tuple(4*int(i==source) for i in range(7)),h,a,d,x,y))-512*logs[source]
            print(dict(source=NAMES[source],float_exact_pair_log2=value,twice_one_step=2*old,
                       repeated_pair_diagnostic=2048*value+constant,seconds=time.monotonic()-start),flush=True)
        for blocked in (16,32,64,128,256,384):
            value=mixed_blocked_selected(blocked,h,a,d,x,y)
            print(dict(blocked=blocked,selected=512-blocked,float_pair_log2_upper=value,
                       repeated_pair_diagnostic=2048*value+constant,seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

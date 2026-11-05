"""Bounded warm-start proof-weight proposals; exact replay decides acceptance.

No construction changes or output files. Each proposal has a wall-time cap;
the local operators have at most 4*81 entries (not a history search tree).
"""
import argparse
import json
import math
import time
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
from layered_outer_locations import hardest,rows,start,objective,evaluate,M,b
from layered_general_components import mixed_dimension,envelopes
import cover_layered_location_regions as cover
from layered_interpolation import LogWeight,exact_tau


def profile_at(R,offsets,refined=False):
    counts=tuple(v+d for v,d in zip(hardest()[1:6],offsets))
    dimension=mixed_dimension(R,(R*M-b)//2,*counts[:4],cover.region_env(R,refined))
    return (R,)+counts+(dimension,0)


def refine(profile,s,m,seconds=25,seed=None,free_m=False,refined=False):
    assert 129<=profile[0]<=256 and 1<=s<=M and 0<=float(m)<=32
    assert 0<seconds<=60
    began=time.monotonic();candidates=[]
    for i,(name,_,two,p,markers) in enumerate(cover.BANK):
        if not two or (seed is not None and name!=seed):continue
        x=np.append(start(markers),math.log2(p/(1-p)))
        candidates.append((evaluate(profile,s,m,x,True),i,x))
    assert candidates,'no matching two-update seed'
    initial,seed_index,latest=min(candidates,key=lambda v:v[0]);latest=latest.copy()
    if free_m:latest=np.append(latest,float(m))
    def weight(x):return x[21] if free_m else m
    rr,ee=rows(latest,s,m,True);latest[14:18]=rr.max(axis=1);latest[18:20]=ee.max(axis=1)
    original=latest.copy();calls=0
    def check():
        if time.monotonic()-began>seconds:raise TimeoutError('bounded warm-start proposal')
    def fn(x):
        nonlocal latest,calls
        check();calls+=1;latest=x.copy()
        return objective(x,profile,s,weight(x),True)
    def cons(x):
        check();rr,ee=rows(x,s,weight(x),True)
        return np.concatenate(((x[14:18,None]-rr).ravel(),(x[18:20,None]-ee).ravel()))
    try:
        result=minimize(fn,latest,method='SLSQP',constraints=[dict(type='ineq',fun=cons)],
            bounds=[(-256,256)]*5+[(0,128)]+[(-128,128)]*8+[(-2048,2048)]*6+[(-16,24)]+([(0,32)] if free_m else []),
            options=dict(maxiter=180,ftol=1e-12))
        latest=result.x;stop=str(result.message)
    except TimeoutError as error:stop=str(error)
    # An infeasible or timed-out numerical iterate is only a proposal.
    # Recompute every row maximum before selecting it; exact replay is below.
    x=min((original,latest),key=lambda v:evaluate(profile,s,weight(v),v,True))
    if free_m:m=LogWeight(F(2.**x[21]))
    markers=[F(2.**(-v if i<6 else v)) for i,v in enumerate(x[:14])]
    colour=F(1/(1+2.**(-x[20])))
    name=f'R{profile[0]}s{s}m{float(m):.6g}'
    index=len(cover.BANK);cover.BANK.append((name,m,True,colour,markers))
    point=profile[1:6]
    exact=cover.location_box(profile[0],point,point,index,s,s,refined)
    tau=exact_tau(m)
    return dict(name=name,profile=profile,s=s,m=float(m),tau=[str(tau.numerator),str(tau.denominator)],two=True,refined_route=refined,seed=cover.BANK[seed_index][0],
        initial_float=initial,float_upper=evaluate(profile,s,m,x,True),exact_upper=exact,
        calls=calls,seconds=time.monotonic()-began,stop=stop,
        colour=[str(colour.numerator),str(colour.denominator)],
        markers=[[str(v.numerator),str(v.denominator)] for v in markers])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--R',type=int,default=192)
    p.add_argument('--s',type=int,required=True);p.add_argument('--m',type=float,default=25)
    p.add_argument('--seconds',type=float,default=25);p.add_argument('--seed');p.add_argument('--free-m',action='store_true');p.add_argument('--refined-route',action='store_true')
    p.add_argument('--matching-moments',action='store_true')
    p.add_argument('--offsets',nargs=5,type=int,default=(128,-128,-128,128,-128))
    a=p.parse_args()
    weight=int(a.m) if float(a.m).is_integer() else LogWeight(F(2.**a.m))
    route='factorial' if a.matching_moments else a.refined_route
    print(json.dumps(refine(profile_at(a.R,a.offsets,route),a.s,weight,a.seconds,a.seed,a.free_m,route)))

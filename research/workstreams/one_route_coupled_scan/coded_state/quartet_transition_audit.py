"""Directed isolated-cell audit of the averaged-burst reference kernel.

No full-distance claim. Serial,90s watchdog,32MiB numeric-array budget,
20GiB free, no files. FLOAT saddlepoints propose rational coefficient tilts;
positive integer lower/upper DPs and directed logs verify each selected cell.
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

import quartet_cut_moment as cut
import quartet_full_weight_kernel as full
from quartet_two_update import log_upper


def root_interval(numerator,denominator,bits=80):
    target=1<<(bits*denominator-numerator);lo=0;hi=1<<bits
    while hi-lo>1:
        middle=(lo+hi)//2
        if middle**denominator>=target:hi=middle
        else:lo=middle
    assert lo**denominator<=target<=hi**denominator
    return F(lo,1<<bits),F(hi,1<<bits)


def tilt(terms,groups,u,v):
    points=np.array(list(terms),dtype=float)
    logs=np.array([math.log(float(x)) for x in terms.values()])
    target=np.array([u/groups,v/groups]);xy=np.zeros(2)
    for _ in range(60):
        weights=logs+points@xy;weights=np.exp(weights-max(weights));weights/=sum(weights)
        mean=weights@points;centered=points-mean
        covariance=(centered.T*weights)@centered
        step=np.linalg.solve(covariance+np.eye(2)*1e-12,mean-target)
        if np.max(abs(mean-target))<1e-10:break
        xy-=step/max(1.,np.max(abs(step))/2.)
    assert np.max(abs(mean-target))<1e-7
    return tuple(F(float(math.exp(z))).limit_denominator(1<<20) for z in xy)


def quota_tilt(weight,rows=2100,logb=8):
    """FLOAT proposal satisfying incoming, outgoing, and label means jointly."""
    gamma=2**(weight*logb/512);phi=2**(-37*logb/50);base=2.**logb;records={}
    for a in range(5):
        before=(1<<a)-1
        for after in range(16):
            v=after.bit_count();h=(before&after).bit_count();d=a-h
            for target in range(16):
                if a==0 and target&after:continue
                if a and 100*(target&~after).bit_count()>=25*a+74*h:continue
                t=target.bit_count();c=t if a else (target&after).bit_count()
                value=math.comb(4,a)*gamma**t*base**(-c)*phi**d
                key=(a,v,t);records[key]=records.get(key,0.)+value
    points=np.array(list(records),dtype=float);logs=np.log(np.array(list(records.values())))
    target=np.array([weight/128,weight/128,rows/2048]);xyz=np.zeros(3)
    for _ in range(80):
        probabilities=logs+points@xyz;probabilities=np.exp(probabilities-max(probabilities));probabilities/=sum(probabilities)
        mean=probabilities@points;centered=points-mean
        covariance=(centered.T*probabilities)@centered
        step=np.linalg.solve(covariance+np.eye(3)*1e-12,mean-target)
        if np.max(abs(mean-target))<1e-10:break
        xyz-=step/max(1.,np.max(abs(step))/2.)
    assert np.max(abs(mean-target))<1e-7
    odds=math.exp(xyz[2]);p=F(float(odds/(1+odds))).limit_denominator(1<<20)
    return p,mean


def coefficient_interval(lower,upper,groups,u,v,x,y):
    scale=1<<31
    tilted={key:value*x**key[0]*y**key[1] for key,value in upper.items()}
    normalization=sum(tilted.values())
    probabilities=[]
    for key,hi in tilted.items():
        a,w=key
        lo=lower[key]*x**a*y**w/normalization*scale
        hi=hi/normalization*scale
        probabilities.append((a,w,lo.numerator//lo.denominator,
                              (hi.numerator+hi.denominator-1)//hi.denominator))
    assert all(0<=lo<=hi<=scale for _,_,lo,hi in probabilities)
    assert (u+1)*(v+1)*8*12<32*(1<<20)
    low=np.zeros((u+1,v+1),dtype=np.uint64);high=np.zeros_like(low)
    low[0,0]=high[0,0]=scale
    for step in range(groups):
        low_next=np.zeros_like(low);high_next=np.zeros_like(high)
        for a,w,lo,hi in probabilities:
            rows=min(4*step+1,u+1-a);cols=min(4*step+1,v+1-w)
            if rows<=0 or cols<=0:continue
            low_next[a:a+rows,w:w+cols]+=(low[:rows,:cols]*np.uint64(lo))>>31
            high_next[a:a+rows,w:w+cols]+=(high[:rows,:cols]*np.uint64(hi)+np.uint64(scale-1))>>31
        low,high=low_next,high_next
        # Bounds ensure every product stays below2^63, and sums below2^64.
        assert int(high.max())<1<<32
    return F(int(low[u,v]),scale),F(int(high[u,v]),scale),normalization


def tests():
    p=F(2,5);b=F(16);canonical=(F(1,4),F(37,50))
    terms=cut.local(p,b,F(1),F(1),F(1,8),canonical,label_marker=F(4))
    exact={}
    for (a,v),first in terms.items():
        for (c,w),second in terms.items():
            exact[a+c,v+w]=exact.get((a+c,v+w),F(0))+first*second
    checks=0;x=F(1,2);y=F(2)
    for u in range(1,8):
        for v in range(1,8):
            lo,hi,z=coefficient_interval(terms,terms,2,u,v,x,y)
            normal=x**u*y**v/z**2
            truth=exact[u,v]*normal
            assert lo<=truth<=hi;checks+=1
    low,high=root_interval(148,25)
    assert low**25<=F(1,1<<148)<=high**25
    print(dict(exact_directed_coefficient_cells=checks,integer_root_checks=2,
               arithmetic='positive integer lower/upper DP',files_written=0),flush=True)


def audit(p=F(525,2048),quota=False,logb=8):
    cells=((128,128),(256,256),(384,384)) if quota else ((128,128),(256,256),(384,384),(128,384),(384,128))
    assert logb in (4,8)
    b=F(2)**logb;tau=F(37,50);canonical=(F(1,4),tau)
    phi_lo,phi_hi=root_interval(37*logb//2,25)
    anchors=[]
    if not quota:
        for gamma in (1,2**(logb//2),2**logb):
            terms=cut.local(float(p),float(b),1.,1.,float(phi_hi),canonical,label_marker=gamma)
            anchors.append(full.power(terms,128))
    for u,v in cells:
        if quota:p,mean=quota_tilt(u,logb=logb)
        k=F(max(u,v),512);assert (max(u,v)*logb)%512==0
        gamma=F(2)**(max(u,v)*logb//512)
        lower=cut.local(p,b,F(1),F(1),phi_lo,canonical,label_marker=gamma)
        upper=cut.local(p,b,F(1),F(1),phi_hi,canonical,label_marker=gamma)
        x,y=tilt(upper,128,u,v)
        lo,hi,z=coefficient_interval(lower,upper,128,u,v,x,y)
        assert lo>0
        ingredients=((128,z),(-u,x),(-v,y),(-1,math.comb(512,u)))
        lower_log=-log_upper(-1,lo)-sum(log_upper(-coefficient,value) for coefficient,value in ingredients)
        upper_log=log_upper(1,hi)+sum(log_upper(coefficient,value) for coefficient,value in ingredients)
        correction=logb*tau*max(u-v,0)
        lower_log+=correction;upper_log+=correction
        chord=None
        if not quota:
            segment=0 if k<=F(1,2) else 1;fraction=float(2*k-segment)
            chord=(1-fraction)*anchors[segment][u,v]+fraction*anchors[segment+1][u,v]+float(correction)
        regional=None
        if quota:
            regional=16*upper_log+log_upper(-1,math.comb(8192,2100))+log_upper(-2100,p)+log_upper(-6092,1-p)
        print(dict(cell=[u,v],proof_p=str(p),reference_fraction=str(k),log2_base=logb,
                   directed_log2_lower=str(lower_log),directed_log2_upper=str(upper_log),
                   lower_float=float(lower_log),upper_float=float(upper_log),
                   directed_width_bits=float(upper_log-lower_log),
                   chord_float=chord,chord_excess_over_verified_upper=None if chord is None else float(chord-upper_log),
                   exact_quota_constant_weight_region_log2_upper=None if regional is None else str(regional),
                   exact_quota_region_upper_float=None if regional is None else float(regional),
                   rational_tilt=[str(x),str(y)],
                   includes_excess_death_correction=True,
                   regional_quota_not_yet_imposed=not quota,
                   fixed_weight_at_every_regional_boundary=quota,
                   whole_code_distance=False,files_written=0),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser();parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--quota',action='store_true')
        parser.add_argument('--logb',type=int,choices=(4,8),default=8)
        args=parser.parse_args();tests()
        if not args.tests_only:audit(quota=args.quota,logb=args.logb)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

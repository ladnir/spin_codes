"""Bounded FLOAT screen of a quota-aware regional read-four envelope.

These numerical values are not directed certificates. Construction fixed.
Serial,90s watchdog,O(r) stored numbers,no files; run on Peach only.
"""
import math
import os
from pathlib import Path
import shutil
import threading
import time


def logchoose(n,k):
    if not 0<=k<=n:return -math.inf
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1))/math.log(2)


def column_moment(p,logb,slope,tau,steps):
    assert 2*slope<=tau<=3*slope
    wa=2**(logb*(3*slope-tau));wl=2**(logb*slope)
    p4=p**4
    state=(1-p,p,0.0,0.0,0.0)
    for _ in range(steps-1):
        total=sum(state)
        state=((1-p)*total,(p-p4)*total+p4*state[0],
               p4*state[1],p4*state[2],p4*(wa*state[3]+wl*state[4]))
    return sum(v*2**(logb*reward) for v,reward in zip(state,(0,slope,2*slope,3*slope,slope+tau)))


def golden_minimum(function,lo,hi,iterations=20):
    ratio=(math.sqrt(5)-1)/2
    x=hi-ratio*(hi-lo);y=lo+ratio*(hi-lo);fx=function(x);fy=function(y)
    for _ in range(iterations):
        if fx<fy:
            hi=y;y=x;fy=fx;x=hi-ratio*(hi-lo);fx=function(x)
        else:
            lo=x;x=y;fx=fy;y=lo+ratio*(hi-lo);fy=function(y)
    return min(fx,fy)


def regional_envelope(width,steps,r,logb,slope,tau,read=4):
    m=width*steps;result=[]
    for c in range(r+1):
        active=m-r+c;unpaid=r-c
        if unpaid==0:
            result.append(logb*slope*m);continue
        def value(p):
            moment=column_moment(p,read*logb,slope,tau,steps)
            # Uniform fixed-unpaid moment <= independent moment / point mass.
            return width/read*math.log2(moment)-active*math.log2(p)-unpaid*math.log2(1-p)
        # Any evaluated p gives an upper bound. Search failure affects
        # sharpness, not its mathematical direction (before rounding).
        best=golden_minimum(value,1e-6,1-1e-6)
        result.append(logchoose(r,c)-logchoose(m,unpaid)+best)
    return result


def coefficient_upper(logcoeff,regions,paid):
    target=paid/regions
    def evaluate(marker):
        values=[value+c*marker for c,value in enumerate(logcoeff)]
        peak=max(values);weights=[2**(value-peak) for value in values]
        total=sum(weights)
        mean=sum(c*weight for c,weight in enumerate(weights))/total
        return regions*(peak+math.log2(total))-paid*marker,mean
    lo=-64.;hi=64.
    for _ in range(64):
        mid=(lo+hi)/2;value,mean=evaluate(mid)
        if mean<target:lo=mid
        else:hi=mid
    upper,mean=evaluate((lo+hi)/2)
    return upper,(lo+hi)/2,mean


def tests():
    # No moments: the conditional bound at p=active/m is deliberately
    # loose, but never below the exact unit moment.
    for p in (.2,.5,.8):
        assert abs(column_moment(p,0,1,2.5,16)-1)<1e-12
    assert column_moment(1,2,1,2.5,16)==4**16
    from quartet_regional_transfer import regional_coefficients
    exact=regional_coefficients(8,2,1,4)
    envelope=regional_envelope(8,2,4,2,1,2.5)
    assert all(upper+1e-9>=math.log2(float(value)) for value,upper in zip(exact,envelope))
    print(dict(float_envelope_small_exact_checks=5,directed_certificate=False),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests()
        slope=707/2560;tau=381/512;B=165556;c=111572
        row_union=logchoose(8192,2100)
        for logb in (1,2,3,4,6,8):
            upper=regional_envelope(512,16,2100,logb,slope,tau)
            coefficient,marker,mean=coefficient_upper(upper,256,c)
            print(dict(log2_base=logb,log2_paid_marker=marker,
                       float_log2_tail_expression=coefficient-logb*B+row_union,
                       paid_mean_per_region=mean,read_degree=4,
                       original_quota_retained_in_each_coefficient=True,
                       free_boundary_per_region=True,directed_certificate=False,
                       whole_code_distance=False),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

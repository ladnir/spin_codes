"""FLOAT exact-coefficient diagnostic for one fixed state-weight family.

Retain incoming weight, outgoing weight, and inactive count. Use complement
weights near full support, so a 23-hole profile needs only 24x24x513 cells.
Exact combinatorial coefficients are evaluated in FLOAT log arithmetic.
Serial, 120s watchdog, 32 MiB numeric-array budget, >=20GiB free, no files.
Not a directed certificate or whole-code coverage.
"""
import argparse
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np

import quartet_full_weight_kernel as full
import quartet_quiet_regions as quiet
from quartet_joint_quota import log_convolve


def terms(p,theta,ell,complement):
    result={}
    for a in range(5):
        before=(1<<a)-1
        for after in range(16):
            v=after.bit_count();d=(before&~after).bit_count()
            for target in range(16):
                coefficients=full.q.choices(4,before,after,target,True)
                if not any(coefficients):continue
                t=target.bit_count();c=t if a else (after&target).bit_count()
                poly=sum(n*ell**j for j,n in enumerate(coefficients))
                value=math.comb(4,a)*p**t*(1-p)**(4-t)*(poly*2**(-32*c))**theta
                value*=ell**(-theta*a)*((1<<32)*ell)**(-theta*d)
                key=(4-a,4-v,t) if complement else (a,v,t)
                result[key]=result.get(key,0.)+value
    return result


def profile(p,theta,ell,weight,groups=128):
    size=groups*4;complement=weight>size//2;cap=min(weight,size-weight)
    assert cap<=32
    local=terms(p,theta,ell,complement)
    assert (cap+1)**2*(size+1)*8*4<32*(1<<20)
    current=np.full((cap+1,cap+1,size+1),-np.inf);current[0,0,0]=0
    for step in range(groups):
        nxt=np.full_like(current,-np.inf);length=4*step+1
        for (a,v,t),value in local.items():
            na=min(cap-a,4*step)+1;nv=min(cap-v,4*step)+1
            if na<=0 or nv<=0:continue
            nxt[a:a+na,v:v+nv,t:t+length]=np.logaddexp2(
                nxt[a:a+na,v:v+nv,t:t+length],current[:na,:nv,:length]+math.log2(value))
        current=nxt
    return current[cap,cap,:]-quiet.log_choose(size,weight)


def tests():
    p=.37;theta=.175;ell=.00377
    kernel=full.power(full.local(p,theta,ell,1,1),2)
    for w in range(1,9):
        coefficients=profile(p,theta,ell,w,2)
        assert abs(np.logaddexp2.reduce(coefficients)-kernel[w,w])<1e-11
        alternate=profile(.6,theta,ell,w,2)
        assert np.array_equal(np.isfinite(coefficients),np.isfinite(alternate))
        valid=np.isfinite(coefficients);t=np.arange(9)[valid]
        raw=coefficients[valid]-t*math.log2(p)-(8-t)*math.log2(1-p)
        other=alternate[valid]-t*math.log2(.6)-(8-t)*math.log2(.4)
        assert np.max(np.abs(raw-other))<1e-11
    print(dict(two_group_diagonal_checks=8,artificial_p_invariance_checks=8),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(120,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser();parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--weight',type=int,default=489)
        parser.add_argument('--p',type=float,default=2100/8192)
        args=parser.parse_args();tests()
        if args.tests_only:return
        w=args.weight;p=args.p;r=2100;theta=7/40;ell=.00377
        assert 1<=w<=512 and min(w,512-w)<=32 and 0<p<1
        coefficients=profile(p,theta,ell,w)
        print(dict(phase='joint profile coefficients',weight=w,seconds=time.monotonic()-start),flush=True)
        birth=np.array([quiet.log_choose(512,t)+quiet.log_choose(512-t,w)
                        +t*math.log2(p)+(512-t)*math.log2(1-p)
                        if t<=512-w else -np.inf for t in range(513)])
        power=np.array([0.]);first=None
        for step in range(16):
            if step==15:first=log_convolve(birth,power,r)
            power=log_convolve(power,coefficients,r)
        pi=quiet.log_choose(8192,r)+r*math.log2(p)+(8192-r)*math.log2(1-p)
        base=quiet.log_choose(8192,r)+theta+32*theta*(943718-1048576+128*r)+theta*943718*math.log2(ell)
        result=base+(first[r]-pi)+255*(power[r]-pi)-w*theta*math.log2(ell)
        raw=coefficients-np.arange(513)*math.log2(p)-(512-np.arange(513))*math.log2(1-p)
        cauchy=[]
        for i in range(9):
            odds=(r/(8192-r))*2**(32*theta*i/8);alt=odds/(1+odds)
            diagonal=np.logaddexp2.reduce(raw+np.arange(513)*math.log2(alt)
                                         +(512-np.arange(513))*math.log2(1-alt))
            altpi=quiet.log_choose(8192,r)+r*math.log2(alt)+(8192-r)*math.log2(1-alt)
            value=base+quiet.log_choose(512,w)+w*math.log2(1-alt)+4095*diagonal-256*altpi-w*theta*math.log2(ell)
            cauchy.append((float(value),alt))
        print(dict(log2_bound=float(result),weight=w,inactive_rows=r,
                   matched_cauchy_grid_best=min(cauchy),
                   first_activation=1,all_turnovers=True,exact_quota_extraction=True,
                   directed=False,whole_code=False,seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

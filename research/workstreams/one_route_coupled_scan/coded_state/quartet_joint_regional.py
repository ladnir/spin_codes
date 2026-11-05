"""Composable two-cohort bucket transfer, FLOAT diagnostic only.

Serial Peach,90s (1200s full band),<8GiB arrays,no files.
Not a distance certificate.
"""
import math
import os
from pathlib import Path
import shutil
import sys
import threading
import time
from functools import lru_cache
import numpy as np
import quartet_fresh_count_screen as fresh
import quartet_fresh_regional as regional
from quartet_history_screen import occupation,logchoose

BIN=32
LOGZ=-2.


@lru_cache(maxsize=64)
def states(a):
    result=[];lookup=np.full((a//BIN+1,a//BIN+1),-1,dtype=int)
    for rlo in range(0,a+1,BIN):
        for slo in range(0,a-rlo+1,BIN):
            lookup[rlo//BIN,slo//BIN]=len(result)
            result.append((rlo,min(a,rlo+BIN-1),slo,min(a-rlo,slo+BIN-1)))
    return result,lookup


def vertices(a,box):
    rl,rh,sl,sh=box
    points={(r,s) for r in (rl,rh) for s in (sl,sh)}
    points.update((r,a-r) for r in (rl,rh));points.update((a-s,s) for s in (sl,sh))
    return [(r,s) for r,s in points if rl<=r<=rh and sl<=s<=sh and r+s<=a]


def old_coefficients(count,second_bins):
    """Uniform over n2 in its bin and all divisions between n3,n4."""
    if count==0:return np.zeros((1,1))
    d=fresh.D;h=fresh.H
    lower=np.array([lo for lo,hi in second_bins])[:,None]
    upper=np.array([hi for lo,hi in second_bins])[:,None]
    degree=np.arange(count+1,dtype=float)[None,:]
    left=np.full((len(second_bins),count+1),-80.);right=np.full_like(left,80.)
    for _ in range(50):
        x=(left+right)/2;z=np.exp2(x)
        l2=np.log2((1+d[2]*z)/h[2]);l3=np.log2((1+d[3]*z)/h[3]);l4=np.log2((1+d[4]*z)/h[4])
        selected_d=np.where(l3>=l4,d[3],d[4]);old=np.maximum(l3,l4)
        s=np.where(l2>=old,upper,lower)
        derivative=s*d[2]*z/(1+d[2]*z)+(count-s)*selected_d*z/(1+selected_d*z)-degree
        left=np.where(derivative<0,x,left);right=np.where(derivative>=0,x,right)
    x=(left+right)/2;z=np.exp2(x)
    l2=np.log2((1+d[2]*z)/h[2]);old=np.maximum(np.log2((1+d[3]*z)/h[3]),np.log2((1+d[4]*z)/h[4]))
    value=np.maximum(lower*l2+(count-lower)*old,upper*l2+(count-upper)*old)-degree*x
    v2=d[2]*z/(1+d[2]*z)**2
    vo=np.minimum(d[3]*z/(1+d[3]*z)**2,d[4]*z/(1+d[4]*z)**2)
    variance=np.minimum(lower*v2+(count-lower)*vo,upper*v2+(count-upper)*vo)
    mass=np.minimum(1.,np.sqrt(math.pi/(8*variance)))
    for theta in (.5,.75,1.,1.25):
        c=1-theta*theta/12
        mass=np.minimum(mass,1/np.sqrt(2*math.pi*variance*c)+np.exp(-variance*theta*theta*c/2))
    return value+np.log2(mass)


def allowed_logs(a,l):
    result=np.full(min(a,l)+1,-np.inf)
    for k,J,logJ in occupation(l):
        for v0 in range(max(0,k-a),min(k,512-a)+1):
            value=logJ-logchoose(512,k)+logchoose(512-a,v0)
            value+=v0*math.log2(fresh.D[0])+(l-k)*math.log2(fresh.H[1])
            result[k-v0]=np.logaddexp2(result[k-v0],value)
    return result-LOGZ*np.arange(len(result))


def allowed_tables(incoming,outgoing):
    """Same A_(a,l,u), via the exact inactive-coordinate recurrence.

    D0=H1 cancels the v0 dependence outside the binomial sum. Adding one
    inactive coordinate gives A_(a-1,l,u)=A_(a,l,u)+H1*A_(a,l,u+1).
    Keep logarithms; no tiny positive summands are explicitly dropped.
    """
    assert fresh.D[0]==fresh.H[1]
    lh=math.log2(fresh.H[1]);table=np.full((len(outgoing),513),-np.inf)
    for j,l in enumerate(outgoing):
        for k,J,logJ in occupation(l):
            table[j,k]=logJ-logchoose(512,k)+(l-k)*lh
    requested=set(incoming);result={}
    for a in range(512,min(incoming)-1,-1):
        if a in requested:
            result[a]={l:table[j,:min(a,l)+1].copy()-LOGZ*np.arange(min(a,l)+1)
                       for j,l in enumerate(outgoing)}
        if a:
            table=np.logaddexp2(table[:,:a],lh+table[:,1:a+1])
    return result


def build(counts):
    start=time.monotonic();sizes=[len(states(a)[0]) for a in counts]
    offsets=np.cumsum([0]+sizes);dimension=int(offsets[-1])
    memory_limit=(8<<30) if BIN==16 and len(counts)>4 else 768*(1<<20)
    assert 6*dimension*dimension*8+4*513*dimension*8+32*(1<<20)<memory_limit
    logs=np.full((dimension,dimension),-np.inf)
    all_allowed=allowed_tables(counts,counts)
    # This coefficient depends on a-r, not on a and r separately.
    # Reuse exactly the same evaluated array across incoming input counts.
    old_by_count={}
    for n in range(max(counts)+1):
        bins=[(s,min(n,s+BIN-1)) for s in range(0,n+1,BIN)]
        old_by_count[n]=old_coefficients(n,bins)+LOGZ*np.arange(n+1)[None,:]
    print(dict(stage='joint_old_coefficients',counts=len(old_by_count),
               seconds=time.monotonic()-start),flush=True)
    for i,a in enumerate(counts):
        boxes,source=states(a);allowed=all_allowed[a]
        for r in range(a+1):
            n=a-r;second_bins=[(s,min(n,s+BIN-1)) for s in range(0,n+1,BIN)]
            old=old_by_count[n]
            B=np.array([logchoose(r,v)+v*math.log2(fresh.D[1])-r*math.log2(fresh.H[1])+LOGZ*v for v in range(r+1)])
            bscale=float(np.max(B));B=np.exp2(B-bscale);assert np.all(B>0)
            for sj,(slo,shi) in enumerate(second_bins):
                escale=float(np.max(old[sj]));E=np.exp2(old[sj]-escale);assert np.all(E>0)
                incoming=int(offsets[i])+source[r//BIN,slo//BIN]
                # The convolution is shared by every outgoing input count.
                convolutions=[]
                for vstart in range(0,r+1,BIN):
                    conv=np.convolve(B[vstart:min(r+1,vstart+BIN)],E)
                    assert np.all(conv>0),'Positive convolution underflow'
                    convolutions.append((vstart,np.log2(conv)+bscale+escale))
                for j,l in enumerate(counts):
                    destination=states(l)[1];row=np.full(sizes[j],-np.inf);A=allowed[l]
                    for vstart,conv in convolutions:
                        length=min(len(conv),len(A)-vstart)
                        if length<=0:continue
                        u=np.arange(vstart,vstart+length);selected=vstart//BIN
                        dest=destination[(l-u)//BIN,selected]
                        values=A[u]+conv[:length]
                        assert np.all(~np.isfinite(values[dest<0])),'Invalid outgoing joint bucket has mass'
                        keep=dest>=0
                        np.logaddexp2.at(row,dest[keep],values[keep])
                    sl=slice(offsets[j],offsets[j+1])
                    logs[incoming,sl]=np.maximum(logs[incoming,sl],row)
        print(dict(stage='joint_rows',input_count=a,rows=sizes[i],seconds=time.monotonic()-start),flush=True)
    scale=float(np.max(logs));matrix=np.exp2(logs-scale)
    assert np.all(matrix[np.isfinite(logs)]>0),'Matrix underflow'
    labels=np.repeat(counts,sizes)
    return matrix,scale,labels,offsets


def initial_values(counts,offsets,scale,x):
    values=np.zeros(int(offsets[-1]))
    for i,l in enumerate(counts):
        j=states(l)[1][l//BIN,0]
        values[offsets[i]+j]=2**(logchoose(512,l)+l*math.log2(fresh.H[1])-scale+x*(l-408))
    return values


def terminal_logs(counts):
    rewards=fresh.LOGB*np.array((fresh.LAMBDA,2*fresh.LAMBDA,3*fresh.LAMBDA,fresh.LAMBDA+fresh.TAU))-np.log2(fresh.H[1:])
    return [max(r*rewards[0]+s*rewards[1]+(a-r-s)*max(rewards[2:]) for r,s in vertices(a,box))
            for a in counts for box in states(a)[0]]


def describe_label(a,index):
    box=states(a)[0][index]
    return dict(input_count=a,fresh_count_bucket=box[:2],second_age_bucket=box[2:])


def main():
    global BIN
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    budget=1200 if '--full-band' in sys.argv else 90
    timer=threading.Timer(budget,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        BIN=16 if '--bin16' in sys.argv else 32;fresh.BIN=BIN;states.cache_clear()
        counts=list(range(392,425)) if '--full-band' in sys.argv else [392,407,408,424]
        result=regional.screen(counts,all_payments=True,kernel=sys.modules[__name__])
        result.update(allowed_input_counts=counts,joint_cohorts=True,seconds=time.monotonic()-start,files_written=0)
        print(result,flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

"""Two-update bounds retaining the intermediate third-age count.

One of the two input counts is the full width. FLOAT upper bounds and
exact toy identities only; no distance certificate. Serial Peach,180s,
<768MiB, no generated files.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_mixed_counts as mixed
import quartet_joint_regional as joint
import quartet_fresh_count_screen as fresh
from quartet_history_screen import occupation,logchoose
from quartet_window_transfer import kernel,gate_counts
from quartet_two_cohort import compositions,elementary,H,D


def entries(width,l):
    if width==512:return occupation(l)
    return [(k,J,math.log2(J)) for (a,k),J in gate_counts(width).items() if a==l]


def prepare_first(cell,width=512,bucket=16,h=None,d=None):
    h=fresh.H if h is None else np.asarray(h,dtype=float)
    d=fresh.D if d is None else np.asarray(d,dtype=float)
    a,ahi,rl,rh,sl,sh=cell;assert a==ahi
    n0=width-a;b3=d[3]/h[4];b4=d[4]/h[4]
    old=max(b3/h[3],b4/h[4])
    hs=np.array([1,h[1],h[2],1/old,1/old])
    ds=np.array([1,d[2],d[3],d[4],d[4]])
    E=mixed.coupled_coefficients(cell,width=a,h=hs,d=ds) if a else np.zeros(1)
    tilt=-2.;E=E+tilt*np.arange(a+1);es=max(E)
    E=mixed.positive_exp(E-es,'pair_first_polynomial')
    B=np.array([logchoose(n0,v)+v*math.log2(d[1])+tilt*v for v in range(n0+1)])
    prepared=[]
    for vstart in range(0,n0+1,bucket):
        b=B[vstart:vstart+bucket];bs=max(b)
        coeff=np.convolve(mixed.positive_exp(b-bs,'pair_first_binomial'),E)
        logs=np.log2(np.maximum(coeff,2**mixed.FLOOR))+bs+es
        prepared.append((vstart,logs))
    return prepared


@lru_cache(maxsize=128)
def first_geometry(width,l,h1):
    terms=entries(width,l)
    return (np.array([k for k,J,logJ in terms]),
            np.array([logJ-logchoose(width,k)+(l-k)*math.log2(h1)+2*k for k,J,logJ in terms]))


def full_first(cell,l,width=512,bucket=16,h=None,d=None,prepared=None):
    h=fresh.H if h is None else np.asarray(h,dtype=float)
    if prepared is None:prepared=prepare_first(cell,width,bucket,h,d)
    out={};ks,geometry=first_geometry(width,l,float(h[1]))
    for vstart,logs in prepared:
        valid=(ks>=vstart)&(ks-vstart<len(logs));k=ks[valid]
        values=geometry[valid]+logs[k-vstart];row=np.full(l//bucket+1,-np.inf)
        np.logaddexp2.at(row,(l-k)//bucket,values)
        for r in np.flatnonzero(np.isfinite(row)):out[int(r),vstart//bucket]=row[r]
    return out


def prepare_second(cell,width=512,h=None,d=None):
    h=fresh.H if h is None else np.asarray(h,dtype=float)
    d=fresh.D if d is None else np.asarray(d,dtype=float)
    a,ahi,rl,rh,sl,sh=cell;assert a==ahi
    b3=d[3]/h[4];b4=d[4]/h[4]
    ds=np.array([h[2],h[3],h[4]*b3,h[4]*b3*b4,h[4]*b4*b4])
    return mixed.coupled_coefficients(cell,width=a,h=h,d=ds) if a else np.zeros(1)


@lru_cache(maxsize=128)
def second_geometry(a,l,width,h2):
    A=np.full(a+1,-np.inf)
    for k,J,logJ in entries(width,l):
        for v0 in range(max(0,k-a),min(k,width-a)+1):
            u=k-v0
            value=logJ-logchoose(width,k)+logchoose(width-a,v0)+(l-u)*math.log2(h2)
            A[u]=np.logaddexp2(A[u],value)
    return A


def full_second(cell,l,width=512,bucket=16,h=None,d=None,prepared=None):
    h=fresh.H if h is None else np.asarray(h,dtype=float)
    E=prepare_second(cell,width,h,d) if prepared is None else prepared
    A=second_geometry(cell[0],l,width,float(h[2]));u=np.flatnonzero(np.isfinite(A))
    values=A[u]+E[u]+(width-l)*math.log2(h[1]);row=np.full(l//bucket+1,-np.inf)
    np.logaddexp2.at(row,(l-u)//bucket,values)
    return {((width-l)//bucket,int(s)):row[s] for s in np.flatnonzero(np.isfinite(row))}


def exact_formula(n,l,first,h=H,d=D):
    width=sum(n);phi=math.prod(hi**ni for hi,ni in zip(h,n));result=defaultdict(F)
    b3=F(d[3],h[4]);b4=F(d[4],h[4])
    if first:
        poly=elementary((n[1],n[2],n[3]+n[4]),(d[2],d[3],d[4]))
        constant=b3**n[3]*b4**n[4]/phi
        for k,J,logJ in entries(width,l):
            for v in range(max(0,k-sum(n[1:])),min(n[0],k)+1):
                value=F(J*math.comb(n[0],v),math.comb(width,k))*d[1]**v*poly[k-v]*h[1]**(l-k)*constant
                result[l-k,v]+=value
    else:
        ds=(h[2],h[3],h[4]*b3,h[4]*b3*b4,h[4]*b4*b4)
        poly=elementary(n[1:],ds[1:])
        for k,J,logJ in entries(width,l):
            for v in range(max(0,k-sum(n[1:])),min(n[0],k)+1):
                value=F(J*math.comb(n[0],v),math.comb(width,k))*ds[0]**v*poly[k-v]*h[2]**(l-k)*h[1]**(width-l)/phi
                result[width-l,l-k+v]+=value
    return dict(result)


def tests():
    exact=bounded=0
    for width in (4,8):
        histograms=list(compositions(width,5)) if width==4 else [(2,2,2,1,1),(0,0,0,0,8),(8,0,0,0,0)]
        cache={}
        for n in histograms:
            a=width-n[0];r=(n[1]//2)*2;s=(n[2]//2)*2
            cell=(a,a,r,min(a,r+1),s,min(a-r,s+1));phi=math.prod(h**v for h,v in zip(H,n))
            for l in range(width+1):
                for first in (False,True):
                    direct=defaultdict(F)
                    for (mid,l1,e1),p1 in kernel(n):
                        if l1!=(width if first else l):continue
                        for (out,l2,e2),p2 in kernel(mid):
                            if l2!=(l if first else width):continue
                            direct[out[1],out[2]]+=p1*p2*(1<<(e1+e2))*F(math.prod(h**v for h,v in zip(H,out)),phi)
                    assert dict(direct)==exact_formula(n,l,first);exact+=len(direct)
                    key=(cell,l,first)
                    if key not in cache:
                        fn=full_first if first else full_second
                        cache[key]=fn(cell,l,width,2,H,D)
                    actual=defaultdict(F)
                    for (r1,s1),v in direct.items():actual[r1//2,s1//2]+=v
                    for key,value in actual.items():
                        assert math.log2(value)<=cache[cell,l,first][key]+1e-10;bounded+=1
    print(dict(exact_full_pair_coefficients=exact,coefficients_vs_uniform_float_pair_bounds=bounded),flush=True)


def logsum(values):
    peak=max(values)
    return peak+math.log2(sum(mixed.positive_exp(np.array(values)-peak,'pair_diagnostic_sum')))


def diagnostic():
    joint.BIN=16;joint.states.cache_clear();counts=(312,408,512)
    base,scale,labels,offsets=joint.build(counts);where={a:i for i,a in enumerate(counts)}
    for a,rl,sl,l in ((312,64,16,312),(408,192,96,312),(408,192,96,408),(512,0,0,312)):
        i=where[a];boxes,lookup=joint.states(a);j=lookup[rl//16,sl//16]
        cell=(a,a,*boxes[j]);source=offsets[i]+j
        for first in (False,True):
            mid=512 if first else l;last=l if first else 512
            mi=where[mid];oi=where[last];ms=slice(offsets[mi],offsets[mi+1]);dest=slice(offsets[oi],offsets[oi+1])
            left=base[source,ms];right=base[ms,dest];ls=max(left);rs=np.max(right)
            row=(left/ls)@(right/rs)
            support=(left>0).astype(int)@(right>0).astype(int)
            assert np.all(row[support>0]>0),'Positive pair-product underflow'
            with np.errstate(divide='ignore'):old=2*scale+math.log2(ls)+math.log2(rs)+np.log2(row)
            coupled=(full_first if first else full_second)(cell,l)
            new=np.full(len(old),-np.inf);dst=joint.states(last)[1]
            for key,value in coupled.items():
                if key[0]<dst.shape[0] and key[1]<dst.shape[1] and dst[key]>=0:new[dst[key]]=value
                else:raise AssertionError('Invalid destination in coupled pair')
            best=np.minimum(old,new)
            print(dict(source=cell,input_counts=[mid,last],ordinary_pair_log2=logsum(old),
                       coupled_pair_log2=logsum(new),intersected_pair_log2=logsum(best),
                       strictly_improved_edges=int(np.sum(best<old-1e-8)),
                       local_only=True,directed_certificate=False),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests();diagnostic();print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

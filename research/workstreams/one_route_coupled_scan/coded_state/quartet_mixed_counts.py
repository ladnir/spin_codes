"""Coarse incoming count cells, exact outgoing count labels.

FLOAT proof diagnostics only. No encoder change or distance certificate.
Serial Peach,300s,2GiB conservative array budget, no generated files.
"""
import math
import os
from pathlib import Path
import shutil
import sys
import threading
import time
import numpy as np
import quartet_joint_regional as joint
import quartet_fresh_count_screen as fresh
import quartet_fresh_regional as regional
from quartet_history_screen import occupation,logchoose

WIDTH=512
BIN=64
LOGZ=joint.LOGZ
FLOOR=-1000.
FLOOR_COUNTS={}


def positive_exp(logs,kind):
    """Raise tiny positive terms to a positive floor; preserve exact zeros.

    This is a conservative representation change, not tail deletion.
    Ordinary floating evaluation errors still need directed verification.
    """
    assert np.all(np.isfinite(logs)|np.isneginf(logs)),'Nonfinite positive coefficient'
    finite=np.isfinite(logs)
    FLOOR_COUNTS[kind]=FLOOR_COUNTS.get(kind,0)+int(np.sum(finite&(logs<FLOOR)))
    return np.where(finite,np.exp2(np.maximum(logs,FLOOR)),0.)


def old_interval_coefficients(nlo,nhi,bins):
    """Uniform over n in [nlo,nhi], s in its bin, s<=n and all old splits."""
    if nhi==0:return np.zeros((1,1))
    d=fresh.D;h=fresh.H
    sl=np.array([lo for lo,hi in bins])[:,None]
    sh=np.array([hi for lo,hi in bins])[:,None]
    degree=np.arange(nhi+1)[None,:]
    low=np.maximum(np.maximum(nlo,degree),sl)+np.zeros_like(sh)
    high=np.full_like(low,nhi)
    corner=np.minimum(high,np.maximum(low,sh))
    ns=np.stack((low,low,high,high,corner))
    ss=np.stack((np.broadcast_to(sl,low.shape),np.minimum(sh,low),
                 np.broadcast_to(sl,low.shape),np.minimum(sh,high),np.minimum(sh,corner)))
    left=np.full(low.shape,-80.);right=np.full_like(left,80.)
    def evaluate(x):
        z=np.exp2(x)
        l2=np.log2((1+d[2]*z)/h[2])
        l3=np.log2((1+d[3]*z)/h[3]);l4=np.log2((1+d[4]*z)/h[4])
        old=np.maximum(l3,l4);selected=np.where(l3>=l4,d[3],d[4])
        values=ss*l2+(ns-ss)*old
        best=np.argmax(values,axis=0)[None,:,:]
        n=np.take_along_axis(ns,best,axis=0)[0]
        s=np.take_along_axis(ss,best,axis=0)[0]
        derivative=s*d[2]*z/(1+d[2]*z)+(n-s)*selected*z/(1+selected*z)-degree
        return z,np.max(values,axis=0)-degree*x,derivative
    for _ in range(48):
        x=(left+right)/2;z,value,derivative=evaluate(x)
        left=np.where(derivative<0,x,left);right=np.where(derivative>=0,x,right)
    z,value,_=evaluate((left+right)/2)
    v2=d[2]*z/(1+d[2]*z)**2
    vo=np.minimum(d[3]*z/(1+d[3]*z)**2,d[4]*z/(1+d[4]*z)**2)
    variance=np.min(ss*v2+(ns-ss)*vo,axis=0)
    with np.errstate(divide='ignore'):
        mass=np.minimum(1.,np.sqrt(math.pi/(8*variance)))
        for theta in (.5,.75,1.,1.25):
            c=1-theta*theta/12
            mass=np.minimum(mass,1/np.sqrt(2*math.pi*variance*c)+np.exp(-variance*theta*theta*c/2))
    return value+np.log2(mass)


def geometry(intervals):
    """Uniform A_(a,l,u) across each interval; all513 outgoing counts."""
    lh=math.log2(fresh.H[1]);table=np.full((513,513),-np.inf)
    for l in range(513):
        for k,J,logJ in occupation(l):table[l,k]=logJ-logchoose(512,k)+(l-k)*lh
    endpoints={lo for lo,hi in intervals};saved={};diagonal=np.full_like(table,-np.inf)
    for a in range(512,-1,-1):
        diagonal[:,a]=table[:,a]
        if a in endpoints:saved[a]=table.copy()
        if a:table=np.logaddexp2(table[:,:a],lh+table[:,1:a+1])
    result=[]
    for lo,hi in intervals:
        A=diagonal[:,:hi+1].copy();A[:,:lo+1]=saved[lo]
        A-=LOGZ*np.arange(hi+1)[None,:]
        result.append(A)
    return result


def layout(intervals):
    cells=[];sources=[];group_for_count=np.full(513,-1,dtype=int)
    for g,(lo,hi) in enumerate(intervals):
        assert np.all(group_for_count[lo:hi+1]<0)
        group_for_count[lo:hi+1]=g
        lookup=np.full((hi//BIN+1,hi//BIN+1),-1,dtype=int)
        for r in range(0,hi+1,BIN):
            for s in range(0,hi-r+1,BIN):
                lookup[r//BIN,s//BIN]=len(cells)
                cells.append((lo,hi,r,min(hi,r+BIN-1),s,min(hi-r,s+BIN-1)))
        sources.append(lookup)
    assert np.all(group_for_count>=0)
    target_map=[];labels=[];edge_lookup=np.full((513,513//BIN+1,513//BIN+1),-1,dtype=int)
    for l in range(513):
        group=group_for_count[l]
        for r in range(0,l+1,BIN):
            for s in range(0,l-r+1,BIN):
                edge_lookup[l,r//BIN,s//BIN]=len(labels)
                target_map.append(sources[group][r//BIN,s//BIN]);labels.append(l)
    return cells,sources,edge_lookup,np.array(labels),np.array(target_map)


def cell_vertices(cell):
    alo,ahi,rl,rh,sl,sh=cell;points=set()
    for a in {alo,ahi,rl+sl,rl+sh,rh+sl,rh+sh}:
        if alo<=a<=ahi:
            points.update((a,r,s) for r,s in joint.vertices(a,(rl,rh,sl,sh)))
    assert points
    return sorted(points)


def coupled_coefficients(cell,width=512,d=None,h=None):
    """Total-coefficient cap retaining n0+n1+n2+n3+n4=width."""
    vertices=np.array(cell_vertices(cell),dtype=float)
    a=vertices[:,0,None];r=vertices[:,1,None];s=vertices[:,2,None]
    counts=np.stack((width-a,r,s,a-r-s))
    d=fresh.D if d is None else np.asarray(d)
    h=fresh.H if h is None else np.asarray(h)
    degree=np.arange(width+1,dtype=float)[None,:]
    left=np.full((1,width+1),-80.);right=np.full_like(left,80.)
    def evaluate(x):
        z=np.exp2(x)
        logs=np.log2((1+d[:,None]*z)/h[:,None])
        old=np.maximum(logs[3],logs[4]);chosen=np.where(logs[3]>=logs[4],d[3],d[4])
        values=counts[0]*logs[0]+counts[1]*logs[1]+counts[2]*logs[2]+counts[3]*old
        best=np.argmax(values,axis=0)[None,:]
        n=np.take_along_axis(counts,np.broadcast_to(best,(4,1,width+1)),axis=1)
        ps=d[:3,None]*z/(1+d[:3,None]*z)
        derivative=np.sum(n[:3,0]*ps,axis=0)+n[3,0]*chosen*z/(1+chosen*z)-degree
        return z,np.max(values,axis=0)[None,:]-degree*x,derivative
    for _ in range(48):
        x=(left+right)/2;z,value,derivative=evaluate(x)
        left=np.where(derivative<0,x,left);right=np.where(derivative>=0,x,right)
    z,value,_=evaluate((left+right)/2)
    vs=d[:,None]*z/(1+d[:,None]*z)**2
    variance=np.min(counts[0]*vs[0]+counts[1]*vs[1]+counts[2]*vs[2]+counts[3]*np.minimum(vs[3],vs[4]),axis=0)
    mass=np.minimum(1.,np.sqrt(math.pi/(8*variance)))
    for theta in (.5,.75,1.,1.25):
        c=1-theta*theta/12
        mass=np.minimum(mass,1/np.sqrt(2*math.pi*variance*c)+np.exp(-variance*theta*theta*c/2))
    return value[0]+np.log2(mass)


def build(intervals,only_cell=None,only_interval=None):
    start=time.monotonic();cells,sources,edge_ids,labels,target_map=layout(intervals)
    D=len(cells);V=len(labels)
    # Rectangular matrices and the largest16-step shifted counter workspace.
    assert (6*D*V+2*8193*(D+V))*8+128*(1<<20)<2*(1<<30)
    logs=np.full((D,V),-np.inf);allowed=geometry(intervals)
    lgrid=np.arange(513)[:,None];rgrid=np.arange(513)[None,:];u=lgrid-rgrid
    padded=((513+BIN-1)//BIN)*BIN
    for g,(alo,ahi) in enumerate(intervals):
        if only_interval is not None and (alo,ahi)!=only_interval:continue
        if only_cell is not None and (alo,ahi)!=only_cell[:2]:continue
        valid=(u>=0)&(u<=ahi)
        A=np.full((513,513),-np.inf)
        li,ri=np.where(valid);A[li,ri]=allowed[g][li,u[li,ri]]
        for r in range(ahi+1):
            if only_cell is not None and not only_cell[2]<=r<=only_cell[3]:continue
            nlo=max(0,alo-r);nhi=ahi-r
            bins=[(s,min(nhi,s+BIN-1)) for s in range(0,nhi+1,BIN)]
            old=old_interval_coefficients(nlo,nhi,bins)+LOGZ*np.arange(nhi+1)[None,:]
            B=np.array([logchoose(r,v)+v*math.log2(fresh.D[1])-r*math.log2(fresh.H[1])+LOGZ*v for v in range(r+1)])
            bs=float(np.max(B));B=positive_exp(B-bs,'binomial')
            for sj,(slo,shi) in enumerate(bins):
                es=float(np.max(old[sj]));E=positive_exp(old[sj]-es,'old_coefficient')
                incoming=sources[g][r//BIN,slo//BIN]
                if only_cell is not None and cells[incoming]!=only_cell:continue
                for vstart in range(0,r+1,BIN):
                    conv=np.convolve(B[vstart:min(r+1,vstart+BIN)],E)
                    FLOOR_COUNTS['convolution']=FLOOR_COUNTS.get('convolution',0)+int(np.sum(conv<2**FLOOR))
                    conv=np.log2(np.maximum(conv,2**FLOOR))+bs+es
                    index=u-vstart;keep=valid&(index>=0)&(index<len(conv))
                    values=np.full((513,padded),-np.inf)
                    li,ri=np.where(keep);values[li,ri]=A[li,ri]+conv[index[li,ri]]
                    peaks=np.max(values,axis=1);safe=np.where(np.isfinite(peaks),peaks,0.)
                    scaled=positive_exp(values-safe[:,None],'bucket_sum')
                    sums=scaled.reshape(513,-1,BIN).sum(axis=2)
                    with np.errstate(divide='ignore'):block=np.log2(sums)+safe[:,None]
                    ids=edge_ids[:,:,vstart//BIN]
                    assert np.all(~np.isfinite(block[ids<0])),'Mass at invalid destination'
                    take=ids>=0
                    logs[incoming,ids[take]]=np.maximum(logs[incoming,ids[take]],block[take])
        print(dict(stage='coarse_count_rows',incoming_count_interval=[alo,ahi],
                   source_states=D,edge_labels=V,seconds=time.monotonic()-start),flush=True)
    scale=float(np.max(logs));base=positive_exp(logs-scale,'matrix')
    return base,scale,cells,edge_ids,labels,target_map


def row_diagnostic(intervals):
    cell=(320,383,0,63,0,63)
    base,scale,cells,edge_ids,labels,target_map=build(intervals,only_cell=cell)
    row=base[cells.index(cell)];coeff=coupled_coefficients(cell)
    results=[]
    for l in (320,384,392,408,424,448,512):
        terms=[logJ-logchoose(512,k)+(l-k)*math.log2(fresh.H[1])+coeff[k]
               for k,J,logJ in occupation(l)]
        peak=max(terms);cap=peak+math.log2(sum(2**(t-peak) for t in terms))
        uncoupled=scale+math.log2(sum(row[labels==l]))
        results.append(dict(input_count=l,float_log2_uncoupled_row=uncoupled,
                            float_log2_coupled_total_cap=cap,saving=uncoupled-cap))
    return dict(incoming_cell=cell,rows=results,local_total_caps_only=True,
                all_count_composition_done=False,directed_certificate=False)


def terminal_logs(cells):
    rewards=fresh.LOGB*np.array((fresh.LAMBDA,2*fresh.LAMBDA,3*fresh.LAMBDA,fresh.LAMBDA+fresh.TAU))-np.log2(fresh.H[1:])
    result=[]
    for alo,ahi,rl,rh,sl,sh in cells:
        candidates={alo,ahi,rl+sl,rl+sh,rh+sl,rh+sh};values=[]
        for a in candidates:
            if alo<=a<=ahi:
                values.extend(r*rewards[0]+s*rewards[1]+(a-r-s)*max(rewards[2:])
                              for r,s in joint.vertices(a,(rl,rh,sl,sh)))
        result.append(max(values))
    return np.array(result)


def screen(intervals):
    base,scale,cells,edge_ids,labels,target_map=build(intervals)
    order=np.argsort(target_map);sorted_targets=target_map[order]
    starts=np.flatnonzero(np.r_[True,sorted_targets[1:]!=sorted_targets[:-1]])
    assert np.array_equal(sorted_targets[starts],np.arange(len(cells)))
    def propose(x):
        with np.errstate(divide='ignore'):logs=np.log2(base)+x*(labels-408)[None,:]
        shift=float(np.max(logs))
        tilted=positive_exp(logs-shift,'tilted_edges')
        # Collapse edges only to propose a cell potential. The regional
        # counter below still uses every exact edge label separately.
        matrix=np.add.reduceat(tilted[:,order],starts,axis=1)
        v=np.ones(len(cells))
        for _ in range(512):
            after=matrix@v;after/=np.max(after);after=np.maximum(after,2**-400)
            residual=float(np.max(np.abs(np.log2(after/v))));v=after
            if residual<1e-10:break
        rho=float(np.max((matrix@v)/v))*(1+2**-30)
        collapsed=matrix*v[None,:]/(rho*v[:,None]);left=np.ones(len(cells))/len(cells)
        for _ in range(512):
            after=left@collapsed;after/=sum(after)
            error=float(np.max(np.abs(after-left)));left=after
            if error<1e-12:break
        edge=((left/v)@tilted)*v[target_map]/rho
        mean=float(edge@labels)/float(sum(edge))
        return mean,tilted,v,rho,residual,shift
    lo=-2.;hi=2.
    for _ in range(14):
        x=(lo+hi)/2;mean,tilted,v,rho,residual,shift=propose(x)
        if mean<408:lo=x
        else:hi=x
    x=(lo+hi)/2;mean,tilted,v,rho,residual,shift=propose(x)
    scale+=shift
    with np.errstate(divide='ignore'):
        P=positive_exp(np.log2(tilted)+np.log2(v[target_map])[None,:]-np.log2(v)[:,None]-math.log2(rho),'normalized_edges')
    assert np.max(P.sum(axis=1))<=1+1e-12
    initial_logs=np.full(len(labels),-np.inf)
    for l in range(513):
        e=edge_ids[l,l//BIN,0]
        initial_logs[e]=logchoose(512,l)+l*math.log2(fresh.H[1])-scale+x*(l-408)
    initial=positive_exp(initial_logs,'initial')
    initial_v=float(initial@v[target_map])/rho
    first=initial*v[target_map]/(rho*initial_v)
    P=np.vstack((P,first));terminal=terminal_logs(cells)-np.log2(v)
    rows=regional.coefficient_rows(P,labels,all_totals=True,target_map=target_map,positive_floor=2**FLOOR)
    qs=np.array([s for s in rows if 0<=436+s<=2100])
    quota=logchoose(6528,436)-logchoose(8192,2100)
    with np.errstate(divide='ignore'):logrows=np.log2(np.array([rows[s] for s in qs]))
    logrows+=np.array([logchoose(6528+s,436+s)-logchoose(6528,436) for s in qs])[:,None]
    def norm(y):
        values=logrows+qs[:,None]*(y-x);peaks=np.max(values,axis=0)
        weights=positive_exp(values-peaks[None,:],'payment_sum');sums=weights.sum(axis=0)
        norms=peaks+np.log2(sums);i=int(np.argmax(norms))
        return float(norms[i]),44+256*float(qs@weights[:,i])/float(sums[i]),i
    lo=x-math.log2(6528/436)-1;hi=lo+2
    for _ in range(10):
        if norm(lo)[1]<=0<=norm(hi)[1]:break
        width=hi-lo
        if norm(lo)[1]>0:lo-=width
        if norm(hi)[1]<0:hi+=width
    assert norm(lo)[1]<=0<=norm(hi)[1]
    for _ in range(40):
        mid=(lo+hi)/2;_,derivative,_=norm(mid)
        if derivative<0:lo=mid
        else:hi=mid
    y=(lo+hi)/2;gamma,_,worst=norm(y)
    # A substochastic labelled kernel propagates the total positive error
    # by at most its l1 coefficient mass. Each upward counter floor adds
    # <=epsilon per total, at <=k*512+1 totals on update k.
    weight_logs=np.array([logchoose(6528+s,436+s)-logchoose(6528,436)+s*(y-x) for s in qs])
    floor_gamma=FLOOR+math.log2(16+136*512)+float(np.max(weight_logs))+math.log2(len(qs))
    value=math.log2(initial_v)+4096*(scale+math.log2(rho))+256*quota+44*y+256*gamma+max(terminal)
    value+=logchoose(8192,2100)-fresh.LOGB*165556
    return dict(float_log2_all_counts_upper=value,all_input_counts=[0,512],
                all_regional_payments=True,inactive_rows=2100,covered_count=111572,
                cells=len(cells),edge_labels=len(labels),bucket_width=BIN,
                marker_log2=x,payment_marker_log2=y,tilted_mean=mean,
                worst_boundary_cell='initial' if worst==len(cells) else cells[worst],
                counter_floor_gamma_gap_bits=gamma-floor_gamma,
                proposal_residual=residual,upward_positive_floor_counts=FLOOR_COUNTS,
                directed_certificate=False,whole_code_distance=False)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(300,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        intervals=[(a,min(512,a+63)) for a in range(0,513,64)]
        result=row_diagnostic(intervals) if '--row-cap' in sys.argv else screen(intervals)
        result.update(seconds=time.monotonic()-start,files_written=0)
        print(result,flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

"""All-count regional composition with coupled total-row caps.

FLOAT diagnostic, not a distance certificate. Serial Peach,1200s,2GiB,
no generated files. The exact cap LP is evaluated by its positive dual.
"""
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_mixed_counts as mixed
import quartet_fresh_count_screen as fresh
from quartet_history_screen import occupation,logchoose


def sorted_cap(beta,values,order,batch=64):
    """Support function for 0<=p<=beta, sum p<=1, for each source row.

    Greedy prefixes only propose a dual threshold. Evaluate the dual using
    positive differences; a wrong threshold would still give an upper bound.
    values and order describe descending continuation values, one per row.
    """
    rows,edges=values.shape;result=np.empty((rows,len(beta)))
    assert beta.shape[1]==edges and np.all(beta>=0) and np.all(values>=0)
    for start in range(0,rows,batch):
        f=values[start:start+batch];idx=order[start:start+batch]
        capacities=beta[:,idx]
        prefix=np.cumsum(capacities,axis=2)
        crossed=prefix>=1.;position=np.argmax(crossed,axis=2)
        theta=np.take_along_axis(f[None,:,:],position[:,:,None],axis=2)[:,:,0]
        theta=np.where(crossed[:,:,-1],theta,0.)
        dual=theta+np.sum(capacities*np.maximum(f[None,:,:]-theta[:,:,None],0.),axis=2)
        result[start:start+len(f)]=np.minimum(dual,f[None,:,0]).T
    return result


def cap_values(beta,values):
    order=np.argsort(-values,axis=1)
    return sorted_cap(beta,np.take_along_axis(values,order,axis=1),order)


def make_parts(base,scale,cells,labels,target_map):
    start=time.monotonic();coeff=np.array([mixed.coupled_coefficients(c) for c in cells])
    part=[]
    for l in range(513):
        take=np.flatnonzero(labels==l);m=base[:,take]
        ks=np.array([k for k,J,logJ in occupation(l)])
        geom=np.array([logJ-logchoose(512,k)+(l-k)*math.log2(fresh.H[1]) for k,J,logJ in occupation(l)])
        logs=coeff[:,ks]+geom[None,:]
        peaks=np.max(logs,axis=1)
        coupled=peaks+np.log2(mixed.positive_exp(logs-peaks[:,None],'cap_sum').sum(axis=1))
        sums=np.sum(m,axis=1);assert np.all(sums>0)
        cap=np.minimum(coupled,scale+np.log2(sums))
        with np.errstate(divide='ignore'):ratios=scale+np.log2(m)-cap[:,None]
        beta=np.minimum(1.,mixed.positive_exp(np.minimum(ratios,0.),'cap_ratio')*(1+2**-30))
        part.append(dict(l=l,targets=target_map[take],beta=beta,cap=cap))
    print(dict(stage='coupled_total_caps',cells=len(cells),seconds=time.monotonic()-start),flush=True)
    return part


def logsum_rows(logs):
    peaks=np.max(logs,axis=1)
    return peaks+np.log2(mixed.positive_exp(logs-peaks[:,None],'positive_sum').sum(axis=1))


def total_map(parts,logv,x):
    logs=[];cache={}
    for p in parts:
        targets=p['targets'];key=tuple(targets)
        if key not in cache:
            peak=float(np.max(logv[targets]));f=mixed.positive_exp(logv[targets]-peak,'relative_potential')[None,:]
            order=np.argsort(-f,axis=1)
            cache[key]=(peak,np.take_along_axis(f,order,axis=1),order)
        peak,f,order=cache[key]
        value=sorted_cap(p['beta'],f,order)[0]
        assert np.all(value>0)
        logs.append(p['cap']+x*(p['l']-408)+peak+np.log2(value))
    return logsum_rows(np.array(logs).T)


def propose(parts,x,iterations=48):
    logv=np.zeros(len(parts[0]['cap']))
    for _ in range(iterations):
        after=total_map(parts,logv,x);after-=np.max(after)
        after=np.maximum(after,-400.)
        residual=float(np.max(np.abs(after-logv)))
        logv=(logv+after)/2
        if residual<1e-9:break
    logrho=float(np.max(total_map(parts,logv,x)-logv))+math.log2(1+2**-30)
    return logrho,logv,residual


def coefficients(parts,logv,logrho,x,steps=16,width=512,center=408,terminal=None):
    """Exact count labels; upper cap operator at every backward transition."""
    state=np.ones((1,len(logv))) if terminal is None else np.asarray(terminal)[None,:].copy()
    start=time.monotonic()
    groups={tuple(p['targets']) for p in parts}
    assert (4*(steps*width+1)*(len(logv)+sum(len(g) for g in groups))
            +10*len(logv)*64*max(len(g) for g in groups))*8<512*(1<<20)
    for step in range(1,steps+1):
        after=np.zeros((len(state)+width,len(logv)));cache={}
        for p in parts:
            l=p['l'];targets=p['targets'];key=tuple(targets)
            if key not in cache:
                peak=float(np.max(logv[targets]))
                with np.errstate(divide='ignore'):logs=np.log2(state[:,targets])+logv[targets][None,:]-peak
                f=mixed.positive_exp(logs,'cap_continuation')
                order=np.argsort(-f,axis=1)
                cache[key]=(peak,np.take_along_axis(f,order,axis=1),order)
            peak,f,order=cache[key]
            value=sorted_cap(p['beta'],f,order)
            with np.errstate(divide='ignore'):logs=np.log2(value)+p['cap'][None,:]+x*(l-center)+peak-logrho-logv[None,:]
            after[l:l+len(state)]+=mixed.positive_exp(logs,'capped_counter')
        state=np.maximum(after,2**mixed.FLOOR)
        assert np.all(np.isfinite(state)) and np.max(state)<=1+1e-7
        print(dict(stage='capped_region',updates=step,totals=len(state),seconds=time.monotonic()-start),flush=True)
    return state


def screen():
    intervals=[(a,min(512,a+63)) for a in range(0,513,64)]
    base,scale,cells,edge_ids,labels,target_map=mixed.build(intervals)
    parts=make_parts(base,scale,cells,labels,target_map)
    # The search only proposes proof weights. Every returned candidate is
    # checked with a maximum ratio; no convergence assumption is used.
    lo=-1.;hi=.5;ratio=(math.sqrt(5)-1)/2;cache={}
    def trial(x):
        if x not in cache:
            cache[x]=propose(parts,x)
            print(dict(stage='capped_potential',marker_log2=x,log2_growth=cache[x][0],residual=cache[x][2]),flush=True)
        return cache[x][0]
    left=hi-ratio*(hi-lo);right=lo+ratio*(hi-lo)
    for _ in range(7):
        if trial(left)<trial(right):hi=right;right=left;left=hi-ratio*(hi-lo)
        else:lo=left;left=right;right=lo+ratio*(hi-lo)
    x=min(cache,key=lambda z:cache[z][0]);logrho,logv,residual=cache[x]
    initial_logs=[]
    for p in parts:
        l=p['l'];edge=edge_ids[l,l//mixed.BIN,0];target=target_map[edge]
        initial_logs.append(logchoose(512,l)+l*math.log2(fresh.H[1])+x*(l-408)+logv[target])
    log_initial_v=float(logsum_rows(np.array(initial_logs)[None,:])[0])-logrho
    terminal=mixed.terminal_logs(cells)-logv
    for p in parts:
        l=p['l'];target=target_map[edge_ids[l,l//mixed.BIN,0]]
        row=(p['targets']==target).astype(float)
        assert sum(row)==1
        p['beta']=np.vstack((p['beta'],row))
        p['cap']=np.append(p['cap'],logchoose(512,l)+l*math.log2(fresh.H[1]))
    logv=np.append(logv,log_initial_v)
    rows=coefficients(parts,logv,logrho,x)
    qs=np.arange(-436,1665)
    quota=logchoose(6528,436)-logchoose(8192,2100)
    logrows=np.log2(rows[6528+qs])
    logrows+=np.array([logchoose(6528+s,436+s)-logchoose(6528,436) for s in qs])[:,None]
    def norm(y):
        values=logrows+qs[:,None]*(y-x);peaks=np.max(values,axis=0)
        weights=mixed.positive_exp(values-peaks[None,:],'capped_payment_sum');sums=weights.sum(axis=0)
        norms=peaks+np.log2(sums);i=int(np.argmax(norms))
        return float(norms[i]),44+256*float(qs@weights[:,i])/float(sums[i]),i
    lo=x-5.;hi=x-2.
    for _ in range(10):
        if norm(lo)[1]<=0<=norm(hi)[1]:break
        width=hi-lo
        if norm(lo)[1]>0:lo-=width
        if norm(hi)[1]<0:hi+=width
    assert norm(lo)[1]<=0<=norm(hi)[1]
    for _ in range(40):
        mid=(lo+hi)/2
        if norm(mid)[1]<0:lo=mid
        else:hi=mid
    y=(lo+hi)/2;gamma,_,worst=norm(y)
    value=log_initial_v+4096*logrho+256*quota+44*y+256*gamma+float(max(terminal))
    value+=logchoose(8192,2100)-fresh.LOGB*165556
    return dict(float_log2_capped_all_counts_upper=value,all_input_counts=[0,512],
                all_regional_payments=True,inactive_rows=2100,covered_count=111572,
                cells=len(cells),edge_labels=len(labels),marker_log2=x,payment_marker_log2=y,
                potential_residual=residual,worst_cell='initial' if worst==len(cells) else cells[worst],
                directed_certificate=False,whole_code_distance=False,
                upward_positive_floor_counts=mixed.FLOOR_COUNTS)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(1200,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        result=screen();result.update(seconds=time.monotonic()-start,files_written=0)
        print(result,flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

"""Product-potential screen for the closed-set quartet kernel.

Valid uniform formula, FLOAT evaluation only. K20,r2100, all paid counts
and histories. Serial Peach180s,<256MiB,no files. Not a directed cert.
"""
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
from quartet_closed_transfer import TYPES,local,terminal
from quartet_history_screen import logchoose


def terms():
    index={t:i for i,t in enumerate(TYPES)};features=[];logs=[];starts=[];rows=[]
    for hist in itertools.combinations_with_replacement(TYPES,4):
        transitions=local(hist)
        if not transitions:continue
        starts.append(len(features));rows.append(hist)
        source=np.bincount([index[t] for t in hist],minlength=12)
        for (dest,unpaid,covered,activity,deaths),ways in transitions:
            delta=np.bincount([index[t] for t in dest],minlength=12)-source
            for c in range(covered+1):
                features.append(list(delta[1:])+[unpaid+c,(707*activity-1905*deaths-3072*c)/2560])
                logs.append(math.log2(ways*math.comb(covered,c)))
    X=np.array(features,dtype=float);logc=np.array(logs);starts=np.array(starts)
    assert X.nbytes+logc.nbytes<64*(1<<20)
    return X,logc,starts,rows


def screen():
    start=time.monotonic();X,logc,starts,rows=terms()
    counts=np.diff(np.r_[starts,len(X)]);group=np.repeat(np.arange(len(rows)),counts)
    accepted=[i for i,t in enumerate(TYPES) if terminal((t,))]
    cache={};evaluations=0
    def kernel(v):
        nonlocal evaluations
        key=v[:13].tobytes()
        if cache.get('key')==key:return cache['values'],cache['jac']
        t=logc+X@v[:13];peak=np.maximum.reduceat(t,starts)
        p=np.exp2(t-peak[group]);totals=np.add.reduceat(p,starts)
        values=peak+np.log2(totals)
        jac=np.add.reduceat(p[:,None]*X,starts,axis=0)/totals[:,None]
        cache.update(key=key,values=values,jac=jac);evaluations+=1
        return values,jac
    # v:11 noninitial type potentials, original-inactive marker x,log2(b).
    linear=np.zeros(13);linear[11]=-537600/(1<<21);linear[12]=-(713983/5-111127)/(1<<21)
    def smooth(v,temp):
        values,jac=kernel(v);peak=max(values);p=np.exp2((values-peak)/temp);p/=sum(p)
        bound=peak+temp*math.log2(sum(np.exp2((values-peak)/temp)))
        h=np.r_[0.,v[:11]];terminal=-h[accepted];tp=max(terminal)
        q=np.exp2((terminal-tp)/temp);q/=sum(q)
        tail=tp+temp*math.log2(sum(np.exp2((terminal-tp)/temp)))
        grad=linear+.25*(p@jac)
        for i,weight in zip(accepted,q):
            if i:grad[i-1]-=512/(1<<21)*weight
        return float(linear@v+.25*bound+512/(1<<21)*tail),grad
    lo=np.array([-40.]*11+[-12.,.001]);hi=np.array([40.]*11+[12.,12.])
    v=np.zeros(13);v[11]=-1.5;v[12]=2;iterations=0
    for temp in (.1,.02,.005,.001):
        inverse=np.eye(13);f,g=smooth(v,temp)
        for _ in range(150):
            direction=-inverse@g
            if g@direction>=0:inverse=np.eye(13);direction=-g
            direction/=max(1.,max(abs(direction))/4)
            rate=1.
            for line in range(32):
                proposal=np.clip(v+rate*direction,lo,hi);step=proposal-v
                nf,ng=smooth(proposal,temp)
                if nf<=f+1e-4*(g@step):break
                rate*=.5
            if nf>f+1e-4*(g@step) or max(abs(step))<1e-10:break
            change=ng-g;curvature=change@step
            if curvature>1e-12:
                rho=1/curvature;update=np.eye(13)-rho*np.outer(step,change)
                inverse=update@inverse@update.T+rho*np.outer(step,step)
            else:inverse=np.eye(13)
            v,f,g=proposal,nf,ng;iterations+=1
        print(dict(stage='closed_local_potential',temperature=temp,iterations=iterations,
                   smooth_objective=f,seconds=time.monotonic()-start),flush=True)
    values,jac=kernel(v);h=np.r_[0.,v[:11]]
    gamma=max(values);tail=max(-h[accepted]);rowset=logchoose(8192,2100)
    value=(1<<19)*gamma-537600*v[11]-(713983/5-111127)*v[12]+512*tail-255*rowset
    worst=int(np.argmax(values))
    print(dict(float_log2_upper=float(value),log2_b=float(v[12]),log2_inactive_marker=float(v[11]),
               quartet_log2_bound=float(gamma),terminal_log2_bound=float(tail),
               log2_type_potentials=[float(x) for x in h],worst_quartet=rows[worst],
               optimizer='bounded smoothed BFGS proposal; no convergence assumption',
               optimizer_iterations=iterations,evaluations=evaluations,
               local_histograms=len(rows),monomials=len(X),all_paid_counts=True,
               all_histories=True,inactive_rows=2100,directed_certificate=False,
               seconds=time.monotonic()-start),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:screen();print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()

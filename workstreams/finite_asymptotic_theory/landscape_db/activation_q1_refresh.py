"""Q1 transfer that retains exact uniform refresh after a zero-input epoch.

Classes: Z=zero, D=arbitrary nonzero, U=uniform nonzero,
L=nonzero density at most 1/(2**s-2). State output precedes update.
This is a new binary64 evaluator; existing producer receipts remain unchanged.
"""
import math
import numpy as np


def epoch_logs(t,s,spectrum,lambdas):
    if s<2 or sum(spectrum.values())!=(1<<s)-1 or any(w<=0 or w>t or n<=0 for w,n in spectrum.items()):
        raise ValueError('complete nonzero inner spectrum required')
    lam=np.asarray(lambdas,dtype=float)
    if lam.ndim!=1 or not np.all(np.isfinite(lam)&(lam>0)):
        raise ValueError('positive finite tilts required')
    m=(1<<s)-1; log_m=math.log(m); kappa=math.log1p(1/(m-1)); survive=math.log1p(-1/m)
    m0=np.full(lam.shape,-np.inf);m1=np.full_like(m0,-np.inf)
    for w,n in spectrum.items():
        mass=math.log(n)-log_m
        np.logaddexp(m0,mass-lam*w,out=m0)
        np.logaddexp(m1,mass+math.log(w/t)-lam*(w-1),out=m1)
        if w<t:
            np.logaddexp(m1,mass+math.log((t-w)/t)-lam*(w+1),out=m1)
    d=min(spectrum); worst0=-lam*d; worst1=-lam*(d-1)
    zero=np.full((len(lam),4,4),-np.inf);one=np.full_like(zero,-np.inf)
    zero[:,0,0]=0
    zero[:,1,2]=worst0;zero[:,2,2]=m0;zero[:,3,2]=np.minimum(worst0,kappa+m0)
    one[:,0,1]=-lam
    for state,moment in ((1,worst1),(2,m1),(3,np.minimum(worst1,kappa+m1))):
        one[:,state,0]=moment-log_m
        one[:,state,3]=moment+survive
    return zero,one


def matrix_product(a,b):
    return np.logaddexp(np.logaddexp(a[:,:,0,None]+b[:,None,0,:],a[:,:,1,None]+b[:,None,1,:]),
                       np.logaddexp(a[:,:,2,None]+b[:,None,2,:],a[:,:,3,None]+b[:,None,3,:]))


def region_logs(zero,one,epochs):
    if epochs<1:
        raise ValueError('positive epoch count required')
    rz=np.full_like(zero,-np.inf);rz[:,0,0]=rz[:,1,1]=rz[:,2,2]=rz[:,3,3]=0
    ra=np.full_like(one,-np.inf);bz,ba=zero,one
    while epochs:
        if epochs&1:
            ra=np.logaddexp(matrix_product(ra,bz),matrix_product(rz,ba));rz=matrix_product(rz,bz)
        epochs>>=1
        if epochs:
            ba=np.logaddexp(matrix_product(ba,bz),matrix_product(bz,ba));bz=matrix_product(bz,bz)
    return rz,ra


def vector_product(v,m,out):
    np.logaddexp(v[:,:,0,None]+m[:,None,0,:],v[:,:,1,None]+m[:,None,1,:],out=out)
    np.logaddexp(out,v[:,:,2,None]+m[:,None,2,:],out=out)
    np.logaddexp(out,v[:,:,3,None]+m[:,None,3,:],out=out)


def coefficient_logs(zero,one,epochs,block):
    rz,ra=region_logs(zero,one,epochs);ra-=math.log(epochs)
    shape=(len(rz),block+1,4)
    current=np.full(shape,-np.inf);updated=np.empty(shape);marked=np.empty(shape)
    current[:,0,0]=0
    for n in range(block):
        vector_product(current[:,:n+1],rz,updated[:,:n+1])
        vector_product(current[:,:n+1],ra,marked[:,:n+1])
        updated[:,n+1]=-np.inf
        np.logaddexp(updated[:,1:n+2],marked[:,:n+1],out=updated[:,1:n+2])
        current,updated=updated,current
    return np.logaddexp.reduce(current,axis=2)-np.array([math.log(math.comb(block,w)) for w in range(block+1)])


def screen(t,s,spectrum,block,rows,counts,log_tilts,refine=False):
    if rows<1 or rows%t:
        raise ValueError('whole epochs required')
    lam=np.exp(np.asarray(log_tilts));cutoff=block*rows//10
    moments=coefficient_logs(*epoch_logs(t,s,spectrum,lam),rows//t,block)
    values=np.minimum(0.,moments+cutoff*lam[:,None]);indices=np.argmin(values,axis=0)
    best=values[indices,np.arange(block+1)];witnesses=np.asarray(log_tilts)[indices]
    weights=sorted(counts);costs=np.array([math.log(rows)+math.log(counts[w]) for w in weights])
    terms=costs+best[weights];evaluated=np.asarray(log_tilts)
    if refine:
        relevant=[w for w,value in zip(weights,terms) if value>=max(terms)-30*math.log(2)]
        extra=np.unique(np.concatenate([witnesses[w]+np.arange(-30,31)/500 for w in relevant]))
        extra_lam=np.exp(extra)
        extra_values=np.minimum(0.,coefficient_logs(*epoch_logs(t,s,spectrum,extra_lam),rows//t,block)
                                +cutoff*extra_lam[:,None])
        indices=np.argmin(extra_values,axis=0);candidate=extra_values[indices,np.arange(block+1)]
        improve=candidate<best;best[improve]=candidate[improve];witnesses[improve]=extra[indices[improve]]
        terms=costs+best[weights];evaluated=np.unique(np.r_[evaluated,extra])
    winner=weights[int(np.argmax(terms))]
    return dict(margin_bits=-float(np.logaddexp.reduce(terms))/math.log(2),dominant_weight=winner,
                dominant_log_surprisal=float(witnesses[winner]),log_tilts_evaluated=evaluated)

"""PROPOSED state-only refresh: Y=X+A Q, Q'=T Q+C X.

Three envelopes retain zero mass, a flat nonzero density cap, and arbitrary
nonzero residual mass. The fresh T acts after emission. This is a different
construction from all native modes; FLOAT screens are not certificates.
"""
from functools import lru_cache
from math import comb
import numpy as np
from hd_bucket_feedback import profiles
from hd_systematic_mds_screen import optimize


@lru_cache(1)
def matrices(z,threshold=8,average=False):
    q=65536.;d=q**16-1
    items=profiles(16,64,4);occ=np.array([row for row,count in items]);degrees=occ.sum(axis=1)
    rho=np.minimum(4,occ).sum(axis=1)
    quiet=np.array([sum(comb(j,c)*z**(j-c)*(1-z)**c*q**(-c-min(4,j-c))
                        for c in range(j+1)) for j in range(17)])
    k=np.prod(quiet[occ],axis=1);t=(z+(1-z)/q)**degrees
    base=q**(-4)
    excess=np.array([sum(comb(n,c)*z**(n-c)*(1-z)**c*(q**(-c)-base)
                         for c in range(min(4,n+1))) for n in range(17)])
    running=np.ones(len(items));nonzero=np.zeros(len(items))
    for b in range(4):
        e=excess[16-occ[:,b]]
        nonzero=nonzero*base+running*e;running*=base+e
    v=t*nonzero/(1-q**(-16))
    if average:
        hyper=np.array([sum(comb(13,c)*comb(3,l-c)*z**(13-c)
                            for c in range(max(0,l-3),min(13,l)+1))/comb(16,l) for l in range(17)])
        live=t*hyper[occ].mean(axis=1)
    else:live=t*z**np.maximum(0,13-occ.max(axis=1))
    birth=np.maximum(0,t-k)
    density=d*q**(-rho.astype(float))*(z+(1-z)/q)**(degrees-rho)
    density[degrees==0]=0
    use_flat=rho>=threshold
    values=np.zeros((len(items),3,3))
    values[:,0,0]=k;values[:,0,1]=np.where(use_flat,density,0);values[:,0,2]=np.where(use_flat,0,birth)
    values[:,1,0]=np.where(degrees>0,v/d,0);values[:,1,1]=v
    values[:,2,0]=np.where(degrees>0,live/d,0);values[:,2,1]=live
    mantissa,exponent=np.frexp(values)
    rounded=np.ldexp(np.ceil(mantissa*8192)/8192,exponent)
    keys=np.column_stack((degrees,rounded.reshape((-1,9))))
    unique,inverse=np.unique(keys,axis=0,return_inverse=True)
    weights=np.zeros(len(unique));np.add.at(weights,inverse,np.array([count for row,count in items],dtype=float))
    with np.errstate(divide='ignore'):
        return unique[:,0].astype(int),np.log2(weights),np.log2(unique[:,1:].reshape((-1,3,3)))


def coefficients(z,alpha,state=16,width=64,threshold=8,average=False):
    if state!=16 or width!=64:raise ValueError('Only four RS[16,4] blocks are screened')
    degrees,weight,logs=matrices(z,threshold,average);result=np.zeros((129,3,3))
    for first in range(0,len(degrees),32):
        left=logs[first:first+32]
        degree=(degrees[first:first+32,None]+degrees[None,:]).reshape(-1)
        counts=weight[first:first+32,None]+weight[None,:]
        for r in range(3):
            for s in range(3):
                entry=np.full(counts.shape,-np.inf)
                for h in range(3):entry=np.logaddexp2(entry,left[:,r,h,None]+logs[None,:,h,s])
                result[:,r,s]+=np.bincount(degree,weights=np.exp2(counts+alpha*entry).reshape(-1),minlength=129)
    with np.errstate(divide='ignore'):return np.log2(result)


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows',type=int,nargs='+',default=(512,2048))
    p.add_argument('--threshold',type=int,default=8)
    p.add_argument('--average',action='store_true')
    p.add_argument('--iterations',type=int,default=45)
    p.add_argument('--regional',action='store_true')
    args=p.parse_args()
    for r in args.rows:
        source=lambda z,a,s,w:coefficients(z,a,s,w,args.threshold,args.average)
        result=optimize(r,16,args.iterations,args.regional,pair_source=source)
        result.update(construction='state-only-diagonal-Cauchy-refresh',state_mixer=False,feedback_refresh=False,threshold=args.threshold,average=args.average)
        print(json.dumps(result),flush=True)

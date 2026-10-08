"""FLOAT proposal: common refresh with global or alternating RS emission.

Keep the scattered rank16 feedback and diagonal-Cauchy common refresh.
Changing A changes the construction. No native timing or certificate is
claimed here; all timing references to mode36 concern its block emission.
"""
from functools import lru_cache
from math import comb
import numpy as np
from hd_bucket_feedback import profiles
from hd_scattered_feedback import maps as scattered_maps
from hd_systematic_mds import maps as global_maps
from hd_blockwise_shared_transfer import mm,rank
from hd_common_refresh_screen import matrices as block_matrices
from hd_systematic_mds_screen import optimize


@lru_cache(None)
def check_global_return_invertible():
    a,_=global_maps(16,64);_,c=scattered_maps();ca=mm(c,a,65536)
    if rank(tuple(tuple(x^int(i==j) for j,x in enumerate(row)) for i,row in enumerate(ca)),65536)!=16:
        raise ArithmeticError('Global return map is singular')


@lru_cache(1)
def global_matrices(z):
    check_global_return_invertible();q=65536.;d=q**16-1
    items=profiles(16,64,4);occ=np.array([row for row,count in items]);j=occ.sum(axis=1)
    rho=np.minimum(4,occ).sum(axis=1)
    quiet=np.array([sum(comb(l,c)*z**(l-c)*(1-z)**c*q**(-c-min(4,l-c))
                        for c in range(l+1)) for l in range(17)])
    k=np.prod(quiet[occ],axis=1);t=(z+(1-z)/q)**j
    inactive=np.zeros((65,17));flat=np.zeros(65)
    for n in range(65):
        for e in range(min(16,n+1)):
            weight=comb(n,e)*z**(n-e)*(1-z)**e
            flat[n]+=weight*(q**(-e)-q**(-16))/(1-q**(-16))
            for r in range(1,17):inactive[n,r]+=weight*(q**(-max(0,r+e-16))-q**(-r))
    v=t*flat[64-j]
    ret=np.minimum(inactive[64-j,rho]*(z+(1-z)/q)**(j-rho),
                   (1-q**(-rho.astype(float)))*z**np.maximum(0,49-j))
    ret=np.minimum.reduce((ret,np.maximum(0,t-k),d*v))
    values=np.stack((k,np.maximum(0,t-k),ret/d,v),axis=1)
    mantissa,exponent=np.frexp(values);rounded=np.ldexp(np.ceil(mantissa*8192)/8192,exponent)
    keys=np.column_stack((j,rounded));unique,inverse=np.unique(keys,axis=0,return_inverse=True)
    weights=np.zeros(len(unique));np.add.at(weights,inverse,np.array([count for row,count in items],dtype=float))
    with np.errstate(divide='ignore'):
        return unique[:,0].astype(int),np.log2(weights),np.log2(unique[:,1:].reshape((-1,2,2)))


def coefficients(z,alpha,state=16,width=64,alternate=False):
    if state!=16 or width!=64:raise ValueError('Only the s16, width64 proposal is supported')
    dl,wl,left=block_matrices(z,scattered=True,joint_rank=True) if alternate else global_matrices(z)
    dr,wr,right=global_matrices(z);result=np.zeros((129,2,2))
    for start in range(0,len(dl),32):
        ll=left[start:start+32];degrees=(dl[start:start+32,None]+dr[None,:]).reshape(-1)
        counts=wl[start:start+32,None]+wr[None,:]
        for r in range(2):
            for s in range(2):
                entry=np.logaddexp2(ll[:,r,0,None]+right[None,:,0,s],ll[:,r,1,None]+right[None,:,1,s])
                result[:,r,s]+=np.bincount(degrees,weights=np.exp2(counts+alpha*entry).reshape(-1),minlength=129)
    with np.errstate(divide='ignore'):return np.log2(result)


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows',type=int,nargs='+',default=(512,2048))
    p.add_argument('--alternate',action='store_true')
    p.add_argument('--iterations',type=int,default=50)
    p.add_argument('--regional',action='store_true')
    args=p.parse_args()
    for r in args.rows:
        source=lambda z,a,s,w:coefficients(z,a,s,w,args.alternate)
        result=optimize(r,16,args.iterations,args.regional,pair_source=source)
        result.update(construction='common-refresh-global-emission',state_mixer=False,alternate_emission=args.alternate)
        print(json.dumps(result),flush=True)

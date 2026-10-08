"""FLOAT active-state-block envelopes for cheaper block-MDS emission.

Component u has a uniform density cap on states with exactly u nonzero
blocks. Keeping these components avoids assigning a one-block quiet cap
to every nonzero state. Fractional powers retain each input profile, but
are currently applied packetwise, not after paired path composition.
"""
from functools import lru_cache
from math import comb,prod
import numpy as np
from hd_bucket_feedback import profiles
from hd_systematic_mds_screen import optimize


@lru_cache(1)
def local_matrices(z,state=16,width=64,h=4,q=65536.):
    b=state//h;length=width//b;d=q**state-1
    items=profiles(state,width,h)
    occ=np.array([row for row,count in items]);degrees=occ.sum(axis=1)
    counts=np.array([count for row,count in items],dtype=float)
    quiet=np.array([sum(comb(l,c)*z**(l-c)*(1-z)**c*q**(-c-min(h,l-c)) for c in range(l+1)) for l in range(length+1)])
    k=np.prod(quiet[occ],axis=1);t=(z+(1-z)/q)**degrees
    # Mark of a uniform nonzero local state on n inactive positions.
    g=np.array([sum(comb(n,c)*z**(n-c)*(1-z)**c*(q**(-c)-q**(-h)) for c in range(min(h,n+1)))/(1-q**(-h)) for n in range(length+1)])
    elementary=np.zeros((len(items),b+1));elementary[:,0]=1
    for block in range(b):
        v=g[length-occ[:,block]]
        for u in range(block+1,0,-1):elementary[:,u]+=v*elementary[:,u-1]
    sizes=np.array([comb(b,u)*(q**h-1)**u for u in range(b+1)])
    values=np.zeros((len(items),b+1,b+1))
    values[:,0,0]=k
    values[:,0,1:]=(t-k)[:,None]*(sizes[1:]/d)[None,:]
    rank=np.minimum(h,occ).sum(axis=1)
    for u in range(1,b+1):
        v=t*elementary[:,u]/comb(b,u)
        floor=np.maximum(0,length-h+1-occ[:,-u:]).sum(axis=1)
        gamma=np.minimum.reduce((k,q**(-rank.astype(float))*z**floor,sizes[u]*v))
        values[:,u,0]=np.where(degrees>0,v/d,0)
        values[:,u,1:]=np.where(degrees>0,v,(q-1)/q*v)[:,None]*(sizes[1:]/d)[None,:]
        values[:,u,u]+=gamma/q
    return degrees,counts,values


def coefficients(z,alpha,state,width,h=4,paired_average=False):
    degrees,counts,values=local_matrices(z,state,width,h)
    size=values.shape[1];single=np.zeros((width+1,size,size))
    if paired_average:
        # Conditional means retain every profile and its multiplicity.
        # Jensen bounds the fractional moment of a two-packet product by
        # the fractional power of these means. This is weaker than retaining
        # all profile pairs, but avoids their Cartesian product.
        for i in range(size):
            for j in range(size):single[:,i,j]=np.bincount(degrees,weights=counts*values[:,i,j],minlength=width+1)
        single/=np.array([comb(width,j) for j in range(width+1)])[:,None,None]
        from hd_systematic_mds_screen import coefficients as pair_matrices
        with np.errstate(divide='ignore'):return pair_matrices(np.log2(single),alpha)
    for i in range(size):
        for j in range(size):single[:,i,j]=np.bincount(degrees,weights=counts*values[:,i,j]**alpha,minlength=width+1)
    paired=np.zeros((2*width+1,size,size))
    for i in range(size):
        for k in range(size):
            for j in range(size):paired[:,i,j]+=np.convolve(single[:,i,k],single[:,k,j])
    with np.errstate(divide='ignore'):return np.log2(paired)


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows',type=int,nargs='+',default=(256,512,2048))
    p.add_argument('--dimension',type=int,default=4)
    p.add_argument('--iterations',type=int,default=60)
    p.add_argument('--regional',action='store_true')
    p.add_argument('--paired-average',action='store_true')
    args=p.parse_args()
    for r in args.rows:
        source=lambda z,a,s,w:coefficients(z,a,s,w,args.dimension,args.paired_average)
        result=optimize(r,16,args.iterations,args.regional,pair_source=source)
        result.update(construction='block-emission-active-state-count',local_dimension=args.dimension,packetwise_fraction=not args.paired_average)
        print(json.dumps(result),flush=True)

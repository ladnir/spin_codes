"""Exact GF256 principal-nullity scan for the local return map I+A A^T.

The determinant lemma reduces every principal submatrix to I+C_D A_D.
This computes all 2^16 subsets, not sampled minors. No distance certificate.
"""
from functools import lru_cache
import numpy as np
from hd_bucket_feedback import maps
from hd_blockwise_shared_transfer import multiply


@lru_cache(None)
def nullities(h=4,beta=1):
    if h!=4 or not 0<beta<256:raise ValueError('Only the sixteen-position local map is scanned')
    _,full_c=maps(local_rank=h);c=np.array([row[:16] for row in full_c[:h]],dtype=np.uint16)
    table=np.array([[multiply(i,j,65536) for j in range(256)] for i in range(256)],dtype=np.uint16)
    inv=np.zeros(256,dtype=np.uint16)
    for i in range(1,256):inv[i]=int(np.flatnonzero(table[i]==1)[0])
    gram=np.zeros((1<<16,h,h),dtype=np.uint16)
    for i in range(16):
        leaf=table[c[:,i,None],c[None,:,i]]
        gram[1<<i:1<<(i+1)]=gram[:1<<i]^leaf
    gram=table[beta,gram]^np.eye(h,dtype=np.uint16)[None,:,:]
    count=len(gram);indices=np.arange(count);ranks=np.zeros(count,dtype=np.int64)
    for column in range(h):
        eligible=(np.arange(h)[None,:]>=ranks[:,None])&(gram[:,:,column]!=0)
        found=eligible.any(axis=1);take=np.flatnonzero(found)
        pivot=np.argmax(eligible[take],axis=1);row=ranks[take]
        old=gram[take,row,:].copy();gram[take,row,:]=gram[take,pivot,:];gram[take,pivot,:]=old
        normalized=table[gram[take,row,:],inv[gram[take,row,column]][:,None]]
        factors=gram[take,:,column].copy();factors[np.arange(len(take)),row]=0
        gram[take]^=table[factors[:,:,None],normalized[:,None,:]]
        gram[take,row,:]=normalized;ranks[take]+=1
    return h-ranks


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--beta',type=int,nargs='+',default=(1,2,3))
    args=p.parse_args();sizes=np.array([i.bit_count() for i in range(1<<16)])
    for beta in args.beta:
        values=nullities(beta=beta)
        print(json.dumps(dict(beta=beta,maximum_by_size=[int(values[sizes==c].max()) for c in range(17)],singular=int(np.count_nonzero(values)),worst=int(values.max()))),flush=True)

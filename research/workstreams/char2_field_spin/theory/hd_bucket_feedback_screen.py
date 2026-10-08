"""FLOAT-only screen for MDS emission and block-diagonal RS feedback.

Each feedback block has local_rank checks. With --profiles, retain
occupancy profiles through paired composition and fractional powers.
Otherwise average them first, giving a weaker Jensen bound.
Neither screen is a complete distance certificate.
"""
from math import comb,log2,factorial,prod
from functools import lru_cache
import numpy as np
from hd_systematic_mds_screen import packet_logs as mds_logs,optimize


def packet_logs(z,state=16,width=64,local_rank=1):
    if state%local_rank or width%(state//local_rank):raise ValueError('Equal bucket sizes required')
    tiles=state//local_rank;q=65536.;b=width//tiles;dim=q**state-1
    quiet=np.array([comb(b,l)*sum(comb(l,c)*z**(l-c)*(1-z)**c*q**(-c-min(local_rank,l-c)) for c in range(l+1)) for l in range(b+1)])
    rank=np.array([comb(b,l)*q**(-min(local_rank,l)) for l in range(b+1)])
    k=np.array([1.]);r=np.array([1.])
    for _ in range(tiles):k=np.convolve(k,quiet);r=np.convolve(r,rank)
    values=np.exp2(mds_logs(z,state,width=width))
    for j in range(width+1):
        kj=k[j]/comb(width,j);rj=r[j]/comb(width,j)
        t=(z+(1-z)/q)**j
        v=t*sum(comb(width-j,d)*z**(width-j-d)*(1-z)**d*(q**(state-d)-1)/dim for d in range(min(state,width+1-j)))
        gamma=min(kj,rj*z**max(0,width+1-j-state),dim*v)
        values[j]=((kj,max(0,t-kj)),(v/dim if j else 0,(v if j else v*(q-1)/q)+gamma/q))
    with np.errstate(divide='ignore'):return np.log2(values)


@lru_cache(None)
def occupancies(state,width,local_rank=1):
    if state%local_rank or width%(state//local_rank):raise ValueError('Equal bucket sizes required')
    tiles=state//local_rank;bucket=width//tiles
    def compositions(total,slots):
        if slots==1:yield (total,);return
        for i in range(total+1):
            for rest in compositions(total-i,slots-1):yield (i,)+rest
    # Buckets with at most local_rank inputs have K=q^-l and rank l.
    # Merge profiles with the same degree and all higher occupancies.
    grouped={}
    for row in compositions(tiles,bucket+1):
        j=sum(i*c for i,c in enumerate(row));key=(j,row[local_rank+1:])
        count=factorial(tiles)//prod(factorial(c) for c in row)*prod(comb(bucket,i)**c for i,c in enumerate(row))
        if key in grouped:grouped[key][1]+=count
        else:grouped[key]=[row,count]
    profiles=np.array([x[0] for x in grouped.values()],dtype=int)
    degrees=profiles@np.arange(bucket+1)
    weight=np.array([log2(x[1]) for x in grouped.values()])
    for j in range(width+1):
        assert sum(x[1] for key,x in grouped.items() if key[0]==j)==comb(width,j)
    return profiles,degrees,weight


@lru_cache(1)
def profile_matrices(z,state,width,cutoff,local_rank=1):
    q=65536.;dim=q**state-1;b=width//(state//local_rank)
    profiles,j,weights=occupancies(state,width,local_rank)
    local=np.array([sum(comb(l,c)*z**(l-c)*(1-z)**c*q**(-c-min(local_rank,l-c)) for c in range(l+1)) for l in range(b+1)])
    k=np.exp2(profiles@np.log2(local));t=(z+(1-z)/q)**j
    per_count=np.array([(z+(1-z)/q)**i*sum(comb(width-i,d)*z**(width-i-d)*(1-z)**d*(q**(state-d)-1)/dim for d in range(min(state,width+1-i))) for i in range(width+1)])
    v=per_count[j];ranks=profiles@np.minimum(local_rank,np.arange(b+1))
    gamma=np.minimum.reduce((k,q**(-ranks.astype(float))*z**np.maximum(0,width+1-j-state),dim*v))
    matrices=np.stack((k,np.maximum(0,t-k),np.where(j>0,v/dim,0),np.where(j>0,v,(q-1)/q*v)+gamma/q),axis=1).reshape((-1,2,2))
    keep=j<=cutoff;tail={}
    # Preserve the frequency of deficient ranks in dense packet classes.
    # A single worst profile for an entire degree is far too pessimistic.
    for i in np.flatnonzero(~keep):
        key=(int(j[i]),int(ranks[i]))
        if key in tail:
            tail[key][0]=np.logaddexp2(tail[key][0],weights[i])
            tail[key][1]=np.maximum(tail[key][1],matrices[i])
        else:tail[key]=[weights[i],matrices[i].copy()]
    grouped=np.concatenate((matrices[keep],np.stack([x[1] for x in tail.values()]))) if tail else matrices[keep]
    degrees=np.concatenate((j[keep],np.array([key[0] for key in tail],dtype=int)))
    counts=np.concatenate((weights[keep],np.array([x[0] for x in tail.values()])))
    # Merge numerically identical upward-rounded majorants before forming
    # pairs. This is only a FLOAT proposal; the certificate independently
    # groups exact rationals using directed dyadic rounding.
    mantissa,exponent=np.frexp(grouped)
    rounded=np.ldexp(np.ceil(mantissa*8192)/8192,exponent)
    keys=np.column_stack((degrees,rounded.reshape((-1,4))))
    unique,inverse=np.unique(keys,axis=0,return_inverse=True)
    merged=np.zeros(len(unique))
    np.add.at(merged,inverse,np.exp2(counts))
    degrees=unique[:,0].astype(int);grouped=unique[:,1:].reshape((-1,2,2));counts=np.log2(merged)
    with np.errstate(divide='ignore'):return degrees,counts,np.log2(grouped)


def paired_coefficients(z,alpha,state,width,cutoff=20,local_rank=1,right_rank=None):
    degrees,weight,logs=profile_matrices(z,state,width,cutoff,local_rank)
    rdegrees,rweight,rlogs=(degrees,weight,logs) if right_rank is None else profile_matrices(z,state,width,cutoff,right_rank)
    result=np.zeros((2*width+1,2,2))
    for first in range(0,len(degrees),32):
        left=logs[first:first+32]
        degree=(degrees[first:first+32,None]+rdegrees[None,:]).reshape(-1)
        counts=weight[first:first+32,None]+rweight[None,:]
        for r in range(2):
            for s in range(2):
                entry=np.logaddexp2(left[:,r,0,None]+rlogs[None,:,0,s],left[:,r,1,None]+rlogs[None,:,1,s])
                result[:,r,s]+=np.bincount(degree,weights=np.exp2(counts+alpha*entry).reshape(-1),minlength=2*width+1)
    with np.errstate(divide='ignore'):return np.log2(result)


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows',type=int,nargs='+',default=(1,16,128,256,512,1024))
    p.add_argument('--state',type=int,default=16);p.add_argument('--width',type=int,default=64)
    p.add_argument('--iterations',type=int,default=90)
    p.add_argument('--profiles',action='store_true')
    p.add_argument('--global-marker',action='store_true')
    p.add_argument('--local-rank',type=int,default=1)
    p.add_argument('--cutoff',type=int,default=20)
    p.add_argument('--right-rank',type=int)
    args=p.parse_args()
    for r in args.rows:
        source=lambda z,s,w:packet_logs(z,s,w,args.local_rank)
        pair_source=lambda z,a,s,w:paired_coefficients(z,a,s,w,cutoff=args.cutoff,local_rank=args.local_rank,right_rank=args.right_rank)
        result=optimize(r,args.state,args.iterations,not args.global_marker,width=args.width,packet_source=source,
                        pair_source=pair_source if args.profiles else None)
        result['retained_collision_profiles']=args.profiles
        result['local_rank']=args.local_rank
        result['right_rank']=args.right_rank
        result['construction']='MDS-emission-bucket-feedback';print(json.dumps(result),flush=True)

"""FLOAT screen for block-MDS emission with global MDS feedback.

This changes the construction, not just its circuit. Uniform-state moments
retain every block occupancy; the arbitrary-state quiet cap uses the true
smaller emission distance. No complete certificate is claimed.
"""
from functools import lru_cache
from fractions import Fraction as F
from math import comb,log2,prod
import numpy as np
from hd_bucket_feedback import profiles as exact_profiles
from hd_systematic_mds_screen import optimize


def uniform_mark(empty,dimension,z,q=65536.):
    return sum(comb(empty,c)*z**(empty-c)*(1-z)**c*q**(-min(dimension,c))
               for c in range(empty+1))


def exact_moments(occupancy,z,state,width,local_dimension,q,local_feedback=False):
    z=F(z);j=sum(occupancy);length=width//len(occupancy)
    def flat(n):
        return sum(F(comb(n,c),q**min(local_dimension,c))*z**(n-c)*(1-z)**c for c in range(n+1))
    def quiet(n,h):
        return sum(F(comb(n,c),q**(c+min(h,n-c)))*z**(n-c)*(1-z)**c for c in range(n+1))
    t=(z+(1-z)/q)**j
    k=prod(quiet(n,local_dimension) for n in occupancy) if local_feedback else quiet(j,state)
    v=t*(q**state*prod(flat(length-n) for n in occupancy)-1)/(q**state-1)
    rank=sum(min(local_dimension,n) for n in occupancy) if local_feedback else min(state,j)
    gamma=min(k,z**max(0,length-local_dimension+1-max(occupancy))/q**rank,(q**state-1)*v)
    return k,t,v,gamma


def paired_feedback_caps(occ,z,local_dimension,length,q=65536.):
    """Worst labelling into two larger feedback blocks, for four emissions."""
    dimension=2*local_dimension
    larger=np.array([sum(comb(j,c)*z**(j-c)*(1-z)**c*q**(-c-min(dimension,j-c))
                          for c in range(j+1)) for j in range(2*length+1)])
    quiet_candidates=[];rank_candidates=[]
    for a,b,c,d in ((0,1,2,3),(0,2,1,3),(0,3,1,2)):
        l=occ[:,a]+occ[:,b];r=occ[:,c]+occ[:,d]
        quiet_candidates.append(larger[l]*larger[r])
        rank_candidates.append(np.minimum(dimension,l)+np.minimum(dimension,r))
    return (np.max(quiet_candidates,axis=0),np.min(rank_candidates,axis=0),
            np.min(quiet_candidates,axis=0))


@lru_cache(1)
def matrices(z,state=16,width=64,local_dimension=4,local_feedback=False,lazy_rank=1,feedback_dimension=None):
    if lazy_rank not in (1,2):raise ValueError('Unsupported refresh rank')
    blocks=state//local_dimension;length=width//blocks;q=65536.;dim=q**state-1
    # Do not use the feedback-only low-occupancy grouping: those profiles
    # have equal quiet moments but can have different emission moments.
    items=exact_profiles(state,width,local_dimension)
    profiles=np.array([np.bincount(row,minlength=length+1) for row,count in items])
    degrees=profiles@np.arange(length+1)
    weight=np.array([log2(count) for row,count in items])
    # Express each block moment as q^-h plus a positive polynomial. Expand
    # product(block moments)-q^-state without catastrophic subtraction.
    base=q**(-local_dimension)
    excess=np.array([sum(comb(n,c)*z**(n-c)*(1-z)**c*(q**(-c)-base)
                         for c in range(min(local_dimension,n+1))) for n in range(length+1)])
    running=np.ones(len(items));nonzero=np.zeros(len(items))
    for b in range(blocks):
        e=excess[np.array([length-row[b] for row,count in items])]
        nonzero=nonzero*base+running*e;running*=base+e
    nonzero/=1-q**(-state)
    t=(z+(1-z)/q)**degrees;v=t*nonzero
    quiet=np.array([sum(comb(j,c)*z**(j-c)*(1-z)**c*q**(-c-min(state,j-c)) for c in range(j+1)) for j in range(width+1)])
    if local_feedback:
        local_quiet=np.array([sum(comb(j,c)*z**(j-c)*(1-z)**c*q**(-c-min(local_dimension,j-c)) for c in range(j+1)) for j in range(length+1)])
        k=np.exp2(profiles@np.log2(local_quiet));ranks=profiles@np.minimum(local_dimension,np.arange(length+1))
    else:k=quiet[degrees];ranks=np.minimum(state,degrees)
    quiet_lower=k
    if feedback_dimension is not None:
        if not local_feedback or blocks!=4 or feedback_dimension!=2*local_dimension:
            raise ValueError('Only paired feedback blocks are supported')
        # A sorted emission profile has three possible pairings into the
        # larger feedback blocks. Bound each matrix entry separately:
        # max(K) bounds the quiet transition, but T-min(K), not T-max(K),
        # bounds the complementary transition. Rank is minimized separately.
        occ=np.array([row for row,count in items])
        k,ranks,quiet_lower=paired_feedback_caps(occ,z,local_dimension,length,q)
    largest=np.array([max(row) for row,count in items])
    gamma=np.minimum.reduce((k,q**(-ranks.astype(float))*z**np.maximum(0,length-local_dimension+1-largest),dim*v))
    epsilon=q**(-lazy_rank)
    values=np.stack((k,np.maximum(0,t-quiet_lower),np.where(degrees>0,v/dim,0),np.where(degrees>0,v,(1-epsilon)*v)+gamma*epsilon),axis=1)
    mantissa,exponent=np.frexp(values)
    rounded=np.ldexp(np.ceil(mantissa*8192)/8192,exponent)
    keys=np.column_stack((degrees,rounded))
    unique,inverse=np.unique(keys,axis=0,return_inverse=True)
    weights=np.zeros(len(unique));np.add.at(weights,inverse,np.exp2(weight))
    with np.errstate(divide='ignore'):
        return unique[:,0].astype(int),np.log2(weights),np.log2(unique[:,1:].reshape((-1,2,2)))


def coefficients(z,alpha,state,width,local_dimension=4,local_feedback=False,lazy_rank=1,alternate_feedback=False):
    degrees,weight,logs=matrices(z,state,width,local_dimension,local_feedback,lazy_rank)
    rd,rw,right=matrices(z,state,width,local_dimension,local_feedback,lazy_rank,2*local_dimension) if alternate_feedback else (degrees,weight,logs)
    if alternate_feedback=='eight':degrees,weight,logs=rd,rw,right
    result=np.zeros((2*width+1,2,2))
    for first in range(0,len(degrees),32):
        left=logs[first:first+32]
        degree=(degrees[first:first+32,None]+rd[None,:]).reshape(-1)
        counts=weight[first:first+32,None]+rw[None,:]
        for r in range(2):
            for s in range(2):
                entry=np.logaddexp2(left[:,r,0,None]+right[None,:,0,s],left[:,r,1,None]+right[None,:,1,s])
                result[:,r,s]+=np.bincount(degree,weights=np.exp2(counts+alpha*entry).reshape(-1),minlength=2*width+1)
    with np.errstate(divide='ignore'):return np.log2(result)


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows',type=int,nargs='+',default=(128,512,2048))
    p.add_argument('--dimensions',type=int,nargs='+',default=(4,8))
    p.add_argument('--iterations',type=int,default=60)
    p.add_argument('--regional',action='store_true')
    p.add_argument('--local-feedback',action='store_true')
    p.add_argument('--lazy-rank',type=int,choices=(1,2),default=1,help='Rank two is a proposed refresh ensemble, not a certificate')
    p.add_argument('--alternate-feedback',action='store_true',help='Alternate feedback dimension h and 2h; pointwise pairing caps')
    p.add_argument('--all-eight-feedback',action='store_true',help='Use feedback dimension 2h in every packet')
    args=p.parse_args()
    for h in args.dimensions:
        for r in args.rows:
            feedback_mode='eight' if args.all_eight_feedback else args.alternate_feedback
            pair_source=lambda z,a,s,w:coefficients(z,a,s,w,h,args.local_feedback,args.lazy_rank,feedback_mode)
            result=optimize(r,16,args.iterations,args.regional,pair_source=pair_source)
            result.update(construction='block-MDS-emission',local_dimension=h,local_feedback=args.local_feedback,lazy_rank=args.lazy_rank,alternate_feedback=feedback_mode)
            print(json.dumps(result),flush=True)

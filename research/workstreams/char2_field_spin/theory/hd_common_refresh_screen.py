"""FLOAT screen for Q'=D2 H D1 (Q+C X), Y=X+A Q.

A is block RS, C=A^T, and C A=0. The diagonal-Cauchy map sends every
fixed nonzero argument to a distribution with the established density cap.
The only return to zero is Q=C X. Its marked mass J is bounded using
the invertibility of I+A C. This is a design screen, not a certificate.
"""
from functools import lru_cache
from math import comb
import numpy as np
from hd_bucket_feedback import profiles
from hd_systematic_mds_screen import optimize


@lru_cache(None)
def principal_defects(h):
    if h not in (1,2,4):return (h,)*(4*h+1)
    from hd_common_refresh_principal import nullities
    values=nullities(h);sizes=np.array([i.bit_count() for i in range(len(values))])
    return tuple(int(values[sizes==c].max()) for c in range(4*h+1))


def local_return_cap(l,z,h,length,q=65536.):
    # Positive expansion over output-zero sets: c active positions and
    # e inactive positions. Nonquiet R cannot have h zero RS evaluations.
    nu=principal_defects(h);rho=min(h,l);result=0.
    for c in range(l+1):
        for e in range(min(h,length-l+1)):
            rank_cap=q**(-max(0,c-min(nu[c],max(0,h-l+c))))
            count_cap=q**(max(0,l-c-h)-l)*(q**min(rho,h-e)-1)
            prob=min(1-q**(-rho),rank_cap,count_cap)
            result+=comb(l,c)*comb(length-l,e)*z**(length-c-e)*(1-z)**(c+e)*prob
    return result


@lru_cache(1)
def matrices(z,state=16,width=64,h=4,cross_blocks=False,exact_small=0,scattered=False,joint_rank=False,field_order=65536):
    if field_order not in (256,65536):raise ValueError('Unsupported screening field')
    if field_order!=65536 and (cross_blocks or exact_small or scattered or joint_rank):
        raise ValueError('Subfield screen requires ordinary common refresh')
    q=float(field_order);d=q**state-1;blocks=state//h;length=width//blocks
    items=profiles(state,width,h)
    occ=np.array([row for row,count in items]);degrees=occ.sum(axis=1)
    quiet=np.array([sum(comb(j,c)*z**(j-c)*(1-z)**c*q**(-c-min(h,j-c))
                        for c in range(j+1)) for j in range(length+1)])
    k=np.prod(quiet[occ],axis=1);t=(z+(1-z)/q)**degrees
    base=q**(-h)
    excess=np.array([sum(comb(n,c)*z**(n-c)*(1-z)**c*(q**(-c)-base)
                         for c in range(min(h,n+1))) for n in range(length+1)])
    running=np.ones(len(items));nonzero=np.zeros(len(items))
    for b in range(blocks):
        e=excess[length-occ[:,b]]
        nonzero=nonzero*base+running*e;running*=base+e
    v=t*nonzero/(1-q**(-state))
    # The image of each F^U block under I+A C is l-dimensional. The quiet
    # words are unchanged; any nonquiet word has at least L-h+1-l outputs.
    # Expand product(K+J)-product(K) positively rather than subtracting.
    local_t=(z+(1-z)/q)**np.arange(length+1)
    local_r=np.minimum(np.maximum(0,local_t-quiet),
        (1-q**(-np.minimum(h,np.arange(length+1)).astype(float)))*z**np.maximum(0,length-h+1-np.arange(length+1)))
    local_r=np.minimum(local_r,[local_return_cap(l,z,h,length,q) for l in range(length+1)])
    if exact_small:
        if h!=4 or length!=16 or cross_blocks or not 1<=exact_small<=4:
            raise ValueError('Exact small-support averages require the unmodified four-check map')
        from hd_common_refresh_enumerator import pinned_counts
        for l in range(1,exact_small+1):
            counts,denominator=pinned_counts(l)
            # C_U is injective for l<=4, so the quiet contribution is only
            # the zero vector. Jensen is applied after paired composition.
            average=sum(counts[w]*z**w for w in range(1,17))/denominator
            local_r[l]=min(local_r[l],average)
    total=np.ones(len(items));ret=np.zeros(len(items))
    for b in range(blocks):
        kb=quiet[occ[:,b]];jb=local_r[occ[:,b]]
        ret=ret*kb+total*jb;total*=kb+jb
    if cross_blocks:
        # Proposed C'=(B tensor I_h) C, with [I;B] MDS on block indices.
        # If w feedback blocks are active, at least b+1-w state blocks
        # are active after B. This changes the construction.
        probabilities=np.zeros((len(items),blocks+1));probabilities[:,0]=1
        for b in range(blocks):
            quiet_probability=q**(-np.minimum(h,occ[:,b]).astype(float))
            for w in range(b+1,0,-1):
                probabilities[:,w]=probabilities[:,w]*quiet_probability+probabilities[:,w-1]*(1-quiet_probability)
            probabilities[:,0]*=quiet_probability
        ret=np.zeros(len(items))
        for w in range(1,blocks+1):
            ret+=probabilities[:,w]*z**np.maximum(0,(length-h+1)*(blocks+1-w)-degrees)
        # A k-row restriction of the block mixer is injective on any k
        # selected input blocks. Choose those with the largest local ranks.
        # Subtract the all-quiet event before marking each output support.
        ranks=np.minimum(h,occ);rho=ranks.sum(axis=1);support_cap=np.zeros(len(items))
        for mask in range(1,1<<blocks):
            zeros=blocks-mask.bit_count()
            rank=ranks[:,-zeros:].sum(axis=1) if zeros else np.zeros(len(items),dtype=int)
            probability=np.maximum(0,q**(-rank.astype(float))-q**(-rho.astype(float)))
            floor=sum(np.maximum(0,length-h+1-occ[:,i]) for i in range(blocks) if mask>>i&1)
            support_cap+=probability*z**floor
        ret=np.minimum(ret,support_cap)
        ret=np.minimum(ret,np.maximum(0,t-k))
    if scattered:
        if cross_blocks or exact_small or h!=4 or state!=16 or width!=64:
            raise ValueError('Scattered feedback is a separate K20 construction')
        from hd_scattered_feedback import return_caps
        ret=np.minimum(return_caps(occ,z),np.maximum(0,t-k))
        if joint_rank:
            from hd_scattered_return_rank import inactive_caps,averaged_inactive_caps
            if joint_rank=='superregular':
                from hd_superregular_scatter import return_caps as superregular_caps
                # Different maps: do NOT reuse mode36's return cap.
                ret=np.minimum(superregular_caps(z),np.maximum(0,t-k))
            elif joint_rank in ('scalar8','scalar16'):
                from hd_scalar_scatter import return_caps as scalar_caps
                ret=np.minimum(scalar_caps(z,256 if joint_rank=='scalar8' else 65536),np.maximum(0,t-k))
            elif joint_rank=='zero-diagnostic':ret=np.zeros_like(ret)
            else:ret=np.minimum(ret,averaged_inactive_caps(z) if joint_rank=='average' else inactive_caps(z))
    elif joint_rank:raise ValueError('Joint-rank cap requires scattered feedback')
    ret=np.minimum(ret,d*v)
    values=np.stack((k,np.maximum(0,t-k),ret/d,v),axis=1)
    mantissa,exponent=np.frexp(values)
    rounded=np.ldexp(np.ceil(mantissa*8192)/8192,exponent)
    keys=np.column_stack((degrees,rounded));unique,inverse=np.unique(keys,axis=0,return_inverse=True)
    weights=np.zeros(len(unique));np.add.at(weights,inverse,np.array([count for row,count in items],dtype=float))
    with np.errstate(divide='ignore'):
        return unique[:,0].astype(int),np.log2(weights),np.log2(unique[:,1:].reshape((-1,2,2)))


def coefficients(z,alpha,state,width,h=4,cross_blocks=False,exact_small=0,scattered=False,joint_rank=False):
    degrees,weight,logs=matrices(z,state,width,h,cross_blocks,exact_small,scattered,joint_rank)
    result=np.zeros((2*width+1,2,2))
    for first in range(0,len(degrees),32):
        left=logs[first:first+32]
        degree=(degrees[first:first+32,None]+degrees[None,:]).reshape(-1)
        counts=weight[first:first+32,None]+weight[None,:]
        for r in range(2):
            for s in range(2):
                entry=np.logaddexp2(left[:,r,0,None]+logs[None,:,0,s],left[:,r,1,None]+logs[None,:,1,s])
                result[:,r,s]+=np.bincount(degree,weights=np.exp2(counts+alpha*entry).reshape(-1),minlength=2*width+1)
    with np.errstate(divide='ignore'):return np.log2(result)


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows',type=int,nargs='+',default=(128,256,512,2048))
    p.add_argument('--dimension',type=int,default=4)
    p.add_argument('--iterations',type=int,default=60)
    p.add_argument('--regional',action='store_true')
    p.add_argument('--cross-blocks',action='store_true',help='Proposed fixed MDS mixing of feedback blocks')
    p.add_argument('--exact-small',type=int,choices=(0,1,2,3,4),default=0,help='Exact local return averages up to this potential count')
    p.add_argument('--scattered',action='store_true',help='Shifted check nodes and a coordinate transpose, native mode36')
    p.add_argument('--joint-rank',action='store_true',help='Retain ranks of inactive output-zero constraints')
    p.add_argument('--average-rank',action='store_true',help='Average labelled natural-degree pivot bounds')
    p.add_argument('--superregular',action='store_true',help='PROPOSED Lagrange-basis maps, not timed mode36')
    p.add_argument('--zero-return-diagnostic',action='store_true',help='INVALID operator with return removed: obstruction diagnostic only')
    p.add_argument('--scalar-bits',type=int,choices=(8,16),help='PROPOSED independent nonzero scalar in each emission')
    args=p.parse_args()
    if args.zero_return_diagnostic and args.regional:
        p.error('The zero-return diagnostic requires the global marker (omit --regional)')
    for r in args.rows:
        rank_mode='zero-diagnostic' if args.zero_return_diagnostic else 'scalar'+str(args.scalar_bits) if args.scalar_bits else 'superregular' if args.superregular else 'average' if args.average_rank else args.joint_rank
        source=lambda z,a,s,w:coefficients(z,a,s,w,args.dimension,args.cross_blocks,args.exact_small,args.scattered,rank_mode)
        result=optimize(r,16,args.iterations,args.regional,pair_source=source)
        result.update(construction='common-diagonal-Cauchy-refresh',local_dimension=args.dimension,cross_blocks=args.cross_blocks,exact_small=args.exact_small,scattered=args.scattered,joint_rank=rank_mode,state_mixer=False)
        if args.zero_return_diagnostic:
            result['diagnostic_value']=result.pop('log2_upper_screen')
            result['valid_upper_bound']=False
        print(json.dumps(result),flush=True)

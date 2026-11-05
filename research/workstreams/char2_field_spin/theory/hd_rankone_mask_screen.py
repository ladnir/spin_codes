"""FLOAT labelled state-mask gate for rank-one refresh after feedback.

States0..15 are arbitrary masses on the exact live-block mask, with0 the
zero state; state16 bounds uniform-live density. Active input coordinates
are uniform over the entire field in the Toeplitz comparison, not forced
nonzero. No native code or certified default is changed.
"""
from functools import lru_cache
from math import comb, log2

import numpy as np
from scipy.special import logsumexp

from hd_rankone_common_screen import local_values
from hd_grouped_packet_screen import LIMIT, log_multinomial
from hd_shifted_outer_density import TOEPLITZ_OUTER_MEASURE, log_density_float


def mask_matrices(occupancies, z, *, q=256, h=4):
    """Label-preserving matrices, vectorized over a bounded profile batch.

    For every fixed old state in a mask, different X blocks are independent.
    The old state blocks themselves need not be independent. The flat row
    uses a positive subtraction of the excluded all-zero old state.
    """
    occ=np.asarray(occupancies)
    if (occ.ndim!=2 or not 1<=occ.shape[1]<=4 or occ.shape[0]>4096
            or not np.issubdtype(occ.dtype,np.integer) or np.any(occ<0)
            or np.any(occ>4*h) or h not in (1,2,4) or q<2 or not 0<z<=1):
        raise ValueError('Invalid or over-budget labelled occupancy batch')
    blocks=occ.shape[1]
    masks=1<<blocks
    bits=np.array([[(mask>>i)&1 for i in range(blocks)] for mask in range(masks)])
    length=4*h
    degree=np.arange(length+1)
    quiet,excess,returns=map(np.asarray,local_values(z,q,h))
    t=z+(1-z)/q
    T=t**degree
    born=np.maximum(0.,T-quiet)
    rho=np.minimum(h,degree)
    delta=np.maximum(0,length-h+1-degree)
    A=T*z**delta
    R=float(q)**-rho*np.minimum(t**(degree-rho),z**np.maximum(1,delta))
    R[0]=0.
    local=np.stack((quiet,born,R,A),axis=1).reshape((-1,2,2))
    tensor=np.ones((len(occ),masks,masks))
    total=np.ones((len(occ),masks))
    for i in range(blocks):
        tensor*=local[occ[:,i,None,None],bits[None,:,i,None],bits[None,None,:,i]]
        total*=np.where(bits[None,:,i],A[occ[:,i,None]],T[occ[:,i,None]])
    eta,gamma=1/q,1-1/q
    out=np.zeros((len(occ),masks+1,masks+1))
    out[:,:masks,:masks]=eta*tensor
    out[:,:masks,0]=tensor[:,:,0]
    nonzero=np.minimum(total,tensor[:,:,1:].sum(axis=2))
    nonzero[:,0]=np.maximum(0.,total[:,0]-tensor[:,0,0])
    out[:,:masks,masks]=gamma*nonzero

    # Full-uniform old-state block: each target mask factor is b*Z+D.
    # Product(b*Z+D)-product(b*Z) removes Q_allzero without cancellation.
    base=float(q)**-h
    all_product=np.ones((len(occ),masks))
    removed_product=np.zeros((len(occ),masks))
    live_product=np.ones(len(occ))
    live_nonzero=np.zeros(len(occ))
    for i in range(blocks):
        l=occ[:,i,None]
        bit=bits[None,:,i]
        Z=np.where(bit,born[l],quiet[l])
        D=np.where(bit,T[l]*excess[length-l],base*returns[l])
        removed_product=removed_product*(base*Z)+all_product*D
        all_product*=base*Z+D
        e=excess[length-occ[:,i]]
        live_nonzero=live_nonzero*base+live_product*e
        live_product*=base+e
    normalization=1-float(q)**(-h*blocks)
    V=t**occ.sum(axis=1)*live_nonzero/normalization
    H=np.minimum(removed_product/normalization,V[:,None])
    out[:,masks,:masks]=eta*H
    out[:,masks,0]=H[:,0]
    out[:,masks,masks]=gamma*np.minimum(V,H[:,1:].sum(axis=1))
    if np.any(out<0) or not np.all(np.isfinite(out)):
        raise ArithmeticError('Invalid positive mask cap')
    return out


def packet_matrix(occupancies,z,*,q=256,h=4):
    return mask_matrices([tuple(occupancies)],z,q=q,h=h)[0]


def project_flat_to_masks(matrices,*,q=256,h=4):
    """Merge flat density into its exact nonzero-mask masses before alpha.

    The last state is left as an unreachable zero row/column, allowing
    identical polynomial machinery. In particular mask0 receives nothing.
    """
    values=np.asarray(matrices)
    if values.ndim<2 or values.shape[-1]!=values.shape[-2]:
        raise ValueError('Square mask matrices required')
    masks=values.shape[-1]-1
    if masks<2 or masks&(masks-1) or q<2 or h<1:
        raise ValueError('Invalid mask count or field geometry')
    blocks=masks.bit_length()-1
    denominator=float(q)**(h*blocks)-1
    mass=np.array([0.]+[(float(q)**h-1)**mask.bit_count()/denominator
                       for mask in range(1,masks)])
    out=np.zeros_like(values)
    out[...,:masks,:masks]=values[...,:masks,:masks]+values[...,:masks,masks,None]*mass
    return out


def log_product(left,right):
    return logsumexp((left[:,:,None]+right[None,:,:])*np.log(2),axis=1)/np.log(2)


def log_power(matrix,exponent):
    result=np.full(matrix.shape,-np.inf)
    np.fill_diagonal(result,0.)
    while exponent:
        if exponent&1:result=log_product(result,matrix)
        exponent>>=1
        if exponent:matrix=log_product(matrix,matrix)
    return result


@lru_cache(3)
def half_geometry(maximum):
    return tuple((a+b,f,comb(16,a)*comb(16-a,b)*comb(16-a-b,f),(a+f,b+f))
        for a in range(17) for b in range(17-a) for f in range(17-a-b)
        if a+b+2*f<=maximum)


def pair_coefficients(z,alpha,maximum=64,*,projection=False):
    """Actual fair-swap pair average, retaining the four block labels.

    A full packet uses ordered left/right half profiles. No sorting of
    occupancies and no fresh source-mask geometry is introduced.
    """
    if type(maximum) is not int or not 0<=maximum<=64 or not 0<alpha<=1:
        raise ValueError('Invalid degree or exponent')
    halves=half_geometry(maximum)
    s=np.array([row[0] for row in halves])
    f=np.array([row[1] for row in halves])
    weights=np.array([row[2]/2**row[0] for row in halves])
    occ=np.array([row[3] for row in halves])
    result=np.zeros((33*33,17,17))
    # Precompute only the distinct ordered occupancy tuples. Full cutoff
    # uses about184MiB, instead of retaining every labelled geometry.
    bound=min(16,maximum)+1
    shape=(bound,)*4
    table=np.zeros((bound**4,17,17))
    for start in range(0,len(table),1024):
        indices=np.arange(start,min(start+1024,len(table)))
        profiles=np.array(np.unravel_index(indices,shape)).T
        keep=profiles.sum(axis=1)<=maximum
        if np.any(keep):
            matrices=mask_matrices(profiles[keep],z)
            if projection:matrices=project_flat_to_masks(matrices)
            table[indices[keep]]=matrices**alpha
    for singleton,full,count,left in halves:
        keep=singleton+s+2*(full+f)<=maximum
        right=occ[keep]
        indices=((left[0]*bound+left[1])*bound+right[:,0])*bound+right[:,1]
        targets=(singleton+s[keep])*33+full+f[keep]
        values=table[indices]*((count/2**singleton)*weights[keep])[:,None,None]
        np.add.at(result,targets,values)
    with np.errstate(divide='ignore'):
        return {(s,f):np.log2(result[s*33+f]) for s in range(33)
                for f in range(33-s) if s+2*f<=maximum}


def polynomial_product(left,right,maximum):
    result={}
    for (s,f),a in left.items():
        for (t,g),b in right.items():
            key=s+t,f+g
            if key[0]+2*key[1]>maximum:continue
            value=log_product(a,b)
            result[key]=np.logaddexp2(result[key],value) if key in result else value
    return result


def exact_sparse_regions(z,alpha,maximum):
    if maximum not in (1,2):raise ValueError('Sparse gate only covers1..2')
    base=pair_coefficients(z,alpha,maximum)
    identity=np.full((17,17),-np.inf)
    np.fill_diagonal(identity,0.)
    result,exponent={(0,0):identity},128
    while exponent:
        if exponent&1:result=polynomial_product(result,base,maximum)
        exponent>>=1
        if exponent:base=polynomial_product(base,base,maximum)
    return {(s+2*f,f):matrix-log2(comb(4096,s)*comb(4096-s,f))
            for (s,f),matrix in result.items()}


def sparse_gate(weight,z,alpha=1.):
    from hd_xor_group_geometry import matching_joint_counts
    regions=exact_sparse_regions(z,alpha,weight)
    prefix=alpha*(weight*log_density_float(TOEPLITZ_OUTER_MEASURE)-log2(255)-LIMIT*log2(z))
    rows=[]
    for (a,b),count in matching_joint_counts(8192,weight).items():
        end=log_power(log_product(regions[weight,a],regions[weight,b]),128)[0]
        value=log2(count)+prefix+float(logsumexp(end*np.log(2))/np.log(2))
        rows.append(dict(a=a,b=b,count=count,log2_bound=value))
    return dict(active_groups=weight,z=z,alpha=alpha,rows=rows,
        log2_union_FLOAT=float(logsumexp(np.array([row['log2_bound'] for row in rows])*np.log(2))/np.log(2)),
        directed=False,whole_code_certificate=False,components=17)


def _symmetric_objective(weight,full,z,alpha,projection):
    from hd_shifted_pairs_count import one_pairing_count
    histogram=(4096-weight+full,weight-2*full,full)
    if min(histogram)<=0:raise ValueError('Interior symmetric histogram required')
    coeff=pair_coefficients(z,alpha,projection=projection)
    keys=np.array(tuple(coeff))
    logs=np.array(tuple(coeff.values()))
    count=one_pairing_count(4096,weight,full)
    prefix=log2(count)+alpha*(weight*log_density_float(TOEPLITZ_OUTER_MEASURE)
        -log2(255)-LIMIT*log2(z))-256*log_multinomial(histogram)
    def objective(odds):
        matrix=logsumexp((logs+np.dot(keys,odds)[:,None,None])*np.log(2),axis=0)/np.log(2)
        end=log_power(matrix,32768)[0]
        return prefix+float(logsumexp(end*np.log(2))/np.log(2))-256*np.dot(histogram[1:],odds)
    return histogram,objective


def fixed_symmetric_gate(weight,full,z,alpha,marker_log2,*,projection=False):
    if len(marker_log2)!=2 or not np.all(np.isfinite(marker_log2)):
        raise ValueError('Two finite fixed coefficient markers required')
    _,objective=_symmetric_objective(weight,full,z,alpha,projection)
    return dict(active_groups=weight,a=full,b=full,z=z,alpha=alpha,marker_log2=list(marker_log2),
        log2_bound_FLOAT=float(objective(marker_log2)),evaluations=1,count_method='single-matching-upper',
        directed=False,whole_code_certificate=False,components=16 if projection else 17,
        flat_projection_before_alpha=projection)


def symmetric_positive_gate(weight,full,z,alpha,*,max_evaluations=100):
    from scipy.optimize import minimize
    if type(max_evaluations) is not int or not 1<=max_evaluations<=100:
        raise ValueError('The bounded gate requires a budget in1..100')
    histogram,objective=_symmetric_objective(weight,full,z,alpha,False)
    fit=minimize(objective,np.log2(np.array(histogram[1:])/histogram[0]),method='Nelder-Mead',
        options={'maxfev':max_evaluations})
    return dict(active_groups=weight,a=full,b=full,z=z,alpha=alpha,marker_log2=fit.x.tolist(),
        log2_bound_FLOAT=float(fit.fun),evaluations=fit.nfev,count_method='single-matching-upper',
        directed=False,whole_code_certificate=False,components=17)

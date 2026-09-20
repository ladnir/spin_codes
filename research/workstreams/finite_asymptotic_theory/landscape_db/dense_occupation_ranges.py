"""Dense occupation interval bounds without materializing L coefficients.

One Bernoulli counting-measure majorant is combined with a positive
coefficient bound on the number of eligible positions. The resulting log
bound is convex in integer Q, so both endpoints cover an entire interval.
These are binary64 diagnostics; a coarse interval can leave substantial slack.
"""
import math

import numpy as np
from scipy.special import gammaln

from composition_occupation import terminal_logs


def log_choose(n,q):
    q=np.asarray(q,dtype=float)
    return gammaln(n+1)-gammaln(q+1)-gammaln(n-q+1)


def epoch_mixture_logs(epoch,t,thetas):
    theta=np.asarray(thetas,dtype=float)
    if epoch.shape!=(t+1,3,3) or theta.ndim!=1 or not np.all((theta>0)&(theta<1)):
        raise ValueError('complete epoch coefficients and interior probabilities required')
    degree=np.arange(t+1)
    # Theta batches share a small exact-integer epoch normalizer table.
    normalizer=np.array([math.log(math.comb(t,j)) for j in range(t+1)])
    probability=(normalizer[None,:]+degree[None,:]*np.log(theta)[:,None]
                 +(t-degree)[None,:]*np.log1p(-theta)[:,None])
    result=np.full((len(theta),3,3),-np.inf)
    for j in range(t+1):
        np.logaddexp(result,probability[:,j,None,None]+epoch[j],out=result)
    return result


def density_cost(counts,block,p):
    if not 0<p<1 or not counts or any(not isinstance(n,int) or n<=0 for n in counts.values()):
        raise ValueError('deterministic integer counts/caps and interior p required')
    return max(math.log(n)-math.log(math.comb(block,w))-w*math.log(p)-(block-w)*math.log1p(-p)
               for w,n in counts.items())


def point_logs(occupations,block,length,log_gamma,eta,log_epoch_moment,cutoff,lam):
    """Cauchy coefficient bound, including original row-position counting."""
    q=np.asarray(occupations,dtype=float)
    if np.any((q<0)|(q>length)):
        raise ValueError('occupation outside the row range')
    return (cutoff*lam+log_epoch_moment+block*length*np.logaddexp(0.,eta)
            +q*(log_gamma-block*eta)-(block-1)*log_choose(length,q))


def interval_log(lo,hi,block,length,log_gamma,eta,log_epoch_moment,cutoff,lam):
    if not 1<=lo<=hi<=length:
        raise ValueError('invalid nonzero occupation interval')
    endpoints=point_logs([lo,hi],block,length,log_gamma,eta,log_epoch_moment,cutoff,lam)
    return math.log(hi-lo+1)+float(max(endpoints))


def region_mixture_upper(epoch,t,length,q,p,eta):
    """For independent Bernoulli(p) bits in Q uniformly eligible positions."""
    if length<t or length%t or not 0<=q<=length:
        raise ValueError('non-native region or occupation')
    eligible=math.exp(-float(np.logaddexp(0.,-eta)))
    theta=p*eligible
    matrix=epoch_mixture_logs(epoch,t,np.array([theta]))
    # Return the entrywise matrix upper, retaining every start/end state.
    from activation_q1 import matrix_product
    current=np.full_like(matrix,-np.inf)
    current[:,0,0]=current[:,1,1]=current[:,2,2]=0.
    remaining=length//t
    while remaining:
        if remaining&1:current=matrix_product(current,matrix)
        remaining>>=1
        if remaining:matrix=matrix_product(matrix,matrix)
    correction=length*np.logaddexp(0.,eta)-q*eta-float(log_choose(length,q))
    return current[0]+correction


def serialized_moments(epoch,t,total_bits,ps,etas):
    if total_bits%t or total_bits<t:
        raise ValueError('non-native serialized length')
    p,eta=np.broadcast_arrays(np.asarray(ps,dtype=float),np.asarray(etas,dtype=float))
    if p.ndim!=1 or not np.all((p>0)&(p<1)) or not np.isfinite(eta).all():
        raise ValueError('invalid dense witnesses')
    eligible=np.exp(-np.logaddexp(0.,-eta))
    matrices=epoch_mixture_logs(epoch,t,p*eligible)
    return terminal_logs(matrices,total_bits//t)

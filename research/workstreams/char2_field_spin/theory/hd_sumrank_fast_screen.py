"""Vectorized FLOAT witness search, never a distance certificate.

Evaluates positive zero-mask moments. Exact rational profile routines remain
the independent reference and must reconstruct every accepted witness.
"""

from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb,log2

import numpy as np
from scipy.optimize import minimize,minimize_scalar

from hd_randomized_outer_state_profile import power,total
from hd_sumrank_transfer import minimum_live_weight


@lru_cache(None)
def profile_data(state,tiles=8):
    from scipy.special import gammaln
    a=np.array(tuple(combinations_with_replacement(range(9),tiles)),dtype=np.int64)
    counts=np.stack([(a==j).sum(axis=1) for j in range(9)],axis=1)
    j=a.sum(axis=1)
    log_mult=(gammaln(tiles+1)-gammaln(counts+1).sum(axis=1))/np.log(2)
    log_mult+=np.log2(np.array([comb(8,k) for k in range(9)]))[a].sum(axis=1)
    log_choose=np.log2(np.array([comb(8*tiles,k) for k in range(8*tiles+1)],dtype=float))
    ranks=np.minimum(state,np.minimum(4,a).sum(axis=1))
    floor=np.array([minimum_live_weight(tuple(row),state) for row in a])
    return a,j,log_mult-log_choose[j],ranks,floor,log_choose


def moments(z,state=4,tiles=8,input_mark=1.0):
    a,j,probability,ranks,floor,choose=profile_data(state,tiles)
    n=len(a)
    local=np.array([[comb(8,b)*z**(8-b)*(1-z)**b*65536.0**(-min(i,b))
                     for b in range(9)] for i in range(9)])
    qlocal=np.zeros((9,5)); flocal=np.zeros((9,5))
    for i in range(9):
        for b in range(9):
            qlocal[i,min(4,max(0,i-b))]+=local[i,b]
            flocal[i,min(4,max(0,b-i))]+=local[i,b]
    if input_mark!=1:
        transform=np.array([[comb(i,b)*input_mark**b*(1-input_mark)**(i-b)*65536.0**(-(i-b)) if b<=i else 0
                             for b in range(9)] for i in range(9)])
        qlocal=transform@qlocal
        flocal=transform@flocal
        free_local=transform@local.sum(axis=1)
    else:
        free_local=local.sum(axis=1)
    quiet=np.zeros((n,state+1)); quiet[:,0]=1
    flat=quiet.copy()
    free=np.ones(n)
    for tile in range(tiles):
        following=np.zeros_like(quiet)
        next_flat=np.zeros_like(flat)
        incoming=a[:,tile]
        for cost in range(5):
            qc=qlocal[incoming,cost]; fc=flocal[incoming,cost]
            if cost<state:
                following[:,cost:state]+=quiet[:,:state-cost]*qc[:,None]
                next_flat[:,cost:state]+=flat[:,:state-cost]*fc[:,None]
                following[:,state]+=quiet[:,state-cost:].sum(axis=1)*qc
                next_flat[:,state]+=flat[:,state-cost:].sum(axis=1)*fc
            else:
                following[:,state]+=quiet.sum(axis=1)*qc
                next_flat[:,state]+=flat.sum(axis=1)*fc
        quiet,flat=following,next_flat
        free*=free_local[incoming]
    powers=65536.0**(-np.arange(state+1,dtype=float))
    k=quiet@powers
    v=flat@((powers-powers[-1])/(1-powers[-1]))
    if not np.all(k>0) or not np.all(v>0):
        raise ArithmeticError('FLOAT local moment underflow')
    return k,free,v


def packet_logs(z,state=4,mixer='feedback_refresh',tiles=8,input_mark=1.0):
    a,j,probability,ranks,floor,choose=profile_data(state,tiles)
    k,t,v=moments(z,state,tiles,input_mark)
    old_k,_,_=moments(z,state,tiles) if input_mark!=1 else (k,t,v)
    factor=input_mark**j
    d=65536.0**state-1
    if mixer=='feedback_refresh':
        gamma=np.minimum.reduce((factor*old_k,factor*z**floor*65536.0**(-ranks),d*v))
        live=v+gamma/65536
        live[j==0]=v[j==0]*(65535/65536)+gamma[j==0]/65536
        zero=v/d
    elif mixer=='identity_feedback':
        # Q'=Q+R CX, Y=FX+AQ. Exact transitive feedback is still present,
        # but no separate rank-one state mixer. On CX=0 the weighted
        # point-density cap is Gamma, without the old 1/q contraction.
        # This is only a FLOAT screen for a different construction.
        gamma=np.minimum.reduce((factor*old_k,factor*z**floor*65536.0**(-ranks),d*v))
        live=v+gamma
        live[j==0]=gamma[j==0]
        zero=v/d
    elif mixer=='double_refresh':
        # H=M Q; Y=F X+A H; Q'=H+R C X, with independent uniform
        # transitive M,R. Both zero/live transitions are upper bounded.
        live=v
        zero=v/d
    elif mixer=='uniform_pre':
        live=v
        zero=np.minimum.reduce((v,factor*old_k*(65536.0**ranks-1)/d,
                                factor*np.maximum(0,(1/65536+65535/65536*z)**j-65536.0**(-j))/d,
                                factor*z**floor*65536.0**(-ranks)*(65536.0**ranks-1)/d))
    else:
        raise ValueError('Unsupported screening mixer')
    zero[j==0]=0
    result=np.stack((k,np.maximum(0,t-k),zero,live),axis=1).reshape((-1,2,2))
    with np.errstate(divide='ignore'):
        return np.log2(result)


def fractional(logs,alpha,state,tiles=8):
    a,j,probability,ranks,floor,choose=profile_data(state,tiles)
    result=np.zeros((8*tiles+1,2,2))
    for r in range(2):
        for s in range(2):
            result[:,r,s]=np.bincount(j,weights=np.exp2(probability+alpha*logs[:,r,s]),minlength=8*tiles+1)
    with np.errstate(divide='ignore'):
        return np.log2(result)


def bound(rows,z,alpha,odds,packets,tiles=8,input_mark=1.0,forbidden_weight=838860):
    from scipy.special import logsumexp
    width=tiles*8
    choose=np.log2(np.array([comb(width,j) for j in range(width+1)],dtype=float))
    local=logsumexp((packets+(choose+np.arange(width+1)*odds)[:,None,None])*np.log(2),axis=0)/np.log(2)
    lc=log2(comb(8192,rows))
    return (lc+alpha*(rows*(4096-128*log2(65535))-log2(65535)-forbidden_weight*log2(z)-129*rows*log2(input_mark))
            +total(power(tuple(tuple(row) for row in local),2**21//width)[0])-256*(lc+rows*odds))


def optimize(rows,state=4,mixer='feedback_refresh',tiles=8,iterations=100,input_mark=1.0):
    calls=0
    cached={}
    best=None
    def objective(theta):
        nonlocal calls,best
        ell,alpha=theta
        if not .0001<=ell<=16 or not .005<=alpha<=1:
            return 1e12+sum(abs(x) for x in theta)*1e9
        z=2.0**(-ell)
        if z not in cached:
            cached.clear()
            cached[z]=packet_logs(z,state,mixer,tiles,input_mark)
        packets=fractional(cached[z],alpha,state,tiles)
        chosen=minimize_scalar(lambda p:bound(rows,z,alpha,p,packets,tiles,input_mark),bounds=(-24,20),method='bounded')
        calls+=1
        if best is None or chosen.fun<best['log2_upper_screen']:
            best=dict(active_rows=rows,state=state,mixer=mixer,tiles=tiles,mark=z,alpha=alpha,
                      marker_log_odds=float(chosen.x),log2_upper_screen=float(chosen.fun),calls=calls,
                      input_mark=input_mark,
                      directed=False,whole_code_certificate=False)
        return chosen.fun
    minimize(objective,(max(.002,rows/512),.0625),method='Nelder-Mead',
             options={'maxiter':iterations,'xatol':1e-5,'fatol':.01})
    return best


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=int,nargs='+',default=(16,128,512,1024))
    parser.add_argument('--state',type=int,default=8)
    parser.add_argument('--mixer',choices=('feedback_refresh','identity_feedback','double_refresh','uniform_pre'),default='feedback_refresh')
    parser.add_argument('--tiles',type=int,default=8)
    parser.add_argument('--iterations',type=int,default=60)
    parser.add_argument('--input-marks',type=float,nargs='+',default=(1.0,))
    args=parser.parse_args()
    for r in args.rows:
        for v in args.input_marks:
            print(json.dumps(optimize(r,args.state,args.mixer,args.tiles,args.iterations,v)),flush=True)

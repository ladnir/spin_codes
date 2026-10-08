"""FLOAT-only screen for systematic sum-rank feedback; not a certificate.

Y=X+A Q, Qnext=M Q+T A^T X. Each local A0 is a self-orthogonal
RS[8,4] evaluation map. The LRS blocks separate its residual subspaces.
Unlike the Cauchy feedthrough, an output-zero mask pays only for its
intersection with the potential input support. All intersections are summed.
"""
from math import comb
import numpy as np
from scipy.optimize import minimize,minimize_scalar
from hd_sumrank_fast_screen import profile_data,fractional,bound


def local_moments(z):
    quiet=np.zeros((9,5));flat=np.zeros((9,5));free=np.zeros(9)
    for a in range(9):
        for c in range(a+1):
            for d in range(9-a):
                b=c+d
                value=comb(a,c)*comb(8-a,d)*z**(8-b)*(1-z)**b*65536.0**(-c)
                quiet[a,min(4,a-c)]+=value
                flat[a,min(4,d)]+=value
                free[a]+=value
    return quiet,flat,free


def moments(z,state=16):
    a,j,probability,ranks,floor,choose=profile_data(state)
    qlocal,flocal,free_local=local_moments(z)
    quiet=np.zeros((len(a),state+1));quiet[:,0]=1
    flat=quiet.copy();free=np.ones(len(a))
    for tile in range(8):
        following=np.zeros_like(quiet);next_flat=np.zeros_like(flat)
        incoming=a[:,tile]
        for cost in range(5):
            qc=qlocal[incoming,cost];fc=flocal[incoming,cost]
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
    return quiet@powers,free,flat@((powers-powers[-1])/(1-powers[-1]))


def packet_logs(z,state=16):
    a,j,probability,ranks,floor,choose=profile_data(state)
    k,t,v=moments(z,state);d=65536.0**state-1
    gamma=np.minimum.reduce((k,z**floor*65536.0**(-ranks),d*v))
    live=v+gamma/65536
    live[j==0]=v[j==0]*(65535/65536)+gamma[j==0]/65536
    zero=v/d;zero[j==0]=0
    values=np.stack((k,np.maximum(0,t-k),zero,live),axis=1).reshape((-1,2,2))
    with np.errstate(divide='ignore'):return np.log2(values)


def optimize(rows,state=16,iterations=60):
    cached={};best=None
    def objective(theta):
        nonlocal best
        ell,alpha=theta
        if not .0001<=ell<=16 or not .005<=alpha<=1:return 1e12
        z=2**(-ell)
        if z not in cached:cached.clear();cached[z]=packet_logs(z,state)
        packets=fractional(cached[z],alpha,state)
        chosen=minimize_scalar(lambda p:bound(rows,z,alpha,p,packets),bounds=(-24,20),method='bounded')
        if best is None or chosen.fun<best['log2_upper_screen']:
            best=dict(active_rows=rows,state=state,mark=z,alpha=alpha,
                      marker_log_odds=float(chosen.x),log2_upper_screen=float(chosen.fun),
                      feedthrough='identity',feedback='ideal-transitive',
                      directed=False,whole_code_certificate=False)
        return chosen.fun
    minimize(objective,(max(.002,rows/512),.0625),method='Nelder-Mead',options={'maxiter':iterations})
    return best


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows',type=int,nargs='+',default=(16,128,512,1024))
    p.add_argument('--state',type=int,default=16)
    p.add_argument('--iterations',type=int,default=60)
    args=p.parse_args()
    for r in args.rows:print(json.dumps(optimize(r,args.state,args.iterations)),flush=True)

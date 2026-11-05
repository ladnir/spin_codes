"""FLOAT screen for cheaper subfield randomizers; not a certificate.

If every encoder coefficient belongs to GF256, extension of scalars preserves
symbol distance over GF65536. Both outer and inner setup must use GF256 coins
for this argument. Grouped outer geometry and the pointwise transfer remain
unchanged, with q=256 instead of 65536.

The diagonal-Cauchy density factor can be charged on transitions to nonzero
state, before fractional powers, rather than after every step globally.
This bounds the actual refresh distribution; it is not an exact enumerator.
"""
from math import comb, log2
from time import perf_counter

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp

from hd_common_refresh_screen import matrices
from hd_randomized_outer_state_profile import power, total

LENGTH=1<<21
OUTER=256
ROWS=LENGTH//OUTER
LIMIT=LENGTH*7//20


def density_bits(q,state):
    # Stable even when q**state overflows or rounds q**state-1 to q**state.
    return state*log2(q/(q-1))+float(np.log1p(-float(q)**(-state)))*log2(np.e)


def coefficients(z,alpha,*,q=256,state=16,width=64,h=2,local_density=True,quiet_inverse=False):
    degrees,weights,logs=matrices(z,state,width,h,field_order=q)
    if local_density:
        logs=logs.copy()
        logs[:,:,1]+=density_bits(q,state)
    if quiet_inverse:
        if not local_density:raise ValueError('Inverse-quiet bound already includes actual refresh')
        from hd_grouped_quiet import quiet_screen
        logs[degrees==0,1,1]=np.minimum(logs[degrees==0,1,1],log2(quiet_screen(z,q=q,state=state,width=width,h=h)))
    result=np.zeros((width+1,2,2))
    for r in range(2):
        for c in range(2):
            result[:,r,c]=np.bincount(degrees,weights=np.exp2(weights+alpha*logs[:,r,c]),minlength=width+1)
    with np.errstate(divide='ignore'):return np.log2(result)


def bound(rows,z,alpha,coeff,*,q=256,state=16,width=64,regional=False,local_density=True):
    lc=log2(comb(ROWS,rows))
    prefix=lc+alpha*(rows*(OUTER*log2(q)-(OUTER//2)*log2(q-1))-log2(q-1)-LIMIT*log2(z))
    def local(odds):
        return logsumexp((coeff+(np.arange(width+1)*odds)[:,None,None])*np.log(2),axis=0)/np.log(2)
    def objective(odds):
        return prefix+total(power(tuple(tuple(row) for row in local(odds)),LENGTH//width)[0])-OUTER*(lc+rows*odds)
    chosen=minimize_scalar(objective,bounds=(-24,20),method='bounded')
    value=float(chosen.fun)
    if regional:
        from hd_sumrank_regions import product
        odds=float(chosen.x);p=1/(1+2**(-odds));degree=min(rows,width)
        shifted=coeff[:degree+1]+(np.arange(degree+1)*log2(p)+(width-np.arange(degree+1))*log2(1-p))[:,None,None]
        peak=np.max(shifted);base=(np.exp2(shifted-peak),peak)
        result=(np.eye(2).reshape((1,2,2)),0.);exponent=ROWS//width
        while exponent:
            if exponent&1:result=product(result,base,rows)
            exponent>>=1
            if exponent:base=product(base,base,rows)
        polynomial,scale=result
        if not np.all(polynomial[rows]>0):raise ArithmeticError('Regional FLOAT underflow')
        mass=lc+rows*log2(p)+(ROWS-rows)*log2(1-p)
        region=np.log2(polynomial[rows])+scale-mass
        value=prefix+total(power(tuple(tuple(row) for row in region),OUTER)[0])
    if not local_density:value+=(LENGTH//width)*density_bits(q,state)
    return value,float(chosen.x)


def proposal(rows,*,q=256,state=16,width=64,h=2,iterations=60,regional=False,local_density=True,quiet_inverse=False):
    best=None;start=perf_counter()
    def objective(theta):
        nonlocal best
        ell,alpha=theta
        if not .0001<=ell<=16 or not .005<=alpha<=1:return 1e12
        z=2**(-ell)
        coeff=coefficients(z,alpha,q=q,state=state,width=width,h=h,local_density=local_density,quiet_inverse=quiet_inverse)
        try:value,odds=bound(rows,z,alpha,coeff,q=q,state=state,width=width,regional=regional,local_density=local_density)
        except ArithmeticError:return 1e12
        if not np.isfinite(value):return 1e12
        if best is None or value<best['log2_upper_screen']:
            best=dict(active_rows=rows,field_order=q,state=state,width=width,local_rank=h,
                mark=z,alpha=alpha,marker_log_odds=odds,log2_upper_screen=value,
                local_density=local_density,quiet_inverse=quiet_inverse,regional=regional,directed=False,whole_code_certificate=False)
        return value
    for alpha in ((.0625,.25,1.) if rows<=16 else (.0625,)):
        minimize(objective,(max(.002,rows/512),alpha),method='Nelder-Mead',options={'maxiter':iterations})
    if best is None:raise ArithmeticError('No usable FLOAT proposal')
    best['elapsed_seconds']=perf_counter()-start
    return best


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows',type=int,nargs='+',default=(1,128,512,2048))
    p.add_argument('--q',type=int,choices=(256,65536),default=256)
    p.add_argument('--state',type=int,default=16)
    p.add_argument('--width',type=int,default=64)
    p.add_argument('--rank',type=int,default=2)
    p.add_argument('--iterations',type=int,default=60)
    p.add_argument('--regional',action='store_true')
    p.add_argument('--global-density',action='store_true')
    p.add_argument('--quiet-inverse',action='store_true')
    args=p.parse_args()
    for rows in args.rows:
        print(json.dumps(proposal(rows,q=args.q,state=args.state,width=args.width,h=args.rank,
            iterations=args.iterations,regional=args.regional,local_density=not args.global_density,
            quiet_inverse=args.quiet_inverse)),flush=True)

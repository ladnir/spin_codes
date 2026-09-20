"""Long-region Q1 model with explicit state cancellation.

Every live interval emits at a deterministic rate mu. Each active input
turns a zero state live, or cancels a live state with probability p.
This is an explanatory model, not a finite encoder certificate.
"""
import math
import numpy as np


def region_logs(cancellation,rate,scaled_tilts):
    if not 0<=cancellation<1 or not 0<rate<=1:
        raise ValueError('invalid cancellation probability or live rate')
    a=np.asarray(scaled_tilts,dtype=float)
    if a.ndim!=1 or not np.all(np.isfinite(a)&(a>0)):
        raise ValueError('positive scaled tilts required')
    v=a*rate;h=np.log(-np.expm1(-v))-np.log(v)
    zero=np.full((len(a),4,4),-np.inf);one=np.full_like(zero,-np.inf)
    # Only classes zero (0) and live (2) are used; pad for the checked kernel.
    zero[:,0,0]=0;zero[:,2,2]=-v
    one[:,0,2]=h
    if cancellation: one[:,2,0]=math.log(cancellation)+h
    one[:,2,2]=math.log1p(-cancellation)-v
    return zero,one


def state_region_logs(state_bits,scaled_tilts):
    if state_bits<2: raise ValueError('state dimension must be at least two')
    mass=(1<<state_bits)-1
    return region_logs(1/mass,(1<<(state_bits-1))/mass,scaled_tilts)


def screen(kernel,block,state_bits,spectra):
    grid=np.arange(-120,141)/20

    def values(logs):
        a=np.exp(logs)
        return np.minimum(0.,kernel.coefficients(*state_region_logs(state_bits,a),1,block)+.1*block*a[:,None])

    initial=values(grid);indices=np.argmin(initial,axis=0)
    best=initial[indices,np.arange(block+1)];witnesses=grid[indices]
    relevant=set()
    for counts in spectra.values():
        terms={w:math.log(n)+best[w] for w,n in counts.items()}
        relevant.update(w for w,v in terms.items() if v>=max(terms.values())-24*math.log(2))
    extra=np.unique(np.concatenate([witnesses[w]+np.arange(-10,11)/200 for w in relevant]))
    candidates=values(extra);indices=np.argmin(candidates,axis=0)
    improved=candidates[indices,np.arange(block+1)]<best
    best[improved]=candidates[indices,np.arange(block+1)][improved]
    witnesses[improved]=extra[indices[improved]]
    from study_engineering_surfaces import summarize_coefficients
    return {family:summarize_coefficients(block,1,counts,best,witnesses) for family,counts in spectra.items()}

"""Bounded FLOAT joint proposal for the 17-state rank-one mask cap.

At fixed positive markers x,y, the half-occupancy weight at (i,j) is
[u^i v^j](1+x(u+v)/2+yuv)^16. The packet weight is the product of the
two ordered half weights. Contracting these weights directly against
M_occupancy entrywise**alpha equals the established polynomial evaluated
at x,y; it avoids building every coefficient. No symmetry is assumed.

Only one medium source profile is explored. This is neither directed
arithmetic nor a whole-code proof. Existing literal negative controls
and construction defaults are untouched.

Recorded bounded control:47 proposal evaluations plus one agreement check
gave +92363.8038867295 bits, versus +99844.04545436724 initially. The best
literal witness is BOUNDED_BEST below. Its Perron mask-size weights were
(23.722,20.859,27.797,19.032,7.446) percent for sizes0..4, with1.143 percent
in the flat component. These are weights of the majorant, not empirical
or true state probabilities. No single small mask or flat tail dominates.
One separate fixed two-refresh gate gave +65853.25259219245 bits. Neither
screen passes, and no two-refresh optimization was performed.

An unimplemented backup alternates rank-one and DHD packets. At DHD
boundaries, keep only zero/flat components. Compose a2-by17 rank-one cap
with a17-by2 DHD cap BEFORE alpha. The DHD cap sends each zero-return cap
to zero and its nonzero total cap to flat density multiplied by
b=(q**16-1)/(q-1)**16. This removes persistent mask tails across pairs,
but still requires an actual paired-placement calculation and a new gate.
"""
from math import log2

import numpy as np
from scipy.special import logsumexp

from hd_grouped_packet_screen import LIMIT, log_multinomial
from hd_rankone_mask_screen import half_geometry, mask_matrices, log_power
from hd_rankone_orbit_screen import FROZEN_MEDIUM
from hd_shifted_outer_density import TOEPLITZ_OUTER_MEASURE, log_density_float
from hd_shifted_pairs_count import one_pairing_count


BOUNDED_BEST=dict(z=.4353679232742893,alpha=.17088472633321333,
                  marker_log2=(-3.3650867183232807,-1.996947766004256))
RECORDED_SCREENS=dict(proposal_evaluations=47,initial_agreement_evaluations=1,
    one_refresh_log2_bound_FLOAT=92363.8038867295,
    two_refresh_log2_bound_FLOAT=65853.25259219245,
    directed=False,whole_code_certificate=False)


def refresh_power(cap, steps, *, q=256):
    """Compose independent rank-one refreshes, after the same feedback.

    The exact kernel is a**steps I+(q**(-steps)-a**steps)L
    +(1-q**(-steps))J. I and L preserve each mask; J is uniform-live.
    Every zero successor stays zero, so its column is left unchanged.
    """
    if type(steps) is not int or steps not in (1,2):
        raise ValueError('Only one or two independent updates are allowed')
    result=np.array(cap,copy=True)
    if steps==2:
        result[...,1:-1]/=q
        result[...,-1]*=(1-q**-2)/(1-q**-1)
    return result


def marked_matrix(z, alpha, marker_log2, *, maximum=64, refresh_steps=1):
    """Log2 marked 17-state matrix, streamed through bounded batches."""
    odds=np.asarray(marker_log2,dtype=float)
    if (not np.isfinite(z) or not 0<z<=1 or not np.isfinite(alpha)
            or not 0<alpha<=1 or odds.shape!=(2,) or not np.all(np.isfinite(odds))
            or type(maximum) is not int or not 0<=maximum<=64
            or type(refresh_steps) is not int or refresh_steps not in (1,2)):
        raise ValueError('Invalid marks, degree, or refresh count')
    log_base=float(logsumexp(np.array((0.,*odds))*np.log(2))/np.log(2))
    bound=min(16,maximum)+1
    half=np.zeros((bound,bound))
    for singleton,full,count,occupancy in half_geometry(maximum):
        half[occupancy]+=np.exp2(log2(count)+singleton*(odds[0]-1)+full*odds[1]-16*log_base)
    result=np.zeros((17,17))
    shape=(bound,)*4
    for start in range(0,bound**4,1024):
        indices=np.arange(start,min(start+1024,bound**4))
        profiles=np.array(np.unravel_index(indices,shape)).T
        profiles=profiles[profiles.sum(axis=1)<=maximum]
        if not len(profiles):continue
        weights=half[profiles[:,0],profiles[:,1]]*half[profiles[:,2],profiles[:,3]]
        caps=refresh_power(mask_matrices(profiles,z),refresh_steps)
        result+=np.einsum('p,pij->ij',weights,caps**alpha,optimize=False)
    with np.errstate(divide='ignore'):
        return np.log2(result)+32*log_base


def perron_diagnostic(log_matrix):
    """Doob/Perron state and edge weights of the MAJORANT, not the code."""
    scale=float(np.max(log_matrix))
    matrix=np.exp2(log_matrix-scale)
    values,right=np.linalg.eig(matrix)
    selected=int(np.argmax(values.real))
    eigenvalue=float(values[selected].real)
    if abs(values[selected].imag)>1e-10 or eigenvalue<=0:
        raise ArithmeticError('No positive Perron eigenvalue')
    v=right[:,selected].real
    left_values,left=np.linalg.eig(matrix.T)
    u=left[:,int(np.argmin(abs(left_values-eigenvalue)))].real
    if v.sum()<0:v=-v
    if u.sum()<0:u=-u
    if min(u.min(),v.min()) < -1e-10:
        raise ArithmeticError('Unexpected signed Perron vectors')
    u=np.maximum(u,0.);v=np.maximum(v,0.)
    inner=float(u@v)
    state=u*v/inner
    edges=u[:,None]*matrix*v[None,:]/(eigenvalue*inner)
    return dict(log2_eigenvalue=scale+log2(eigenvalue),
        state_weights=state.tolist(),flat_weight=float(state[16]),
        live_mask_weight_by_size=[float(sum(state[s] for s in range(16) if s.bit_count()==k))
                                  for k in range(5)],
        flat_self_edge=float(edges[16,16]),zero_self_edge=float(edges[0,0]),
        nonflat_to_flat=float(edges[:16,16].sum()),flat_to_nonflat=float(edges[16,:16].sum()),
        interpretation='Perron weights of the positive majorant, not actual state probabilities')


def fixed_gate(z,alpha,marker_log2,*,refresh_steps=1,diagnostic=False):
    weight,full=1024,400
    histogram=(4096-weight+full,weight-2*full,full)
    matrix=marked_matrix(z,alpha,marker_log2,refresh_steps=refresh_steps)
    end=log_power(matrix,32768)[0]
    prefix=(log2(one_pairing_count(4096,weight,full))
        +alpha*(weight*log_density_float(TOEPLITZ_OUTER_MEASURE)-log2(255)-LIMIT*log2(z))
        -256*log_multinomial(histogram)-256*np.dot(histogram[1:],marker_log2))
    result=dict(active_groups=weight,a=full,b=full,z=z,alpha=alpha,
        marker_log2=list(marker_log2),refresh_steps=refresh_steps,
        log2_bound_FLOAT=prefix+float(logsumexp(end*np.log(2))/np.log(2)),
        directed=False,whole_code_certificate=False,components=17,
        count_method='single-matching-upper')
    if diagnostic:result['perron']=perron_diagnostic(matrix)
    return result


def bounded_proposal(*,max_evaluations=48,progress=None):
    """Eight deterministic seeds, then one bounded four-variable simplex.

    Coordinates are (-log2(z), log2(alpha), log2(x), log2(y)). Invalid
    proposals return infinity without building any packet matrices. The
    maximum counts every uncached finite objective evaluation, including
    the frozen agreement control, and is never enlarged by the optimizer.
    """
    from scipy.optimize import minimize
    if type(max_evaluations) is not int or not 8<=max_evaluations<=48:
        raise ValueError('An explicit objective budget in8..48 is required')
    old=FROZEN_MEDIUM
    cache={}
    def evaluate(point):
        point=tuple(float(v) for v in point)
        if point in cache:return cache[point]['log2_bound_FLOAT']
        if (not np.all(np.isfinite(point)) or not .001<=point[0]<=8
                or not -8<=point[1]<=0 or min(point[2:]) < -24 or max(point[2:])>12):
            return float('inf')
        if len(cache)>=max_evaluations:return float('inf')
        value=fixed_gate(2**-point[0],2**point[1],point[2:])
        cache[point]=value
        if progress is not None:progress(len(cache),value)
        return value['log2_bound_FLOAT']
    base=(-log2(old['z']),log2(old['alpha']),*old['marker_log2'])
    seeds=[base]
    seeds.extend((base[0],log2(a),base[2],base[3]) for a in (.25,.5,1.))
    seeds.extend((-log2(z),log2(a),base[2],base[3])
                 for z,a in ((.7,.25),(.8,.5),(.9,1.),(.7,1.)))
    for seed in seeds:evaluate(seed)
    best=min(cache,key=lambda p:cache[p]['log2_bound_FLOAT'])
    simplex=np.tile(best,(5,1))
    for i,step in enumerate((.15,-.75,.75,.75)):
        simplex[i+1,i]+=step
    minimize(evaluate,np.array(best),method='Nelder-Mead',
        options={'maxfev':max_evaluations-len(cache)+1,'initial_simplex':simplex,
                 'xatol':1e-5,'fatol':1e-4})
    best=min(cache,key=lambda p:cache[p]['log2_bound_FLOAT'])
    return dict(best=cache[best],evaluations=len(cache),max_evaluations=max_evaluations,
                rows=list(cache.values()),directed=False,whole_code_certificate=False)


if __name__=='__main__':
    import argparse
    import json
    parser=argparse.ArgumentParser()
    choice=parser.add_mutually_exclusive_group()
    choice.add_argument('--search',action='store_true')
    choice.add_argument('--best',action='store_true')
    parser.add_argument('--max-evaluations',type=int,default=47)
    parser.add_argument('--two-refresh',action='store_true')
    args=parser.parse_args()
    if args.two_refresh and not args.best:
        parser.error('--two-refresh requires --best; it is not a search mode')
    if args.search:
        def progress(n,row):
            print(json.dumps(dict(evaluation=n,**row),allow_nan=False),flush=True)
        print(json.dumps(bounded_proposal(max_evaluations=args.max_evaluations,progress=progress),
                         allow_nan=False),flush=True)
    elif args.best:
        print(json.dumps(fixed_gate(**BOUNDED_BEST,refresh_steps=2 if args.two_refresh else 1,
                                   diagnostic=True),indent=2,allow_nan=False))
    else:
        old=FROZEN_MEDIUM
        print(json.dumps(fixed_gate(old['z'],old['alpha'],old['marker_log2'],diagnostic=True),
                         indent=2,allow_nan=False))

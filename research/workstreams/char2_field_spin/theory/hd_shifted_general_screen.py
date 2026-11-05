"""FLOAT screening of arbitrary oriented supports for alternating pair routes.

Support profile = (active even groups, active odd groups, full even pairs,
full odd pairs). Unlike the first shifted screen, neither matching must be
entirely full. This uses nonnegative coefficient bounds, not a certificate.
"""
from functools import lru_cache
from math import comb, log, log2
import argparse
import json
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from hd_grouped_packet_screen import LIMIT, OCC, local_logs, log_multinomial
from hd_randomized_outer_state_profile import multiply, power, total
from hd_shifted_general_count import (joint_count_upper, exact_joint_count,
                                     matching_histograms)
from hd_shifted_pairs_screen import three_type_geometry, positive_alias_power
from hd_sumrank_regions import product

PAIRS = 4096


@lru_cache(128)
def support_bits(profile,exact):
    count = (exact_joint_count if exact else joint_count_upper)(PAIRS,*profile)
    if count<=0:
        raise ValueError('Empty support profile')
    return log2(count)


@lru_cache(None)
def oriented_geometry():
    halves = np.array([(a,b) for a in range(17) for b in range(17)])
    entries = np.full((289,17), -np.inf)
    degrees = np.zeros((289,17,3), dtype=int)
    for index,(a,b) in enumerate(halves):
        for full in range(max(0,a+b-16),min(a,b)+1):
            counts = (16-a-b+full,a-full,b-full,full)
            entries[index,full] = log_multinomial(counts)
            degrees[index,full] = counts[1:]
    combined = np.sort(np.concatenate((np.repeat(halves,289,axis=0),
                                      np.tile(halves,(289,1))),axis=1),axis=1)
    radix = np.array([17**3,17**2,17,1])
    lookup = np.full(17**4,-1,dtype=int)
    lookup[OCC@radix] = np.arange(len(OCC))
    indices = lookup[combined@radix]
    if np.any(indices<0):
        raise ArithmeticError('Missing oriented occupancy profile')
    return entries,degrees,indices


def general_step(z,alpha,odds):
    entries,degrees,indices = oriented_geometry()
    halves = logsumexp((entries+degrees@np.asarray(odds))*log(2),axis=1)/log(2)
    weights = (halves[:,None]+halves[None,:]).ravel()
    logs = local_logs(z)[indices]
    return logsumexp((weights[:,None,None]+alpha*logs)*log(2),axis=0)/log(2)


def general_coefficients(z,alpha,odds):
    counts,weights,degrees,indices = three_type_geometry()
    halves = weights+counts@np.asarray(odds)
    weights = (halves[:,None]+halves[None,:]).ravel()
    logs = local_logs(z).reshape((-1,4))
    result = np.empty((4,33**3))
    for entry in range(4):
        result[entry] = np.bincount(degrees,weights=np.exp2(weights+alpha*logs[indices,entry]),minlength=33**3)
    return result.reshape((4,33,33,33))


def homogeneous_coefficients(z,alpha,category,steps=2):
    """Compose state paths before alpha when all source packets share a type."""
    if category not in (0,1,2) or steps not in (1,2):
        raise ValueError('Homogeneous type/composition unsupported')
    half = np.array([(a,b) for a in range(17) for b in range(17)])
    degree = half.sum(axis=1)
    weights = np.array([log2(comb(16,int(a)))+log2(comb(16,int(b))) for a,b in half])
    if category==2:
        occ = np.column_stack((half[:,0],half[:,0],half[:,1],half[:,1]))
    else:
        occ = np.column_stack((half[:,0],np.zeros(289,dtype=int),half[:,1],np.zeros(289,dtype=int)))
    logs = local_logs(z,occupancies=occ)
    if steps==2:
        logs = np.logaddexp2(logs[:,None,:,0,None]+logs[None,:,None,0,:],
                            logs[:,None,:,1,None]+logs[None,:,None,1,:]).reshape((-1,2,2))
        weights = (weights[:,None]+weights[None,:]).ravel()
        degree = (degree[:,None]+degree[None,:]).ravel()
    result = np.empty((32*steps+1,2,2))
    for row in range(2):
        for col in range(2):
            result[:,row,col] = np.bincount(degree,weights=np.exp2(weights+alpha*logs[:,row,col]),minlength=32*steps+1)
    return result


def regional(z,alpha,histogram,odds,method='global',grid=64):
    target = np.array(histogram[1:],dtype=int)
    normalizer = log_multinomial(histogram)+target@np.asarray(odds)
    active = np.flatnonzero(target)
    if method=='global' or (method=='hybrid' and len(active)!=1):
        step = general_step(z,alpha,odds)
        return np.array(power(tuple(map(tuple,step)),128))-normalizer, {}
    if method not in ('alias','hybrid'):
        raise ValueError('Unknown regional coefficient method')
    if len(active)==1:
        axis = int(active[0])
        values = homogeneous_coefficients(z,alpha,axis,steps=2)
        values *= np.exp2(np.arange(len(values))*odds[axis])[:,None,None]
        degree = int(target[axis])
        peak = float(np.max(values))
        base = (values/peak,log2(peak))
        run = (np.eye(2).reshape((1,2,2)),0.)
        exponent = 64
        while exponent:
            if exponent&1:
                run = product(run,base,degree)
            exponent >>= 1
            if exponent:
                base = product(base,base,degree)
        values,scale = run
        if not np.all(values[degree]>0):
            raise ArithmeticError('Univariate regional FLOAT underflow')
        return np.log2(values[degree])+scale-normalizer, {'method':'univariate-polynomial','composition_steps':2}
    coefficients = general_coefficients(z,alpha,odds)
    # Coefficients beyond a target in any coordinate cannot contribute to
    # that coefficient of a positive power. Removing them is exact here.
    for axis,maximum in enumerate(target,1):
        if maximum<32:
            index = [slice(None)]*4
            index[axis] = slice(maximum+1,None)
            coefficients[tuple(index)] = 0
    logs,diagnostic = positive_alias_power(coefficients,target,grid)
    return logs-normalizer, diagnostic


def bound(profile,z,alpha,odds,*,method='global',grid=64,exact_count=False):
    profile = tuple(map(int,profile))
    count_bits = support_bits(profile,exact_count)
    if not 0<z<1 or not 0<alpha<=1 or np.asarray(odds).shape!=(2,3):
        raise ValueError('Invalid mark, alpha, or marker array')
    histograms = matching_histograms(PAIRS,*profile)
    first,d0 = regional(z,alpha,histograms[0],odds[0],method,grid)
    second,d1 = regional(z,alpha,histograms[1],odds[1],method,grid)
    run = power(multiply(tuple(map(tuple,first)),tuple(map(tuple,second))),128)
    r = profile[0]+profile[1]
    prefix = count_bits+alpha*(r*(2048-128*log2(255))-log2(255)-LIMIT*log2(z))
    return dict(profile=list(profile),active_groups=r,mark=z,alpha=alpha,
                marker_log_odds=np.asarray(odds).tolist(),log2_support_count=count_bits,
                exact_support_count=exact_count,coefficient_method=method,grid=grid,
                log2_upper_screen=float(prefix+total(run[0])),regional_diagnostics=[d0,d1],
                directed=False,whole_code_certificate=False)


def proposal(profile,iterations=100,exact_count=False,start=None,method='global'):
    histograms = matching_histograms(PAIRS,*profile)
    odds = np.array([[log2(n/max(1,hist[0])) if n else -40. for n in hist[1:]]
                     for hist in histograms])
    live = np.array([n>0 for hist in histograms for n in hist[1:]])
    tied = None
    if tuple(histograms[0])==tuple(histograms[1]):
        tied = (0,1,2)
    elif tuple(histograms[0])==tuple(histograms[1][i] for i in (0,2,1,3)):
        # The local h4 bound is invariant under exchanging the two lanes
        # of every pair. This is an equality of these operators, not a
        # claim that arbitrary stateful region operators commute.
        tied = (1,0,2)
    if tied is not None:
        live[3:] = False
    if method=='hybrid':
        # A coefficient extracted exactly as a polynomial is independent
        # of its positive tilt. Searching this flat direction only exposes
        # FLOAT noise and can push intermediate coefficients to underflow.
        for parity,hist in enumerate(histograms):
            if sum(n>0 for n in hist[1:])==1:
                live[3*parity:3*parity+3] = False
    initial = [max(.002,sum(profile[:2])/512),.08]+odds.ravel()[live].tolist()
    if start is not None:
        initial[:2] = [-log2(start['mark']),start['alpha']]
        initial[2:] = np.array(start['marker_log_odds']).ravel()[live].tolist()
    best = None
    def objective(theta):
        nonlocal best
        ell,alpha,*markers = theta
        if not .00001<=ell<=16 or not .001<=alpha<=1 or any(abs(x)>32 for x in markers):
            return 1e12
        full = odds.ravel().copy()
        full[live] = markers
        if tied is not None:
            full[3:] = full[np.array(tied)]
        try:
            result = bound(profile,2**-ell,alpha,full.reshape((2,3)),exact_count=exact_count,method=method)
        except ArithmeticError:
            return 1e12
        value = result['log2_upper_screen']
        if np.isfinite(value) and (best is None or value<best['log2_upper_screen']):
            best = result
        return value if np.isfinite(value) else 1e12
    simplex = np.tile(initial,(len(initial)+1,1))
    for axis,value in enumerate(initial):
        simplex[axis+1,axis] += .25 if axis>=2 else max(.002,.05*abs(value))
    minimize(objective,initial,method='Nelder-Mead',
             options={'maxiter':iterations,'initial_simplex':simplex})
    if best is None:
        raise ArithmeticError('No finite general profile screen')
    return best


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile',nargs=4,type=int,required=True,metavar=('EVEN','ODD','EVEN_FULL','ODD_FULL'))
    parser.add_argument('--iterations',type=int,default=100)
    parser.add_argument('--alias',action='store_true')
    parser.add_argument('--grid',type=int,default=64)
    parser.add_argument('--exact-count',action='store_true')
    parser.add_argument('--hybrid',action='store_true',help='Extract homogeneous regional coefficients during search')
    args = parser.parse_args()
    result = proposal(args.profile,args.iterations,args.exact_count,method='hybrid' if args.hybrid else 'global')
    print(json.dumps(result),flush=True)
    if args.alias:
        print(json.dumps(bound(args.profile,result['mark'],result['alpha'],result['marker_log_odds'],
                               method='alias',grid=args.grid,exact_count=args.exact_count)),flush=True)

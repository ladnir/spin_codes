"""FLOAT h4 screen for shifted pairings, retaining singleton orientations.

Classes consist of unions of `active_pairs` even matching edges, with `runs`
active runs in the contracted cycle. The odd matching has runs of each
singleton orientation and active_pairs-runs full pairs. Independent packet
permutations remain mandatory in every region. No directed certificate.
"""
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, log, log2

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from hd_grouped_packet_screen import LIMIT, OCC, coefficients, local_logs, log_multinomial
from hd_randomized_outer_state_profile import multiply, power, total
from hd_shifted_pairs_count import circular_run_count
from hd_sumrank_regions import product


@lru_cache(None)
def half_geometry():
    halves = np.array(tuple(combinations_with_replacement(range(17),2)))
    terms = []
    for a,b in halves:
        terms.append(tuple((a+b-2*c,c,log_multinomial((16-a-b+c,a-c,b-c,c)))
                           for c in range(max(0,a+b-16),min(a,b)+1)))
    combined = np.sort(np.concatenate((np.repeat(halves,len(halves),axis=0),
                                       np.tile(halves,(len(halves),1))),axis=1),axis=1)
    code = lambda value:value@np.array([17**3,17**2,17,1])
    lookup = np.full(17**4,-1,dtype=int)
    lookup[code(OCC)] = np.arange(len(OCC))
    indices = lookup[code(combined)]
    if np.any(indices<0):
        raise ArithmeticError('Missing h4 occupancy profile')
    return halves,tuple(terms),indices


def half_weights(singleton_odds,full_odds):
    halves,terms,_ = half_geometry()
    return np.array([logsumexp(np.array([n+singles*singleton_odds+full*full_odds
                                         for singles,full,n in row])*log(2))/log(2)
                     +(1 if a!=b else 0) for (a,b),row in zip(halves,terms)])


def mixed_step(z,alpha,singleton_odds,full_odds):
    """One-step polynomial evaluated at x10=x01, keeping their true counts."""
    weights = half_weights(singleton_odds,full_odds)
    _,_,indices = half_geometry()
    profiles = local_logs(z)[indices]
    exponents = (weights[:,None]+weights[None,:]).ravel()[:,None,None]+alpha*profiles
    return logsumexp(exponents*log(2),axis=0)/log(2)


@lru_cache(None)
def three_type_geometry():
    """All labeled half profiles, for three-variable coefficients."""
    counts = np.array([(a,b,c) for a in range(17) for b in range(17-a)
                       for c in range(17-a-b)],dtype=np.int32)
    weights = np.array([log_multinomial((16-a-b-c,a,b,c)) for a,b,c in counts])
    half_codes = counts@np.array([33**2,33,1],dtype=np.int32)
    degree_codes = (half_codes[:,None]+half_codes[None,:]).ravel()
    occ = counts[:,:2]+counts[:,2,None]
    combined = np.sort(np.concatenate((np.repeat(occ,len(occ),axis=0),
                                       np.tile(occ,(len(occ),1))),axis=1),axis=1)
    radix = np.array([17**3,17**2,17,1])
    lookup = np.full(17**4,-1,dtype=np.int32)
    lookup[OCC@radix] = np.arange(len(OCC))
    indices = lookup[combined@radix]
    return counts,weights,degree_codes,indices


def mixed_coefficients(z,alpha,singleton_odds,full_odds):
    counts,weights,degrees,indices = three_type_geometry()
    half_weights = weights+(counts[:,0]+counts[:,1])*singleton_odds+counts[:,2]*full_odds
    weight = (half_weights[:,None]+half_weights[None,:]).ravel()
    logs = local_logs(z).reshape((-1,4))
    result = np.empty((4,33**3))
    for entry in range(4):
        result[entry] = np.bincount(degrees,weights=np.exp2(weight+alpha*logs[indices,entry]),minlength=33**3)
    return result.reshape((4,33,33,33))


def positive_alias_power(coefficients,target,grid,*,steps=128,cushion=1e-12):
    """FLOAT root-of-unity alias sum after Perron similarity and scaling.

    The exact alias sum majorizes its target coefficient by positivity.
    The optional cushion is a numerical diagnostic margin, not a proved
    floating error enclosure. A certificate requires directed arithmetic.
    """
    if grid<33 or grid&(grid-1) or steps<1 or steps&(steps-1):
        raise ValueError('Use a power-of-two grid >=64 and power-of-two steps')
    positive = coefficients.sum(axis=(1,2,3)).reshape((2,2))
    a,b,c,d = positive.ravel()
    eigenvalue = (a+d+np.sqrt((a-d)**2+4*b*c))/2
    vector = np.array([b,eigenvalue-a] if eigenvalue-a>=eigenvalue-d else [eigenvalue-d,c])
    if not eigenvalue>0 or not np.all(vector>0):
        raise ArithmeticError('Degenerate Perron scaling')
    vector /= np.max(vector)
    factors = (vector[None,:]/vector[:,None]/eigenvalue).ravel()
    normalized = coefficients*factors[:,None,None,None]
    row_sums = normalized.sum(axis=(1,2,3)).reshape((2,2)).sum(axis=1)
    if np.max(abs(row_sums-1))>1e-10:
        raise ArithmeticError('Perron-scaled row sums do not normalize')
    values = np.fft.fftn(normalized,s=(grid,grid,grid),axes=(1,2,3))
    for _ in range(steps.bit_length()-1):
        a,b,c,d = values
        trace,bc = a+d,b*c
        values = np.stack((a*a+bc,b*trace,c*trace,d*d+bc))
    aliases = np.fft.ifftn(values,axes=(1,2,3))
    index = tuple(int(n)%grid for n in target)
    selected = aliases[(slice(None),)+index]
    if np.max(abs(selected.imag))>1e-10 or np.min(selected.real)<-1e-10:
        raise ArithmeticError('FLOAT alias cancellation diagnostic failed')
    diagnostic_floor = 100*max(1e-12,cushion)
    if np.min(selected.real)<diagnostic_floor:
        raise ArithmeticError('FLOAT alias target coefficient is below the normalized cancellation floor')
    with np.errstate(divide='ignore'):
        logs = np.log2(np.maximum(0,selected.real)+cushion).reshape((2,2))
    logs += steps*log2(eigenvalue)+np.log2(vector[:,None]/vector[None,:])
    return logs,dict(grid=grid,perron_value=float(eigenvalue),perron_vector=vector.tolist(),
        perron_row_sum_error=float(np.max(abs(row_sums-1))),
        minimum_normalized_alias=float(np.min(selected.real)),
        maximum_target_imaginary=float(np.max(abs(selected.imag))),
        normalized_numerical_cushion=cushion,normalized_diagnostic_floor=diagnostic_floor,directed=False)


def alias_screen(active_pairs,runs,z,alpha,odds,*,grid=64,mode='alternating',even_steps=2):
    even,global_odd = regional_matrices(active_pairs,runs,z,alpha,odds,even_regional=True,even_steps=even_steps)
    _,singleton_odds,full_odds = odds
    coefficients = mixed_coefficients(z,alpha,singleton_odds,full_odds)
    target = (runs,runs,active_pairs-runs)
    odd,diagnostics = positive_alias_power(coefficients,target,grid)
    normalizer = log_multinomial((4096-active_pairs-runs,)+target)+2*runs*singleton_odds+(active_pairs-runs)*full_odds
    odd -= normalizer
    if mode=='alternating':
        run = power(multiply(tuple(map(tuple,even)),tuple(map(tuple,odd))),128)
    elif mode=='random':
        run = power(tuple(map(tuple,np.logaddexp2(even,odd)-1)),256)
    else:
        raise ValueError('Unknown pairing schedule')
    prefix = log2(circular_run_count(4096,active_pairs,runs))+alpha*(2*active_pairs*(2048-128*log2(255))-log2(255)-LIMIT*log2(z))
    return dict(active_pairs=active_pairs,active_groups=2*active_pairs,contracted_runs=runs,
        mark=z,alpha=alpha,marker_log_odds=list(odds),mode=mode,
        log2_upper_screen=prefix+total(run[0]),even_composition_steps=even_steps,
        odd_entry_log2_improvements=(global_odd-odd).tolist(),alias_diagnostics=diagnostics,
        full_under_even_matching_only=True,directed=False,whole_code_certificate=False)


def regional_matrices(active_pairs,runs,z,alpha,odds,*,even_regional=False,even_steps=1):
    """Global coefficient bounds on each parity's regional matrix."""
    full_odds,singleton_odds,mixed_full_odds = odds
    degree,coeff = coefficients(z,alpha,(2,),h=4,route_width=2,steps=even_steps)
    weighted = coeff+degree[:,0,None,None]*full_odds
    full_norm = log2(comb(4096,active_pairs))+active_pairs*full_odds
    if even_regional:
        weighted = weighted[:min(active_pairs,32*even_steps)+1]
        peak = np.max(weighted)
        base = (np.exp2(weighted-peak),peak)
        result = (np.eye(2).reshape((1,2,2)),0.)
        exponent = 128//even_steps
        while exponent:
            if exponent&1:
                result = product(result,base,active_pairs)
            exponent >>= 1
            if exponent:
                base = product(base,base,active_pairs)
        polynomial,scale = result
        if not np.all(polynomial[active_pairs]>0):
            raise ArithmeticError('FLOAT even-region coefficient underflow')
        even = np.log2(polynomial[active_pairs])+scale-full_norm
    else:
        full_step = logsumexp(weighted*log(2),axis=0)/log(2)
        even = np.array(power(tuple(map(tuple,full_step)),128//even_steps))-full_norm
    odd_step = mixed_step(z,alpha,singleton_odds,mixed_full_odds)
    counts = (4096-active_pairs-runs,runs,runs,active_pairs-runs)
    odd_norm = log_multinomial(counts)+2*runs*singleton_odds+(active_pairs-runs)*mixed_full_odds
    odd = np.array(power(tuple(map(tuple,odd_step)),128))-odd_norm
    return even,odd


def bound(active_pairs,runs,z,alpha,odds,*,mode='alternating',even_regional=False,even_steps=1):
    if not 0<active_pairs<4096 or not 1<=runs<=min(active_pairs,4096-active_pairs):
        raise ValueError('Invalid contracted-cycle run profile')
    even,odd = regional_matrices(active_pairs,runs,z,alpha,odds,even_regional=even_regional,even_steps=even_steps)
    if mode=='alternating':
        run = power(multiply(tuple(map(tuple,even)),tuple(map(tuple,odd))),128)
    elif mode=='random':
        region = np.logaddexp2(even,odd)-1
        run = power(tuple(map(tuple,region)),256)
    else:
        raise ValueError('Pairing schedule must be alternating or random')
    support_bits = log2(circular_run_count(4096,active_pairs,runs))
    prefix = support_bits+alpha*(2*active_pairs*(2048-128*log2(255))-log2(255)-LIMIT*log2(z))
    return prefix+total(run[0])


def proposal(active_pairs,runs,*,mode='alternating',iterations=150,even_regional=False,even_steps=1):
    n0 = max(1,4096-active_pairs-runs)
    initial = [max(.001,2*active_pairs/512),.08,
               log2(active_pairs/(4096-active_pairs)),log2(runs/n0),
               log2(max(1,active_pairs-runs)/n0)]
    fixed_full_odds = initial[2]
    if even_regional:
        initial.pop(2)
    best = None
    def objective(theta):
        nonlocal best
        ell,alpha,*odds = theta
        if even_regional:
            odds = [fixed_full_odds]+odds
        if not 1e-5<ell<16 or not .002<alpha<=1 or any(abs(v)>64 for v in odds):
            return 1e12
        value = float(bound(active_pairs,runs,2**-ell,alpha,odds,mode=mode,even_regional=even_regional,even_steps=even_steps))
        if np.isfinite(value) and (best is None or value<best['log2_upper_screen']):
            best = dict(active_pairs=active_pairs,active_groups=2*active_pairs,contracted_runs=runs,
                odd_singletons_each_orientation=runs,odd_full_pairs=active_pairs-runs,
                mark=2**-ell,alpha=alpha,marker_log_odds=odds,mode=mode,
                log2_support_count=log2(circular_run_count(4096,active_pairs,runs)),
                log2_upper_screen=value,full_under_even_matching_only=True,
                coefficient_method='global-three-type-markers-with-equal-singleton-markers',
                even_regional=even_regional,even_composition_steps=even_steps,
                directed=False,whole_code_certificate=False)
        return value if np.isfinite(value) else 1e12
    simplex = np.tile(initial,(len(initial)+1,1))
    for coordinate in range(len(initial)):
        simplex[coordinate+1,coordinate] += .5 if coordinate>=2 else .05*initial[coordinate]
    minimize(objective,initial,method='Nelder-Mead',options={'maxiter':iterations,'initial_simplex':simplex})
    if best is None:
        raise ArithmeticError('No finite shifted-pair FLOAT screen')
    return best


if __name__=='__main__':
    import argparse,json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active-pairs',type=int,default=512)
    parser.add_argument('--runs',type=int,nargs='+',default=(128,256,448))
    parser.add_argument('--mode',choices=('alternating','random'),default='alternating')
    parser.add_argument('--iterations',type=int,default=150)
    parser.add_argument('--even-regional',action='store_true')
    parser.add_argument('--even-steps',type=int,choices=(1,2),default=1)
    args = parser.parse_args()
    for runs in args.runs:
        print(json.dumps(proposal(args.active_pairs,runs,mode=args.mode,iterations=args.iterations,even_regional=args.even_regional,even_steps=args.even_steps)),flush=True)

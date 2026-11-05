"""FLOAT bounds for independently permuted four-group source packets.

The inner input is X[16*lane+j] from lane `lane` of routed source packet j.
Potential source groups are classified by the number (1..4) in each packet.
Local lane occupancies use an entrywise envelope over exact Gale--Ryser
constraints. This is a new route bound, not the scalar-route certificate.
Search output is not directed arithmetic and is never a whole-code receipt.
"""
from functools import lru_cache
from itertools import combinations_with_replacement, product as cartesian
from math import comb, factorial, lgamma, log, log2

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from hd_common_refresh_screen import principal_defects
from hd_grouped_quiet import quiet_screen
from hd_randomized_outer_state_profile import power, total
from hd_sumrank_regions import product

PACKETS = 2048
REGIONS = 256
STEPS = 128
LIMIT = (1 << 21)*7//20
Q = 256.
OCC = np.array(tuple(combinations_with_replacement(range(17), 4)), dtype=int)
DEGREE = OCC.sum(axis=1)
LARGEST = np.cumsum(OCC[:, ::-1], axis=1)


def log_multinomial(counts):
    return (lgamma(sum(counts)+1)-sum(lgamma(n+1) for n in counts))/log(2)


@lru_cache(None)
def return_bernstein(h):
    """Cache the z-independent positive coefficients of the return cap."""
    length = 4*h
    nu = principal_defects(h)
    coefficients = np.zeros((length+1,length+1))
    for l in range(length+1):
        rho = min(h,l)
        for c in range(l+1):
            for e in range(min(h,length-l+1)):
                rank_cap = Q**(-max(0,c-min(nu[c],max(0,h-l+c))))
                count_cap = Q**(max(0,l-c-h)-l)*(Q**min(rho,h-e)-1)
                probability = min(1-Q**(-rho),rank_cap,count_cap)
                coefficients[l,c+e] += comb(l,c)*comb(length-l,e)*probability
    return coefficients


@lru_cache(None)
def size_profiles(sizes):
    """Return step size counts, ordering counts, and feasible lane profiles."""
    result = []
    for counts in cartesian(range(17), repeat=len(sizes)):
        if sum(counts)>16:
            continue
        scalar_count = sum(s*n for s,n in zip(sizes, counts))
        caps = np.array([sum(min(s,k)*n for s,n in zip(sizes, counts))
                         for k in range(1,5)])
        feasible = np.flatnonzero((DEGREE == scalar_count)
                                 & np.all(LARGEST <= caps, axis=1))
        if not len(feasible):
            raise ArithmeticError('Empty graphical degree sequence')
        full = (16-sum(counts),)+counts
        multiplicity = factorial(16)
        for n in full:
            multiplicity //= factorial(n)
        result.append((counts, multiplicity, feasible))
    return tuple(result)


def local_logs(z,h=4,occupancies=None):
    """Common-refresh operator; h8/h16 use the conservative nullity <= h."""
    if h not in (4,8,16):
        raise ValueError('Supported local ranks are 4, 8, and 16')
    occ = OCC if occupancies is None else np.asarray(occupancies)
    length,blocks = 4*h,16//h
    if occ.shape[1]!=blocks or np.any(occ<0) or np.any(occ>length):
        raise ValueError('Occupancies do not fit the emission blocks')
    degree = occ.sum(axis=1)
    quiet = np.array([sum(comb(j,c)*z**(j-c)*(1-z)**c*Q**(-c-min(h,j-c))
                          for c in range(j+1)) for j in range(length+1)])
    k = np.prod(quiet[occ], axis=1)
    t = (z+(1-z)/Q)**degree
    base = Q**-h
    excess = np.array([sum(comb(n,c)*z**(n-c)*(1-z)**c*(Q**(-c)-base)
                           for c in range(min(h,n+1))) for n in range(length+1)])
    running = np.ones(len(occ))
    nonzero = np.zeros(len(occ))
    for b in range(blocks):
        e = excess[length-occ[:,b]]
        nonzero = nonzero*base+running*e
        running *= base+e
    v = t*nonzero/(1-Q**-16)
    local_t = (z+(1-z)/Q)**np.arange(length+1)
    local_r = np.minimum(np.maximum(0,local_t-quiet),
        (1-Q**(-np.minimum(h,np.arange(length+1)).astype(float)))
        *z**np.maximum(0,length-h+1-np.arange(length+1)))
    degrees = np.arange(length+1)
    local_r = np.minimum(local_r,return_bernstein(h)@(z**(length-degrees)*(1-z)**degrees))
    running = np.ones(len(occ))
    ret = np.zeros(len(occ))
    for b in range(blocks):
        kb,jb = quiet[occ[:,b]],local_r[occ[:,b]]
        ret = ret*kb+running*jb
        running *= kb+jb
    ret = np.minimum(ret,(Q**16-1)*v)
    density = (Q**16-1)/(Q-1)**16
    entries = np.stack((k,density*np.maximum(0,t-k),ret/(Q**16-1),density*v),axis=1)
    entries[degree==0,3] = np.minimum(entries[degree==0,3],quiet_screen(z,q=256,state=16,width=64,h=h))
    with np.errstate(divide='ignore'):
        return np.log2(entries.reshape((-1,2,2)))


def coefficients(z,alpha,sizes,h=4,steps=1,route_width=4,exact_four=False):
    if exact_four:
        if steps!=4 or route_width!=2 or sizes!=(2,):
            raise ValueError('Four-step profile composition requires full pairs')
        degree = np.arange(33)
        weights = np.array([log2(comb(32,j)) for j in degree])
        if h==16:
            logs = local_logs(z,h,2*degree[:,None])
        else:
            # Jensen conditional on each step's full-pair count. Retain the
            # four individual step counts until after matrix composition.
            _,single = coefficients(z,1.,sizes,h,1,route_width)
            logs = single-weights[:,None,None]
        pairs = np.logaddexp2(logs[:,None,:,0,None]+logs[None,:,None,0,:],
                             logs[:,None,:,1,None]+logs[None,:,None,1,:]).reshape((-1,2,2))
        pair_weights = (weights[:,None]+weights[None,:]).ravel()
        pair_degree = (degree[:,None]+degree[None,:]).ravel()
        full_degree = (pair_degree[:,None]+pair_degree[None,:]).ravel()
        full_weights = pair_weights[:,None]+pair_weights[None,:]
        result = np.full((129,2,2),-np.inf)
        for r in range(2):
            for c in range(2):
                entry = np.logaddexp2(pairs[:,r,0,None]+pairs[None,:,0,c],
                                     pairs[:,r,1,None]+pairs[None,:,1,c])
                exponent = full_weights+alpha*entry
                peak = np.max(exponent)
                value = np.bincount(full_degree,weights=np.exp2(exponent-peak).ravel(),minlength=129)
                with np.errstate(divide='ignore'):
                    result[:,r,c] = np.log2(value)+peak
        return np.arange(129)[:,None],result
    if steps in (4,8,16):
        if route_width!=2 or sizes!=(2,) or h!=16:
            raise ValueError('Larger composed blocks currently require full pairs and h16')
        _,single = coefficients(z,1.,sizes,h,1,route_width)
        peak = np.max(single)
        base = (np.exp2(single-peak),peak)
        result = (np.eye(2).reshape((1,2,2)),0.)
        exponent = steps
        while exponent:
            if exponent&1:
                result = product(result,base,32*steps)
            exponent >>= 1
            if exponent:
                base = product(base,base,32*steps)
        polynomial,scale = result
        weights = np.array([log2(comb(32*steps,j)) for j in range(32*steps+1)])
        with np.errstate(divide='ignore'):
            value = (1-alpha)*weights[:,None,None]+alpha*(np.log2(polynomial)+scale)
        return np.arange(32*steps+1)[:,None],value
    if route_width==2:
        if sizes!=(2,) and h==16 and steps in (1,2):
            if any(s not in (1,2) for s in sizes):
                raise ValueError('Source pairs have sizes one or two')
            counts = np.array([row for row in cartesian(range(33),repeat=len(sizes)) if sum(row)<=32])
            degree = counts@np.array(sizes)
            logs = local_logs(z,16,np.arange(65)[:,None])[degree]
            weights = np.array([log_multinomial((32-sum(row),)+tuple(row)) for row in counts])
            if steps==1:
                return counts,weights[:,None,None]+alpha*logs
            pair = np.logaddexp2(logs[:,None,:,0,None]+logs[None,:,None,0,:],
                                logs[:,None,:,1,None]+logs[None,:,None,1,:])
            weighted = weights[:,None,None,None]+weights[None,:,None,None]+alpha*pair
            degree = counts[:,None,:]+counts[None,:,:]
            if len(sizes)==1:
                result = np.array([logsumexp(weighted[degree[:,:,0]==j]*log(2),axis=0)/log(2)
                                   for j in range(65)])
                return np.arange(65)[:,None],result
            return degree.reshape((-1,len(sizes))),weighted.reshape((-1,2,2))
        if sizes!=(2,) or h not in (4,8,16) or steps not in (1,2):
            raise ValueError('Two-group routing currently requires full pairs')
        if h==16:
            degree = np.arange(33)
            logs = local_logs(z,h,2*degree[:,None])
            weights = np.array([log2(comb(32,j)) for j in degree])
        else:
            halves = np.array(tuple(combinations_with_replacement(range(17),2)))
            occ = np.repeat(halves,2,axis=1) if h==4 else 2*halves
            logs = local_logs(z,h,occ)
            weights = np.array([log2(comb(16,a)*comb(16,b)*(1 if a==b else 2)) for a,b in halves])
            degree = halves.sum(axis=1)
        if steps==2:
            pair = np.logaddexp2(logs[:,None,:,0,None]+logs[None,:,None,0,:],
                                logs[:,None,:,1,None]+logs[None,:,None,1,:])
            weighted = weights[:,None,None,None]+weights[None,:,None,None]+alpha*pair
            degree = degree[:,None]+degree[None,:]
        else:
            weighted = weights[:,None,None]+alpha*logs
        result = np.array([logsumexp(weighted[degree==j]*log(2),axis=0)/log(2)
                           for j in range(32*steps+1)])
        return np.arange(32*steps+1)[:,None],result
    if route_width!=4:
        raise ValueError('Supported source packet widths are two and four')
    if steps==2:
        if len(sizes)!=1:
            raise ValueError('Paired composition currently requires a homogeneous size class')
        _,single = coefficients(z,1.,sizes,h)
        weights = np.array([log2(comb(16,j)) for j in range(17)])
        single -= weights[:,None,None]
        pair = np.logaddexp2(single[:,None,:,0,None]+single[None,:,None,0,:],
                            single[:,None,:,1,None]+single[None,:,None,1,:])
        weighted = weights[:,None,None,None]+weights[None,:,None,None]+alpha*pair
        degrees = np.arange(17)[:,None]+np.arange(17)[None,:]
        result = np.full((33,2,2),-np.inf)
        for degree in range(33):
            result[degree] = logsumexp(weighted[degrees==degree]*log(2),axis=0)/log(2)
        return np.arange(33)[:,None],result
    if steps!=1:
        raise ValueError('Supported composition lengths are 1 and 2')
    if h!=4:
        if sizes!=(4,):
            raise ValueError('h8/h16 currently screen fully occupied source packets only')
        occ = np.repeat((np.arange(17)*(h//4))[:,None],16//h,axis=1)
        logs = local_logs(z,h,occ)
        return np.arange(17)[:,None],np.array([log2(comb(16,j)) for j in range(17)])[:,None,None]+alpha*logs
    logs = local_logs(z)
    profiles = size_profiles(sizes)
    counts = np.array([row for row,_,_ in profiles])
    coefficient = np.array([log2(n)+alpha*np.max(logs[indices],axis=0)
                            for _,n,indices in profiles])
    return counts,coefficient


def bound(histogram,z,alpha,marker_logs,*,regional=False,h=4,steps=1,route_width=4,exact_four=False):
    """Bound all active-group supports with the specified n1..n4 histogram."""
    histogram = tuple(histogram)
    source_packets = PACKETS*4//route_width
    if len(histogram)!=4 or min(histogram)<0 or not 0<sum(histogram)<=source_packets or any(histogram[route_width:]):
        raise ValueError('Invalid source packet histogram')
    sizes = tuple(s for s in range(1,5) if histogram[s-1])
    n = np.array([histogram[s-1] for s in sizes])
    counts,coeff = coefficients(z,alpha,sizes,h,steps,route_width,exact_four)
    active_groups = sum(s*n for s,n in enumerate(histogram,1))
    log_assignments = log_multinomial((source_packets-sum(histogram),)+histogram)
    log_supports = log_assignments+sum(n*log2(comb(route_width,s)) for s,n in enumerate(histogram,1) if n)
    prefix = log_supports+alpha*(active_groups*(2048-128*log2(255))
                                -log2(255)-LIMIT*log2(z))
    weighted = coeff+(counts@np.asarray(marker_logs))[:,None,None]
    if not regional:
        matrix = logsumexp(weighted*log(2),axis=0)/log(2)
        return prefix+total(power(tuple(map(tuple,matrix)),STEPS*REGIONS//steps)[0])\
            -REGIONS*(log_assignments+n@marker_logs)
    if len(sizes)!=1:
        raise ValueError('Regional extraction currently requires a homogeneous size class')
    target = int(n[0])
    degree = min(target,(64//route_width)*steps)
    coeff = weighted[:degree+1]
    peak = np.max(coeff)
    base = (np.exp2(coeff-peak),peak)
    run = (np.eye(2).reshape((1,2,2)),0.)
    exponent = STEPS//steps
    while exponent:
        if exponent&1:
            run = product(run,base,target)
        exponent >>= 1
        if exponent:
            base = product(base,base,target)
    polynomial,scale = run
    if not np.all(polynomial[target]>0):
        raise ArithmeticError('Regional FLOAT underflow')
    region = np.log2(polynomial[target])+scale-log_assignments-n@marker_logs
    return prefix+total(power(tuple(map(tuple,region)),REGIONS)[0])


def proposal(histogram,*,regional=False,iterations=80,robust=False,h=4,steps=1,route_width=4,exact_four=False,initial_witness=None):
    histogram = tuple(histogram)
    nonempty = sum(histogram)
    sizes = tuple(s for s in range(1,5) if histogram[s-1])
    denominator = max(1,PACKETS*4//route_width-nonempty)
    initial_odds = [log2(histogram[s-1]/denominator) for s in sizes]
    best = None

    def objective(theta):
        nonlocal best
        ell,alpha,*odds = theta
        if regional:
            odds = initial_odds
        if not 1e-5<=ell<=16 or not .002<=alpha<=1 or any(abs(x)>64 for x in odds):
            return 1e12
        try:
            value = float(bound(histogram,2**-ell,alpha,np.array(odds),regional=regional,h=h,steps=steps,route_width=route_width,exact_four=exact_four))
        except ArithmeticError:
            return 1e12
        if not np.isfinite(value):
            return 1e12
        if best is None or value<best['log2_upper_screen']:
            best = dict(histogram=list(histogram),active_groups=sum(s*n for s,n in enumerate(histogram,1)),
                mark=2**-ell,alpha=alpha,marker_log_odds=odds,regional=regional,
                local_rank=h,composition_steps=steps,source_packet_width=route_width,
                internal_profile_average=steps>2 and not exact_four,exact_four=exact_four,
                conditional_half_profile_average=exact_four and h!=16,
                log2_upper_screen=value,directed=False,whole_code_certificate=False)
        return value

    alphas = (.05,.125,.25,.5,.9) if robust else ((.125,1.) if nonempty<=16 else (.0625,))
    factors = (.0625,.25,1.) if robust else (1.,)
    if initial_witness is not None:
        if robust:
            raise ValueError('An explicit start and robust multistart are separate modes')
        alphas,factors = (initial_witness[1],),(1.,)
    for alpha in alphas:
        for factor in factors:
            initial = [max(.001,sum(s*n for s,n in enumerate(histogram,1))/512*factor),alpha]
            if initial_witness is not None:
                initial[:2] = initial_witness[:2]
            if not regional:
                initial += initial_odds
            options = {'maxiter':iterations}
            if not regional:
                simplex = np.tile(initial,(len(initial)+1,1))
                for coordinate in range(len(initial)):
                    simplex[coordinate+1,coordinate] += (.5 if coordinate>=2 else max(.001,.05*abs(initial[coordinate])))
                options['initial_simplex'] = simplex
            minimize(objective,initial,method='Nelder-Mead',options=options)
    if best is None:
        raise ArithmeticError('No usable FLOAT candidate')
    return best


if __name__=='__main__':
    import argparse,json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--histogram',type=int,nargs=4,required=True,metavar=('N1','N2','N3','N4'))
    parser.add_argument('--regional',action='store_true')
    parser.add_argument('--iterations',type=int,default=80)
    parser.add_argument('--robust',action='store_true',help='Use fifteen dispersed search starts')
    parser.add_argument('--rank',type=int,choices=(4,8,16),default=4)
    parser.add_argument('--steps',type=int,choices=(1,2,4,8,16),default=1)
    parser.add_argument('--route-width',type=int,choices=(2,4),default=4)
    parser.add_argument('--exact-four',action='store_true')
    parser.add_argument('--initial-witness',type=float,nargs=2,metavar=('NEG_LOG2_MARK','ALPHA'))
    args = parser.parse_args()
    print(json.dumps(proposal(args.histogram,regional=args.regional,iterations=args.iterations,robust=args.robust,h=args.rank,steps=args.steps,route_width=args.route_width,exact_four=args.exact_four,initial_witness=args.initial_witness)),flush=True)

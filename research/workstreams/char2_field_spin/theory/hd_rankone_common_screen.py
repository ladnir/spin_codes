"""FLOAT gate for rank-one refresh of Q+CX, not M Q+CX.

The unchanged packet map emits X+A Q, with block RS[16,4] A and C=A^T.
Only the refresh changes: M=I+u v^T, u uniform nonzero and v uniform
conditional on v^T u != 1. Components are zero, arbitrary-live mass,
and uniform-live density. This is neither a directed nor whole-code proof.
"""
from functools import lru_cache
from math import comb, log2

import numpy as np
from scipy.special import logsumexp

from hd_common_refresh_screen import local_return_cap
from hd_shifted_swap_sparse import sparse_geometry
from hd_shifted_outer_density import log_density_float, TOEPLITZ_OUTER_MEASURE
from hd_grouped_packet_screen import LIMIT, log_multinomial


def packet_matrix(occupancies, z, *, q=256, h=4, flat_mode='density'):
    """A positive three-component cap for each fixed routed support.

    Random X is uniform on its potential support, including zero values.
    Profile bounds apply before route averaging and fractional powers.
    In particular no diagonal-Cauchy empty-packet cap is reused.
    """
    occ = tuple(occupancies)
    if not occ or h not in (1, 2, 4) or q < 2 or not 0 < z <= 1:
        raise ValueError('Invalid field, local rank, or weight mark')
    if any(type(v) is not int or not 0 <= v <= 4*h for v in occ):
        raise ValueError('Invalid occupancy')
    if flat_mode not in ('density','mass'):
        raise ValueError('Unknown flat-component decomposition')
    state, length, degree = h*len(occ), 4*h, sum(occ)
    quiet, excess, returns = local_values(z, q, h)
    t = z+(1-z)/q
    K, total, J = 1., 1., 0.
    live_total, live_nonzero = 1., 0.
    base = float(q)**-h
    for value in occ:
        k, j = quiet[value], returns[value]
        J = J*k+total*j
        total *= k+j
        K *= k
        e = excess[length-value]
        live_nonzero = live_nonzero*base+live_total*e
        live_total *= base+e
    d = float(q)**state-1
    T = t**degree
    V = T*live_nonzero/(1-float(q)**-state)
    J = min(J, d*V)
    rho = sum(min(h, value) for value in occ)
    delta = max(0, length-h+1-max(occ))
    arbitrary = T*z**delta
    returned = (float(q)**-rho*min(t**(degree-rho), z**max(1, delta))) if degree else 0.
    pointwise = min(T, z**delta)
    eta, gamma = 1/q, 1-1/q
    born = max(0., T-K)
    flat_row = (J/d,0.,eta*pointwise+gamma*V) if flat_mode=='density' else (J/d,eta*V,gamma*V)
    return np.array(((K, eta*born, gamma*born),
                     (returned, eta*arbitrary, gamma*arbitrary),flat_row))


@lru_cache(32)
def local_values(z, q, h):
    length = 4*h
    quiet = tuple(sum(comb(j,c)*z**(j-c)*(1-z)**c*float(q)**(-c-min(h,j-c))
                      for c in range(j+1)) for j in range(length+1))
    base = float(q)**-h
    excess = tuple(sum(comb(n,c)*z**(n-c)*(1-z)**c*(float(q)**-c-base)
                       for c in range(min(h,n+1))) for n in range(length+1))
    t = z+(1-z)/q
    returns = tuple(min(max(0., t**j-quiet[j]),
                        (1-float(q)**-min(h,j))*z**max(0,length-h+1-j),
                        local_return_cap(j,z,h,length,float(q))) for j in range(length+1))
    return quiet, excess, returns


def log_product(left, right):
    return logsumexp((left[:,:,None]+right[None,:,:])*np.log(2), axis=1)/np.log(2)


def log_power(matrix, exponent):
    result = np.full((3,3), -np.inf)
    np.fill_diagonal(result, 0.)
    while exponent:
        if exponent&1:
            result = log_product(result, matrix)
        exponent >>= 1
        if exponent:
            matrix = log_product(matrix, matrix)
    return result


def pair_coefficients(z, alpha, maximum=64, *, flat_mode='density'):
    if not 0 < alpha <= 1 or type(maximum) is not int or not 0 <= maximum <= 64:
        raise ValueError('Invalid fractional exponent or packet cutoff')
    packet = {}
    result = {}
    for s, f, count, occ in sparse_geometry(maximum):
        if occ not in packet:
            packet[occ] = packet_matrix(occ, z, flat_mode=flat_mode)**alpha
        result.setdefault((s,f), np.zeros((3,3)))[:] += (count/2**s)*packet[occ]
    with np.errstate(divide='ignore'):
        return {key: np.log2(value) for key,value in result.items()}


def polynomial_product(left, right, maximum):
    result = {}
    for (s,f), a in left.items():
        for (t,g), b in right.items():
            key = s+t,f+g
            if key[0]+2*key[1] > maximum:
                continue
            value = log_product(a,b)
            result[key] = np.logaddexp2(result[key], value) if key in result else value
    return result


def exact_sparse_regions(z, alpha, maximum=2):
    """FLOAT positive extraction, hard-limited to the initial sparse gate."""
    if maximum not in (1,2):
        raise ValueError('Sparse gate is limited to weights one and two')
    base = pair_coefficients(z,alpha,maximum)
    identity = np.full((3,3),-np.inf)
    np.fill_diagonal(identity,0.)
    result, exponent = {(0,0):identity}, 128
    while exponent:
        if exponent&1:
            result = polynomial_product(result,base,maximum)
        exponent >>= 1
        if exponent:
            base = polynomial_product(base,base,maximum)
    return {(s+2*f,f): matrix-log2(comb(4096,s)*comb(4096-s,f))
            for (s,f),matrix in result.items()}


def profile_bound(weight, a, b, count, z, alpha, regions):
    prefix = log2(count)+alpha*(weight*log_density_float(TOEPLITZ_OUTER_MEASURE)
                                -log2(255)-LIMIT*log2(z))
    end = log_power(log_product(regions[weight,a],regions[weight,b]),128)[0]
    return prefix+float(logsumexp(end*np.log(2))/np.log(2))


def sparse_gate(weight, z, alpha=1.):
    from hd_xor_group_geometry import matching_joint_counts
    if weight not in (1,2):
        raise ValueError('Sparse gate is limited to weights one and two')
    regions = exact_sparse_regions(z,alpha,weight)
    rows = [dict(a=a,b=b,count=count,log2_bound=profile_bound(weight,a,b,count,z,alpha,regions))
            for (a,b),count in matching_joint_counts(8192,weight).items()]
    union = float(logsumexp(np.array([row['log2_bound'] for row in rows])*np.log(2))/np.log(2))
    return dict(active_groups=weight,z=z,alpha=alpha,rows=rows,log2_union_FLOAT=union,
                directed=False,whole_code_certificate=False,refresh='rank-one-after-feedback')


def arbitrary_path_lower_log(weight, z, alpha, *, source_count=1):
    """Lower bound on THIS majorant, never a lower bound on failure.

    A first active packet births the arbitrary component; its self-loop
    can then charge only thirteen additional emitted coordinates/packet.
    The event that the first routed potential coordinate is active has
    probability weight/8192 under every fixed source matching profile.
    """
    if (type(weight) is not int or not 1 <= weight <= 8192 or not 0 < z <= 1
            or not 0 < alpha <= 1 or type(source_count) is not int or source_count < 1):
        raise ValueError('Invalid path-lower-bound parameters')
    packets = 32768
    return (log2(source_count)+log2(weight/8192)+alpha*(
        weight*log_density_float(TOEPLITZ_OUTER_MEASURE)-log2(255)
        -8*packets+log2(255/256)
        +(256*weight+13*(packets-1)-LIMIT)*log2(z)))


def symmetric_positive_gate(weight, full, z, alpha, *, max_evaluations=100, flat_mode='density'):
    """One balanced matching class; optimize only two positive markers.

    Single-matching count upper bounds the indexed supports in the joint
    class. Every packet uses the actual independent fair-swap pair law.
    """
    from scipy.optimize import minimize
    from hd_shifted_pairs_count import one_pairing_count
    if type(max_evaluations) is not int or not 1 <= max_evaluations <= 200:
        raise ValueError('An explicit budget in1..200 is required')
    histogram = (4096-weight+full,weight-2*full,full)
    if min(histogram) <= 0:
        raise ValueError('This diagnostic requires an interior symmetric histogram')
    coeff = pair_coefficients(z,alpha,flat_mode=flat_mode)
    keys = np.array(tuple(coeff))
    logs = np.array(tuple(coeff.values()))
    count = one_pairing_count(4096,weight,full)
    prefix = log2(count)+alpha*(weight*log_density_float(TOEPLITZ_OUTER_MEASURE)
                                -log2(255)-LIMIT*log2(z))-256*log_multinomial(histogram)
    def objective(odds):
        matrix = logsumexp((logs+np.dot(keys,odds)[:,None,None])*np.log(2),axis=0)/np.log(2)
        end = log_power(matrix,32768)[0]
        return prefix+float(logsumexp(end*np.log(2))/np.log(2))-256*np.dot(histogram[1:],odds)
    initial = np.log2(np.array(histogram[1:])/histogram[0])
    fit = minimize(objective,initial,method='Nelder-Mead',options={'maxfev':max_evaluations})
    return dict(active_groups=weight,a=full,b=full,z=z,alpha=alpha,
                marker_log2=fit.x.tolist(),log2_bound_FLOAT=float(fit.fun),
                evaluations=fit.nfev,count_method='single-matching-upper',directed=False,
                whole_code_certificate=False,refresh='rank-one-after-feedback',flat_mode=flat_mode)

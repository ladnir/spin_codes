"""Bounded FLOAT active-row screens for smaller GRS outers, state8/width16.

The packet cover is exactly the stronger positive joint-return cover used
by state8_density_certificate.py. Regional placement averaging retains the
finite without-replacement normalization. Sampled active counts are not a
complete sparse certificate, and a failed upper bound does not reject a code.
No native implementation or runtime measurement is performed here.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal
from state8_joint_return import joint_return_bound, minimum_return_weight
from superregular_return_screen import local_matrix
from rate_design_coset import mul, power


ACTIVE_SAMPLES = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512)

# FLOAT witness seeds only; these are not interval or probability certificates.
WITNESSES_128_87 = {
    1:(6.582489861693686,5.541263545158426),
    2:(6.5757395035101665,5.541263545158426),
    4:(6.5723381159152545,5.541263545158426),
    8:(6.7446187347929625,5.320554324478209),
    16:(6.505363389863845,2.9768305126928345),
    32:(6.439330040564994,2.2792591335024888),
    64:(6.377801453275847,2.003847457724925),
    96:(6.3468322152020455,1.93374094857465),
    128:(6.332834733713911,1.9037272141625097),
    160:(6.325998197527945,1.8892625599192094),
    192:(6.32276022733745,1.8825577098849262),
    256:(6.32168880233887,1.8807830309646731),
    512:(6.339939394192531,1.9230030131189555),
}


@lru_cache(maxsize=None)
def shortening_coefficients(dimension=87):
    """Value at point k from source evaluations at points0,...,k-1."""
    if not 1 <= dimension < 256:
        raise ValueError('invalid consecutive-point interpolation dimension')
    result = []
    for j in range(dimension):
        numerator = denominator = 1
        for h in range(dimension):
            if h != j:
                numerator = mul(numerator, dimension ^ h)
                denominator = mul(denominator, j ^ h)
        result.append(mul(numerator,power(denominator,254)))
    return tuple(result)


def shortening_cost(dimension=87):
    """Exact scalar coefficient-basis count, not an instruction schedule."""
    coefficients = shortening_coefficients(dimension)
    basis, span = [1], {0,1}
    for value in sorted(set(coefficients)):
        if value not in span:
            span |= {x ^ value for x in tuple(span)}
            basis.append(value)
    return dict(dimension=dimension,missing_point=dimension,
                basis_products=len(basis)-1,basis=basis,
                coefficient_multiplicities=dict(sorted(Counter(coefficients).items())),
                runtime_claim=False)


def packet(weight, y, tau, joint=True):
    result = local_matrix(weight, y, tau, state=8, width=16)
    if weight:
        if joint:
            result[1][0] = joint_return_bound(weight, y, tau, 16)
        else:
            # Match the improved pure return, rather than the obsolete 9-w.
            old = 9 if weight == 1 else max(1, 9-weight)
            result[1][0] *= y**(minimum_return_weight(weight, 16)-old)
    return result


def log_coefficients(active, y, tau, rows, joint=True):
    if not 0 <= active <= rows or rows <= 0 or rows % 16:
        raise ValueError('invalid packet-aligned region geometry')
    local = np.asarray([packet(w, y, tau, joint) for w in range(min(16, active)+1)], dtype=float)
    with np.errstate(divide='ignore'):
        factor = np.log(local)
    for w in range(len(local)):
        factor[w] += math.log(math.comb(16, w))
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    packets = rows//16
    while packets:
        if packets & 1:
            result = log_polynomial_product(result, factor, active)
        packets >>= 1
        if packets:
            factor = log_polynomial_product(factor, factor, active)
    return result


def bound(active, point, k_log=20, length=128, dimension=88, joint=True):
    g = geometry(k_log, length, dimension)
    if not 1 <= active <= g['rows'] or not 0 < dimension < length <= 256:
        raise ValueError('invalid outer or active-row count')
    logtilt, logtau = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    tau = min(255., math.exp(logtau))
    coefficients = log_coefficients(active, math.exp(logy), tau, g['rows'], joint)
    regional = regional_from_coefficients(coefficients, active, g['rows'], math.log(tau))
    return (log_choose(g['rows'], active)-((length-dimension)*active+1)*math.log(tau)
            -g['threshold']*logy+log_terminal(regional, length))/math.log(2)


def search(active, k_log=20, length=128, dimension=88, iterations=45, start=None, joint=True):
    first = (math.log(5*length), math.log(255)) if start is None else tuple(start)
    starts = (first, (math.log(5*length), math.log(16)))
    options = [minimize(lambda point: bound(active,point,k_log,length,dimension,joint), point,
                        method='L-BFGS-B', bounds=[(-2,10),(0,math.log(255))],
                        options={'maxiter':iterations,'ftol':1e-10,'maxls':16})
               for point in starts]
    best = min(options, key=lambda value: value.fun)
    g = geometry(k_log, length, dimension)
    coarse = bound(active, best.x, k_log, length, dimension, False)
    return dict(outer=[length,dimension],k_log=k_log,rows=g['rows'],n=g['n'],
                output_threshold=g['threshold'],active=active,state=8,width=16,
                margin_bits=-float(best.fun),logtilt=float(best.x[0]),logtau=float(best.x[1]),
                pure_margin_same_marks=-coarse,joint_gain_bits=coarse-float(best.fun),
                optimizer_success=bool(best.success),iterations=int(best.nit),
                joint_return=joint,rigorous=False,whole_code=False,
                scope='one complete active-row class, evaluated in floating point')


def generic_six_density_bound(point, weight, k_log=20, length=128, dimension=86):
    """Finite global comparison with ONLY the generic state-six packet law."""
    g = geometry(k_log,length,dimension)
    if not 0 < weight < g['n']:
        raise ValueError('input weight must lie strictly inside encoded length')
    alpha = weight/g['n']
    logx, logy, logtau, logtilt = map(float,point)
    y, tau = math.exp(logy), min(255.,math.exp(logtau))
    values = np.asarray([local_matrix(w,y,tau,state=6,width=16) for w in range(17)])
    with np.errstate(divide='ignore'):
        logs = np.log(values)
    for w in range(17):
        logs[w] += math.log(math.comb(16,w))+w*logx
    logarithms = logsumexp(logs,axis=0)
    row = float(logsumexp([length*math.log1p(-alpha),
                          -(length-dimension)*logtau+length*np.logaddexp(
                              math.log1p(-alpha),logtau+logtilt+math.log(alpha))]))
    shift = float(logarithms.max())
    packets = g['n']//16
    inner = packets*shift+log_terminal_power(np.exp(logarithms-shift),packets)
    return (length*math.log(g['rows']+1)+g['rows']*row+inner
            -weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def generic_six_density_search(alpha, k_log=20, length=128, dimension=86, iterations=45, start=None):
    g = geometry(k_log,length,dimension)
    weight = max(1,round(alpha*g['n']))
    starts = [start] if start is not None else [(-2,-.2,math.log(16),0),
                                               (0,-1,math.log(255),-1)]
    candidates = [minimize(generic_six_density_bound,point,args=(weight,k_log,length,dimension),
                           method='L-BFGS-B',bounds=[(-30,12),(-25,0),(0,math.log(255)),(-20,20)],
                           options={'maxiter':iterations,'ftol':1e-11,'maxls':20}) for point in starts]
    best = min(candidates,key=lambda item:item.fun)
    return dict(outer=[length,dimension],k_log=k_log,width=16,state=6,input_weight=weight,
                input_density=weight/g['n'],margin_bits=-float(best.fun),
                log_parameters=best.x.tolist(),iterations=int(best.nit),
                optimizer_success=bool(best.success),rigorous=False,whole_code=False,
                local_cover='generic coarse superregular; no state8 law',
                finite_factors_included=True,scope='one exact input-weight class')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', type=int, nargs=2, default=(128,88))
    parser.add_argument('--active', type=int, nargs='+', default=ACTIVE_SAMPLES)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--iterations', type=int, default=45)
    parser.add_argument('--start', type=float, nargs=2)
    parser.add_argument('--pure', action='store_true')
    parser.add_argument('--shortening-cost', action='store_true')
    parser.add_argument('--generic-six-density', type=float, nargs='+')
    args = parser.parse_args()
    if args.shortening_cost:
        print(json.dumps(shortening_cost(args.outer[1])),flush=True)
        raise SystemExit(0)
    if args.generic_six_density:
        previous = None
        for alpha in args.generic_six_density:
            result = generic_six_density_search(alpha,args.k_log,*args.outer,args.iterations,previous)
            previous = result['log_parameters']
            print(json.dumps(result),flush=True)
        raise SystemExit(0)
    previous = args.start
    for active in args.active:
        result = search(active,args.k_log,*args.outer,args.iterations,previous,not args.pure)
        previous = (result['logtilt'],result['logtau'])
        print(json.dumps(result),flush=True)

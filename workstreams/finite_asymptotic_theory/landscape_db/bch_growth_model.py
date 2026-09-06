"""Persistent-state Q1 models for explaining BCH size and message growth.

These idealizations exclude cancellations and finite-epoch output fluctuations.
They are explanatory models, not bounds on the actual encoder's failure rate.
"""
import math

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp


def onset_probability(block, weight, delta=.1):
    """Exact late-first-activation probability in the continuous idealization."""
    if not 1 <= weight <= block or not 0 < delta < .5:
        raise ValueError('invalid onset geometry')
    tail = 2*delta*block
    whole = math.floor(tail); fraction = tail-whole
    count = (math.comb(whole, weight) if weight <= whole else 0)
    count += fraction*(math.comb(whole, weight-1) if weight-1 <= whole else 0)
    return count/math.comb(block, weight)


def onset_log_moment(block, weight, scaled_tilt):
    """log E exp(-a Y/L), with Y/L=(H-U)/2, U uniform on [0,1]."""
    if not 1 <= weight <= block or scaled_tilt <= 0:
        raise ValueError('invalid moment geometry')
    h = np.arange(weight, block+1, dtype=float)
    log_mass = np.array([math.log(math.comb(int(v)-1,weight-1))
                         -math.log(math.comb(block,weight)) for v in h])
    v = scaled_tilt/2
    return (math.log(-math.expm1(-v))-math.log(v)
            +float(logsumexp(log_mass-v*(h-1))))


def onset_chernoff(block, weight, delta=.1):
    if (weight-1)/2 >= delta*block:
        return dict(log_upper=-math.inf,scaled_tilt=None)
    def objective(log_a):
        a = math.exp(log_a)
        return a*delta*block+onset_log_moment(block,weight,a)
    fit = minimize_scalar(objective,bounds=(-12,12),method='bounded',options={'xatol':1e-11})
    if not fit.success or min(fit.x+12,12-fit.x)<.01:
        raise ArithmeticError('unresolved onset tilt optimum')
    return dict(log_upper=min(0.,float(fit.fun)),scaled_tilt=math.exp(fit.x))


def spectrum_model(block, counts, delta=.1):
    direct = []; chernoff = []; shells = []
    for weight, count in sorted(counts.items()):
        if not weight or not count:
            continue
        probability = onset_probability(block,weight,delta)
        bound = onset_chernoff(block,weight,delta)
        if probability:
            direct.append(math.log(count)+math.log(probability))
        if math.isfinite(bound['log_upper']):
            chernoff.append(math.log(count)+bound['log_upper'])
        shells.append(dict(weight=weight,count=count,onset_probability=probability,**bound))
    return dict(block_bits=block,delta=delta,
                onset_intercept_bits=-float(logsumexp(direct))/math.log(2) if direct else None,
                chernoff_intercept_bits=-float(logsumexp(chernoff))/math.log(2) if chernoff else None,
                shells=shells)

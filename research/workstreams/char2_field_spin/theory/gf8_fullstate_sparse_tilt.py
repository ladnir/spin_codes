"""Floating sparse support screen with the outer's minimum-weight mark.

For A active MDS rows, every admissible raw input support has W >= d*A.
Multiplying its nonnegative summand by lambda**(W-d*A), lambda >= 1,
therefore gives an upper bound before the usual source relaxation.
This is an additional deterministic support restriction, not randomness.
The local rank mark tau is unchanged; only regional degree coefficients
are multiplied by lambda**h. These screens are not certificates.
"""
from __future__ import annotations

import argparse
import json
import math
import numpy as np
from scipy.optimize import minimize_scalar

from cauchy_firstmoment_screen import geometry, log_choose
from gf8_fullstate_screen import WIDTH, local_matrix
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


FROZEN = {
    (8, 64): (6.254072327017321, 1.7442547256123069),
    (8, 128): (6.211329347215399, 1.6425502027360765),
    (16, 64): (6.272717183550212, 2.030030715348006),
    (16, 128): (6.253583178599222, 1.7385402919814381),
}


def fixed_coefficients(active, point, state=8, dimension=96, k_log=20, length=128):
    if not 0 < dimension < length <= 256:
        raise ValueError('invalid GF256 MDS dimensions')
    g = geometry(k_log, length, dimension)
    if not 0 < active <= g['rows']:
        raise ValueError('invalid active-row count')
    logtilt, logtau = map(float, point)
    if not 0 <= logtau <= math.log(255):
        raise ValueError('rank mark tau must be in [1,255]')
    logy = -math.exp(logtilt)*active/g['n']
    local = np.asarray([local_matrix(w, math.exp(logy), math.exp(logtau), state)
                        for w in range(min(WIDTH, active)+1)])
    with np.errstate(divide='ignore'):
        factor = np.log(local)
    for weight in range(len(factor)):
        factor[weight] += math.log(math.comb(WIDTH, weight))
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    packets = g['rows']//WIDTH
    while packets:
        if packets & 1:
            result = log_polynomial_product(result, factor, active)
        packets >>= 1
        if packets:
            factor = log_polynomial_product(factor, factor, active)
    source = log_choose(g['rows'], active)-((length-dimension)*active+1)*logtau
    return result, g, logtau, source-g['threshold']*logy


def evaluate(coefficients, g, active, logtau, base, loglambda):
    if loglambda < 0:
        raise ValueError('minimum-weight mark lambda must be at least one')
    regional = regional_from_coefficients(coefficients, active, g['rows'], logtau+loglambda)
    minimum = (g['length']-g['dimension']+1)*active
    return (base+log_terminal(regional, g['length'])-minimum*loglambda)/math.log(2)


def screen(active, point, state=8, dimension=96, k_log=20, length=128):
    coefficients, g, logtau, base = fixed_coefficients(active, point, state, dimension, k_log, length)
    objective = lambda mark: evaluate(coefficients, g, active, logtau, base, float(mark))
    control = objective(0.)
    derivative = (objective(1e-5)-control)/1e-5
    result = minimize_scalar(objective, bounds=(0., 12.), method='bounded',
                             options={'xatol': 1e-10})
    loglambda, bound = (float(result.x), float(result.fun)) if result.fun < control else (0., control)
    return dict(outer=[length, dimension], state=state, width=WIDTH, active=active,
                n=g['n'], rows=g['rows'], k_log=k_log, log_parameters=list(point),
                loglambda=loglambda, lambda_value=math.exp(loglambda),
                original_margin_bits=-control, marked_margin_bits=-bound,
                derivative_bits_at_one=derivative,
                minimum_weight=(length-dimension+1)*active,
                rigorous=False, whole_code=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', type=int, nargs='+', default=[8,16])
    parser.add_argument('--active', type=int, nargs='+', default=[64,128])
    args = parser.parse_args()
    for state in args.state:
        for active in args.active:
            print(json.dumps(screen(active, FROZEN[state,active], state)), flush=True)

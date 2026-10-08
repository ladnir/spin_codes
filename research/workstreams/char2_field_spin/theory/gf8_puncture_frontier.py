"""Isolated floating screens for fixed or periodic local puncturing.

These are proof-bound evaluations, not distance certificates or evidence
that a negative bound describes an actual bad code. The local cover is
imported unchanged from gf8_punctured_screen. Periods compose in order;
their matrices are never averaged. No native work is performed here.
"""
from __future__ import annotations

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from gf8_punctured_screen import local_matrices
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


def profile(dimension, schedule=(6,), width=8, state=8, k_log=20):
    schedule = tuple(schedule)
    if (not schedule or any(not 1 <= m <= width for m in schedule)
            or width not in (4, 8, 16, 32) or not 1 <= state <= 32):
        raise ValueError('invalid packet or puncture schedule')
    g = geometry(k_log, 128, dimension)
    if g['rows'] % (width * len(schedule)):
        raise ValueError('schedule must divide each regional packet count')
    output_n = g['n'] * sum(schedule) // (width * len(schedule))
    return dict(outer=[128, dimension], width=width, state=state,
                outputs_period=list(schedule), k_log=k_log, requested_k=2**k_log,
                padded_k=g['k'], rows=g['rows'], input_n=g['n'], output_n=output_n,
                threshold=(output_n+9)//10-1, nominal_rate=g['k']/output_n,
                rigorous=False, whole_code=False)


def period_matrix(point, schedule, width=8, state=8):
    logx, logy, logtau = map(float, point[:3])
    y, tau = math.exp(logy), min(255., math.exp(logtau))
    multiplicity = [math.comb(width, w)*math.exp(w*logx) for w in range(width+1)]
    product, logarithm = np.eye(2), 0.
    for outputs in schedule:
        local = np.asarray(local_matrices(y, tau, width, outputs, state), dtype=float)
        packet = np.einsum('w,wij->ij', multiplicity, local)
        product = product @ packet
        scale = float(product.max())
        product /= scale
        logarithm += math.log(scale)
    return product, logarithm


def middle_bound(point, density, dimension, schedule=(6,), width=8, state=8,
                 k_log=20, asymptotic=False):
    g = profile(dimension, schedule, width, state, k_log)
    weight = max(1, min(g['input_n']-1, round(density*g['input_n'])))
    alpha = density if asymptotic else weight/g['input_n']
    logx, logy, logtau, logtilt = map(float, point)
    matrix, scale = period_matrix(point, schedule, width, state)
    source = source_row_log(alpha, logtau, logtilt, 128, dimension, True)
    if asymptotic:
        a, b, c, d = matrix.ravel()
        radius = ((a+d)+math.sqrt((a-d)**2+4*b*c))/2
        inner = (math.log(radius)+scale)/(width*len(schedule))
        kept_fraction = sum(schedule)/(width*len(schedule))
        return (source/128+inner-alpha*(logx+logtilt)-.1*kept_fraction*logy)/math.log(2)
    periods = g['input_n']//(width*len(schedule))
    inner = log_terminal_power(matrix, periods)+periods*scale
    return (128*modal_log_penalty(g['rows'])+g['rows']*source+inner
            -weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def sparse_bound(point, active, dimension, schedule=(6,), width=8, state=8, k_log=20):
    g = profile(dimension, schedule, width, state, k_log)
    logtilt, logtau = map(float, point)
    y = math.exp(-math.exp(logtilt)*active/g['output_n'])
    tau = min(255., math.exp(logtau))
    identity = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    factor = identity
    for outputs in schedule:
        local = np.asarray(local_matrices(y, tau, width, outputs, state)[:min(width, active)+1], dtype=float)
        with np.errstate(divide='ignore'):
            local = np.log(local)
        for w in range(len(local)):
            local[w] += math.log(math.comb(width, w))
        factor = log_polynomial_product(factor, local, active)
    result = identity
    steps = g['rows']//(width*len(schedule))
    while steps:
        if steps & 1:
            result = log_polynomial_product(result, factor, active)
        steps >>= 1
        if steps:
            factor = log_polynomial_product(factor, factor, active)
    regional = regional_from_coefficients(result, active, g['rows'], logtau)
    return (log_choose(g['rows'], active)-((128-dimension)*active+1)*logtau
            -g['threshold']*math.log(y)+log_terminal(regional, 128))/math.log(2)


def search(value, dimension, schedule=(6,), width=8, state=8, sparse=False,
           iterations=60, start=None, asymptotic=False):
    g = profile(dimension, schedule, width, state)
    if sparse:
        target = int(value)
        objective = lambda p: sparse_bound(p, target, dimension, schedule, width, state)
        starts = [(6.5, 2.), (7., 3.)]
        bounds = [(-2, 10), (0, math.log(255))]
    else:
        target = float(value)
        objective = lambda p: middle_bound(p, target, dimension, schedule, width, state,
                                            asymptotic=asymptotic)
        starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
        bounds = [(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)]
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(objective, p, method='L-BFGS-B', bounds=bounds,
                        options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25}) for p in starts]
    best = min(choices, key=lambda r: r.fun if math.isfinite(r.fun) else math.inf)
    return dict(**g, mode='sparse' if sparse else 'asymptotic-middle' if asymptotic else 'middle',
                target=target, margin_bits=-float(best.fun),
                log_parameters=best.x.tolist(), optimizer_success=bool(best.success),
                iterations=int(best.nit))


def dense_point(density, dimension, schedule=(6,), width=8, state=8, k_log=20,
                ball_exponent=1.27, asymptotic=False):
    """Same all-completions dense bound, with a floating tilt optimizer."""
    g = profile(dimension, schedule, width, state, k_log)
    if not 0 < density <= 1:
        raise ValueError('density must lie in (0,1]')
    alpha = density
    parity = 128-dimension
    empty = 128*math.log1p(-alpha) if alpha < 1 else -math.inf
    def per_input(logtilt):
        occupied = -8*parity*math.log(2)+128*math.log1p((256/255*math.exp(logtilt)-1)*alpha)
        bracket = float(np.logaddexp(empty, occupied))
        return (alpha*logtilt-bracket/128)/math.log(2)
    optimum = minimize_scalar(lambda x: -per_input(x), bounds=(0, 40), method='bounded')
    information = max(per_input(0), -float(optimum.fun))
    kept_fraction = g['output_n']/g['input_n']
    rate_margin = information-ball_exponent*kept_fraction-8*(1-kept_fraction)
    total = rate_margin if asymptotic else g['input_n']*rate_margin-128*math.log2(g['rows']+1)
    return dict(**g, mode='asymptotic-dense' if asymptotic else 'dense', target=density,
                margin_bits=total, log_tilt=float(optimum.x), ball_exponent=ball_exponent)


def dense_cutoff(dimension, schedule=(6,), width=8, state=8, target_bits=40):
    endpoint = dense_point(1., dimension, schedule, width, state)
    if endpoint['margin_bits'] < target_bits:
        return dict(endpoint, approximate_cutoff=None)
    lower, upper = .01, 1.
    for _ in range(50):
        midpoint = (lower+upper)/2
        if dense_point(midpoint, dimension, schedule, width, state)['margin_bits'] >= target_bits:
            upper = midpoint
        else:
            lower = midpoint
    return dict(endpoint, approximate_cutoff=upper)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension', type=int, default=80)
    parser.add_argument('--schedule', type=int, nargs='+', default=[6])
    parser.add_argument('--width', type=int, default=8)
    parser.add_argument('--state', type=int, default=8)
    parser.add_argument('--sparse', type=int, nargs='*', default=[])
    parser.add_argument('--density', type=float, nargs='*', default=[.03, .10, .50, .98])
    parser.add_argument('--iterations', type=int, default=60)
    parser.add_argument('--asymptotic', action='store_true')
    args = parser.parse_args()
    for sparse, values in ((True, args.sparse), (False, args.density)):
        previous = None
        for value in values:
            receipt = search(value, args.dimension, tuple(args.schedule), args.width, args.state,
                             sparse, args.iterations, previous, args.asymptotic)
            previous = receipt['log_parameters']
            print(json.dumps(receipt), flush=True)
    print(json.dumps(dense_cutoff(args.dimension, tuple(args.schedule), args.width, args.state)), flush=True)

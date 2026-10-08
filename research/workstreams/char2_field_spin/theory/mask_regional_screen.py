"""Floating regional rank screen with sixteen actual masks and one pooled state.

Pooled nonzero mass is never redistributed among the sixteen actual masks.
Regional support averaging is normalized without replacement, not Bernoulli.
The numerical output is a screen, not a directed whole-code certificate.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_choose
from cauchy_joint_return import rank_return_moment
from mask_state_rank import local_mask_transfer
from sparse_multirow_zero_set import local_rank_transfer


POOLED = 16
STATES = 17
WITNESSES = {
    4: (5.32479700865, 4.96215263668),
    8: (5.13883585678, 2.88863774319),
    16: (5.07867455120, 2.33096175675),
}


def packet_matrix(weight, y, tau, cutoff=2, phase=0):
    """Uniform common-history bound, averaged over C(16,weight) placements.

    Known masks use the exact-mask rank kernel only up to cutoff. Unsupported
    nonzero transitions forget their mask. P uses only a uniform nonzero
    continuation; a return to zero is its sole exit. The joint return API is
    the selected-zero RANK potential, never its fixed-word counterpart.
    """
    if (cutoff not in (1, 2) or phase not in range(4)
            or not 0 <= weight <= 16 or not math.isfinite(float(y))
            or not math.isfinite(float(tau)) or not 0 < y <= 1
            or not 1 <= tau <= 255):
        raise ValueError('invalid packet parameters')
    old = np.asarray(local_rank_transfer(weight, y, tau), dtype=float)
    if 1 <= weight <= 3:
        old[1, 0] = min(old[1, 0], float(rank_return_moment(
            weight, y, tau, phase=phase, average_inputs=True)))
    result = np.zeros((STATES, STATES))
    if weight <= cutoff:
        result[:POOLED, :POOLED] = local_mask_transfer(weight, y, tau, phase)
    else:
        result[0, 0], result[0, POOLED] = old[0]
        result[1:POOLED, 0] = old[1, 0]
        result[1:POOLED, POOLED] = old[1, 1]
    result[POOLED, 0], result[POOLED, POOLED] = old[1]
    if not np.isfinite(result).all() or (result < 0).any():
        raise ArithmeticError('nonfinite or negative packet bound')
    return result


def packet_coefficients(active, y, tau, cutoff=2):
    """Entrywise phase maximum, with NO phase averaging or binomial factor."""
    return np.asarray([np.maximum.reduce([
        packet_matrix(w, y, tau, cutoff, phase) for phase in range(4)])
        for w in range(min(active, 16)+1)])


@lru_cache(maxsize=256)
def hypergeometric_weights(nleft, nright, degree):
    """Exact rational weights; Python integer binomials cannot overflow."""
    return tuple(tuple(Fraction(math.comb(nleft, i)*math.comb(nright, h-i),
                                math.comb(nleft+nright, h))
                       for i in range(max(0, h-nright), min(h, nleft)+1))
                 for h in range(min(degree, nleft+nright)+1))


def normalized_product(left, nleft, right, nright, degree, exact=False):
    """Compose F_h=(sum over h-subsets of a block)/C(block length,h).

    Matrix order is left then right. Truncation loses no coefficient of
    degree at most degree. exact=True supports tiny Fraction replay tests.
    """
    left, right = np.asarray(left), np.asarray(right)
    if (min(nleft, nright, degree) < 0 or left.ndim != 3 or right.ndim != 3
            or left.shape[1] != left.shape[2] or left.shape[1:] != right.shape[1:]
            or len(left) != min(nleft, degree)+1
            or len(right) != min(nright, degree)+1):
        raise ValueError('inconsistent normalized block shapes')
    dtype = object if exact else float
    result = np.zeros((min(degree, nleft+nright)+1, *left.shape[1:]), dtype=dtype)
    weights = hypergeometric_weights(nleft, nright, degree)
    for h, row in enumerate(weights):
        low, high = max(0, h-nright), min(h, nleft)
        indices = np.arange(low, high+1)
        products = np.matmul(left[indices], right[h-indices])
        if exact:
            for weight, value in zip(row, products):
                result[h] += weight*value
        else:
            result[h] = np.einsum('i,ijk->jk', np.asarray(row, dtype=float), products)
    if not exact and (not np.isfinite(result).all() or (result < 0).any()):
        raise ArithmeticError('nonfinite or negative normalized convolution')
    return result


def normalized_power(coefficients, width, packets, degree, exact=False):
    if packets < 0 or width <= 0 or degree < 0:
        raise ValueError('invalid block geometry')
    coefficients = np.asarray(coefficients, dtype=object if exact else float)
    result = np.eye(coefficients.shape[1], dtype=coefficients.dtype)[None, :, :]
    nresult, nfactor = 0, width
    factor = coefficients
    while packets:
        if packets & 1:
            result = normalized_product(result, nresult, factor, nfactor, degree, exact)
            nresult += nfactor
        packets >>= 1
        if packets:
            factor = normalized_product(factor, nfactor, factor, nfactor, degree, exact)
            nfactor *= 2
    return result


def region_coefficients(active, y, tau, rows, cutoff=2):
    if not 0 <= active <= rows or rows <= 0 or rows % 16:
        raise ValueError('invalid active count or regional capacity')
    local = packet_coefficients(active, y, tau, cutoff)
    return normalized_power(local, 16, rows//16, active)


def scaled_region(active, y, tau, rows, cutoff=2):
    coefficients = region_coefficients(active, y, tau, rows, cutoff)
    logs = np.asarray([log_choose(active, w)+w*math.log(tau) for w in range(active+1)])
    scale = float(logs.max())
    result = np.einsum('i,ijk->jk', np.exp(logs-scale), coefficients)
    maximum = float(result.max())
    if not maximum > 0 or not math.isfinite(maximum):
        raise ArithmeticError('invalid regional scale')
    return result/maximum, scale+math.log(maximum)


def bound(active, y, tau, k_log=20, cutoff=2):
    """Return log2 of the relaxed exact-A-row failure bound (FLOATING)."""
    if active <= 0:
        raise ValueError('positive active-row count required')
    g = geometry(k_log, 32, 22)
    region, region_scale = scaled_region(active, y, tau, g['rows'], cutoff)
    vector = np.zeros(STATES)
    vector[0] = 1
    accumulated = 0.
    for _ in range(32):
        vector = vector @ region
        maximum = float(vector.max())
        if not maximum > 0 or not math.isfinite(maximum):
            raise ArithmeticError('terminal contraction underflow or overflow')
        vector /= maximum
        accumulated += region_scale+math.log(maximum)
    terminal = accumulated+math.log(float(vector.sum()))
    return (log_choose(g['rows'], active)-(10*active+1)*math.log(tau)
            -g['threshold']*math.log(y)+terminal)/math.log(2)


def witness(active, point, k_log=20, cutoff=2):
    logtilt, logtau = map(float, point)
    g = geometry(k_log, 32, 22)
    logy = -math.exp(logtilt)*active/g['n']
    if -g['rows']*logy > 600:
        raise ValueError('regional attenuation exceeds this floating screen guard')
    y, tau = math.exp(logy), math.exp(logtau)
    # Roundoff at the closed upper endpoint must not turn 255 into >255.
    tau = min(tau, 255.) if abs(tau-255.) < 1e-10 else tau
    value = bound(active, y, tau, k_log, cutoff)
    return dict(active=active, k_log=k_log, cutoff=cutoff, margin_bits=-value,
                logtilt=logtilt, logtau=logtau, y=y, tau=tau,
                phase_rule='entrywise maximum', rigorous=False,
                method='17-state masks plus irreversible pooled nonzero fallback')


def bounded_search(active, start, k_log=20, cutoff=2, maxiter=30):
    g = geometry(k_log, 32, 22)
    upper = min(13., math.log(500*g['n']/(active*g['rows'])))
    def objective(point):
        return -witness(active, point, k_log, cutoff)['margin_bits']
    result = minimize(objective, start, method='L-BFGS-B',
                      bounds=[(-2., upper), (0., math.log(255.))],
                      options={'maxiter': maxiter, 'maxls': 12, 'ftol': 1e-10})
    value = witness(active, result.x, k_log, cutoff)
    value.update(optimized=True, iterations=int(result.nit), evaluations=int(result.nfev),
                 optimizer_success=bool(result.success), optimizer_message=str(result.message))
    return value


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active', type=int, nargs='+', default=[4, 8, 16])
    parser.add_argument('--cutoff', type=int, choices=(1, 2), default=2)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--optimize', action='store_true')
    parser.add_argument('--maxiter', type=int, default=30)
    args = parser.parse_args()
    for active in args.active:
        start = WITNESSES.get(active, (5., 2.5))
        print(json.dumps(witness(active, start, args.k_log, args.cutoff)), flush=True)
        if args.optimize:
            print(json.dumps(bounded_search(active, start, args.k_log, args.cutoff, args.maxiter)), flush=True)

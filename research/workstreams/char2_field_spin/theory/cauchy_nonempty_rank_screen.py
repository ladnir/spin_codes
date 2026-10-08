"""Exploratory rank cover retaining the MDS minimum nonzero row weight.

The Poisson-binomial row envelope is homogenized only when tau*lambda>=2.
This makes its count function discretely convex, so pairwise balancing of
regional densities increases its expectation. Finite costs are omitted:
reported negative exponents are screens, not finite certificates.
"""

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from regional_refresh_screen import (log_period_radius, log_packet, support_types,
                                    search as previous_search)
from diagonal_refresh_rank import root_factors
from sparse_multirow_zero_set import transfer as zero_set_transfer
from cauchy_joint_return import rank_return_moment


def log_row(alpha, log_tau, log_r, length=32, parity=10):
    if not 0 < alpha < 1 or log_r < math.log(2)-1e-12:
        raise ValueError('need alpha in (0,1) and tau*lambda>=2')
    return float(logsumexp([length*math.log1p(-alpha)]+[
        math.log(math.comb(length, w))-parity*log_tau
        +w*(log_r+math.log(alpha))+(length-w)*math.log1p(-alpha)
        for w in range(parity+1, length+1)]))


def joint_return_radius(point):
    """Candidate return-first rank replacement for the w=1 coefficient.

    On a return the input is nonzero and all four incoming state coordinates
    must be nonzero. Condition on its input diagonal, then use each of the
    four state diagonals once. The returned output has exact weight16.
    This is an alternate exposure, not an extra factor on the old transfer.
    """
    log_x, log_y, log_tau = point[:3]
    table = log_packet(log_x, log_y, log_tau, True, refresh_law='diagonal')
    coeff, z, only_i, both, w, v = support_types()
    factors = root_factors(math.exp(log_tau))
    terms = (coeff+w*log_x+v*log_y-(only_i+np.minimum(both, 4))*log_tau
             +np.log(factors)[z.astype(int)])
    valid = (w > 1) & (v > 0) & (both > 0) & ((v != 1) | (w == 16))
    table[1, 0] = logsumexp(list(terms[valid])+[
        math.log(16)+log_x+16*log_y-4*log_tau])
    shift = float(table.max())
    a, b, c, d = np.exp(table-shift).ravel()
    return shift+math.log((a+d+math.sqrt((a-d)**2+4*b*c))/2)


def zero_set_radius(point,joint_returns=False):
    x,y,tau=map(math.exp,point[:3])
    table=np.asarray(zero_set_transfer(x,y,tau),dtype=float)
    if joint_returns:
        from sparse_multirow_zero_set import local_rank_transfer
        # Rebuild this entry positively; subtracting a dominant old term
        # could lose the much smaller new return coefficient numerically.
        table[1,0]=sum(math.comb(16,w)*x**w*(
            min(local_rank_transfer(w,y,tau)[1][0],
                rank_return_moment(w,y,tau,average_inputs=True)) if w<=3 else
            local_rank_transfer(w,y,tau)[1][0]) for w in range(1,17))
    scale=float(table.max())
    if not math.isfinite(scale) or scale<=0:
        return math.inf
    a,b,c,d=(table/scale).ravel()
    return math.log(scale)+math.log((a+d+math.sqrt((a-d)**2+4*b*c))/2)


def objective(point, alpha, length=32, parity=10, joint_return=False,zero_set=False):
    log_x, log_y, log_tau, log_r = point
    radius = zero_set_radius(point,joint_return) if zero_set else (
        joint_return_radius(point) if joint_return else log_period_radius(point, 1, 'diagonal'))
    return (radius/16
            -alpha*(log_x+log_r-log_tau)-.1*log_y
            +log_row(alpha, log_tau, log_r, length, parity)/length)/math.log(2)


def search(alpha, joint_return=False,zero_set=False):
    prior = previous_search(alpha, 1, refresh_law='diagonal')
    x, y, t, lam = prior['log_parameters']
    starts = [[x, y, t, max(math.log(2), t+lam)],
              [-5, -.02, math.log(255), math.log(255)],
              [-1, -1, math.log(255), math.log(255)]]
    bounds = [(-35, 15), (-25, 0), (0, math.log(255)), (math.log(2), 25)]
    candidates = [minimize(objective, start, args=(alpha,32,10,joint_return,zero_set), method='L-BFGS-B',
                           bounds=bounds, options={'maxiter':300}) for start in starts]
    best = min(candidates, key=lambda result: result.fun)
    return dict(alpha=alpha, prior_exponent=prior['log2_upper_per_N'],
                nonempty_exponent=float(best.fun), parameters=best.x.tolist(),
                rigorous=False, finite_costs_included=False, joint_return=joint_return,
                zero_set=zero_set)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alpha', type=float, nargs='+',
                        default=[.0001, .001, .01, .04, .08, .10])
    parser.add_argument('--joint-return', action='store_true')
    parser.add_argument('--zero-set',action='store_true')
    args = parser.parse_args()
    for alpha in args.alpha:
        print(json.dumps(search(alpha, args.joint_return,args.zero_set)), flush=True)

"""Binary64 screen of two-group distance using positive generating functions.

This includes same-epoch cancellation and all message ranks, but uses a coarse
support grid and floating point. It is not an outward certificate.
"""
import argparse
from math import comb, log

import numpy as np
from scipy.special import gammaln, logsumexp, xlogy
from scipy.optimize import minimize

from bch_joint_support import authenticated_caps, support_caps
from shortened_bound import dimension_caps
from group_rank_one_verify import TILTS
from two_column_moment import census
from two_group_moment import collision_census, worst_regions, self_test


def log_power_moment(matrix, length=128, terminal=None):
    """Scaled positive matrix power; returns log(e_zero M^length 1)."""
    scale = 0.
    value = matrix.copy()
    exponent = length
    row = np.zeros(len(matrix))
    row[0] = 1.
    total = 0.
    while exponent:
        if exponent & 1:
            row = row@value
            norm = row.max()
            if norm == 0:
                return -np.inf
            row /= norm
            total += scale+log(norm)
        exponent >>= 1
        if exponent:
            value = value@value
            norm = value.max()
            value /= norm
            scale = 2*scale+log(norm)
    return total+log(row.sum() if terminal is None else row@terminal)


def log_binomial_mass(n, u, p):
    return gammaln(n+1)-gammaln(u+1)-gammaln(n-u+1)+xlogy(u,p)+xlogy(n-u,1-p)


def log_bound(region, u, v, tilt, p=None, q=None):
    # Bernoulli support is only a coefficient-bound device. Conditioning on
    # both support sizes recovers independent uniform fixed-size subsets.
    p = u/256 if p is None else p
    q = v/256 if q is None else q
    width=max(a for a,b in region)
    pa = [comb(width,a)*p**a*(1-p)**(width-a) for a in range(width+1)]
    pb = [comb(width,b)*q**b*(1-q)**(width-b) for b in range(width+1)]
    matrix = sum(pa[a]*pb[b]*region[a,b] for a in range(width+1) for b in range(width+1))
    return log_power_moment(matrix,256//width)+tilt*209715-log_binomial_mass(256,u,p)-log_binomial_mass(256,v,q)


def optimized_bound(region, u, v, tilt, return_witness=False):
    # Optimize the positive coefficient witness, not the actual distribution.
    # Any p,q in (0,1) is valid; numerical optimization affects tightness only.
    free = [i for i,w in enumerate((u,v)) if w < 256]
    if not free:
        value = log_bound(region,u,v,tilt,1.,1.)
        return (value,(1.,1.)) if return_witness else value
    def probabilities(x):
        result = [1.,1.]
        for i,z in zip(free,x):
            result[i] = 1/(1+np.exp(-z))
        return result
    def evaluate(x):
        return log_bound(region,u,v,tilt,*probabilities(x))
    start = [log(w/(256-w)) for w in (u,v) if w < 256]
    optimum = minimize(evaluate,start,method='L-BFGS-B',bounds=[(-12.,18.)]*len(free),
                       options={'ftol':1e-11,'maxiter':50})
    point = start if evaluate(start) < optimum.fun else optimum.x
    return (evaluate(point),probabilities(point)) if return_witness else evaluate(point)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step', type=int, default=8)
    parser.add_argument('--optimize', action='store_true')
    parser.add_argument('--memory', action='store_true')
    parser.add_argument('--columns',type=int,choices=(1,2),default=2)
    args = parser.parse_args()
    assert args.step > 0
    self_test()
    spectrum = authenticated_caps()
    caps = support_caps(spectrum,g=4,dimensions=dimension_caps())
    # For a coarse diagnostic, use CDF caps themselves as upper counts of
    # each disjoint support bucket. Never use differences as shell counts.
    grid = sorted(set([38,57,67,72,256]+list(range(40,257,args.step))))
    assert not args.memory or args.columns==2
    single = census(args.columns)
    pair = collision_census(args.columns)
    best = np.zeros((len(grid),len(grid)))
    for text in TILTS:
        tilt = float(text)
        if args.memory:
            from two_group_memory import raw_regions
            from two_group_weight import weighted_regions
            region = weighted_regions(raw_regions(single,pair,tilt),1.,1.)
        else:
            region = worst_regions(single,pair,tilt)
        for i,u in enumerate(grid):
            for j in range(i,len(grid)):
                v = grid[j]
                value = (optimized_bound(region,u,v,tilt) if args.optimize else log_bound(region,u,v,tilt))
                best[i,j] = best[j,i] = min(best[i,j],value)
        print('BINARY64 two-group tilt',text,'sample log2 probabilities',
              [(grid[i],round(best[i,i]/log(2),2)) for i in range(0,len(grid),4)],flush=True)
    # This grid is exploratory, not an upper bound between its sampled points.
    total_caps = [sum(row[u] for row in caps) for u in grid]
    terms = [(best[i,j]+log(a)+log(b),grid[i],grid[j])
             for i,a in enumerate(total_caps) for j,b in enumerate(total_caps) if a and b]
    value = log(comb(2048,2))+logsumexp([x[0] for x in terms])
    print('BINARY64 coarse CDF/grid union diagnostic log2:',value/log(2),flush=True)
    print('Dominant grid pairs:',[(u,v,round((z+log(comb(2048,2)))/log(2),3))
                                  for z,u,v in sorted(terms,reverse=True)[:10]],flush=True)
    for h in range(4):
        for k in range(h,4):
            terms = [best[i,j]+log(caps[h][u])+log(caps[k][v])
                     for i,u in enumerate(grid) for j,v in enumerate(grid) if caps[h][u] and caps[k][v]]
            print('Rank pair',h+1,k+1,'grid diagnostic log2',
                  (log(comb(2048,2)*(1 if h==k else 2))+logsumexp(terms))/log(2),flush=True)
    print('No certificate: support grid, coefficient relaxation, and binary64 screen only.')


if __name__ == '__main__':
    main()

"""Diagnostic: retain total row weight for high-rank two-group messages.

Selected support/weight points, binary64 arithmetic: no distance certificate.
The auxiliary weight variables change the bound, not the encoder distribution.
"""
import argparse
from math import comb,log

import numpy as np
from scipy.optimize import minimize
from flint import fmpz_poly

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from group_rank_one_verify import TILTS
from two_column_moment import census
from two_group_moment import collision_census,worst_regions
from two_group_screen import log_bound,optimized_bound


def weighted_regions(raw,rho,sigma,statistic='weight'):
    empty,one,pairs,shapes,ia,ib = raw
    assert statistic in ('weight','all-ones')
    weights = np.array([sum(s) if statistic=='weight' else s.count(4) for s in shapes])
    occupancy = np.array([sum(w!=0 for w in s) for s in shapes])
    result = {(0,0):empty}
    width=len(shapes[0])
    for b in range(1,width+1):
        mask = occupancy==b
        result[b,0] = np.max(one[mask]*(rho**weights[mask])[:,None,None],axis=0)
        result[0,b] = np.max(one[mask]*(sigma**weights[mask])[:,None,None],axis=0)
    weighted = pairs*(rho**weights[ia]*sigma**weights[ib])[:,None,None]
    for a in range(1,width+1):
        for b in range(1,width+1):
            result[a,b] = np.max(weighted[(occupancy[ia]==a)&(occupancy[ib]==b)],axis=0)
    return result


def optimize_symmetric(raw,u,weight,tilt,statistic='weight'):
    # A symmetric witness is enough to bound this symmetric point; it need
    # not be the global optimum. Always retain rho=1 as a fallback.
    unweighted = weighted_regions(raw,1.,1.)
    fallback,(p,q) = optimized_bound(unweighted,u,u,tilt,True)
    if u == 256:
        def objective(x):
            rho = np.exp(x[0])
            return log_bound(weighted_regions(raw,rho,rho,statistic),u,u,tilt,1.,1.)-2*weight*x[0]
        start = [0.]
        bounds = [(-4.,4.)]
    else:
        def objective(x):
            p = 1/(1+np.exp(-x[0]))
            rho = np.exp(x[1])
            return log_bound(weighted_regions(raw,rho,rho,statistic),u,u,tilt,p,p)-2*weight*x[1]
        start = [log(p/(1-p)),0.]
        bounds = [(-12.,18.),(-4.,4.)]
    best = (fallback,0.)
    for rho0 in (0.,-.1,.1,-.5,.5):
        point = start.copy()
        point[-1] = rho0
        fit = minimize(objective,point,method='Powell',bounds=bounds,
                       options={'xtol':1e-6,'ftol':1e-9,'maxiter':80})
        if fit.fun < best[0]:
            best = (fit.fun,float(fit.x[-1]))
    return best[0],np.exp(best[1]),fallback


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supports',nargs='+',type=int,default=[160,192,216,240,256])
    parser.add_argument('--weights',nargs='+',type=int,default=[448,512,576])
    args = parser.parse_args()
    spectrum = authenticated_caps()
    caps = support_caps(spectrum,g=4,dimensions=dimension_caps())[3]
    single,pair = census(),collision_census()
    counts = {u:(fmpz_poly([0]+spectrum[1:u+1])**4) for u in args.supports}
    best = {(u,w):(0.,0.,0.) for u in args.supports for w in args.weights}
    for text in TILTS:
        tilt = float(text)
        raw = worst_regions(single,pair,tilt,return_shapes=True)
        for u,w in best:
            value,rho,fallback = optimize_symmetric(raw,u,w,tilt)
            if value < best[u,w][0]:
                best[u,w] = value,tilt,rho
        print('BINARY64 total-row-weight tilt',text,flush=True)
    for (u,w),(value,tilt,rho) in best.items():
        count = min(caps[u],int(counts[u][w]))
        score = (value+2*log(count)+log(comb(2048,2)))/log(2) if count else -np.inf
        print('point u,W',u,w,'log2 probability',round(value/log(2),4),
              'log2 capped point contribution',round(score,4),'tilt,rho',tilt,round(rho,6),flush=True)
    print('Selected points only; no interval coverage or outward replay; no certificate.')


if __name__=='__main__':
    main()

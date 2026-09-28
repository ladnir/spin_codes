"""Floating-point screen for one-column multi-group support bounds.

Every operator includes collisions within an epoch. Selected support points
and numerical optimization are diagnostics, not a certificate.
"""
from itertools import combinations_with_replacement
from math import comb,log
import argparse

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from flint import ctx

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from group_rank_one_verify import TILTS
from occupancy_model import local_data,build,self_test
from two_group_screen import log_power_moment,log_binomial_mass


def matrix_for_probabilities(region,probabilities):
    masses=np.array([1.])
    for p in probabilities:
        masses=np.convolve(masses,[1-p,p])
    return sum(mass*matrix for mass,matrix in zip(masses,region))


def point_bound(region,supports,tilt,probabilities,terminal=None):
    moment=log_power_moment(matrix_for_probabilities(region,probabilities),256,terminal)
    return moment+tilt*209715-sum(log_binomial_mass(256,u,p) for u,p in zip(supports,probabilities))


def optimize(region,supports,tilt,terminal=None):
    free=[i for i,u in enumerate(supports) if u<256]
    if not free:
        probabilities=[1.]*len(supports)
        return point_bound(region,supports,tilt,probabilities,terminal),probabilities
    def probabilities(x):
        ps=[1.]*len(supports)
        for i,z in zip(free,x):
            ps[i]=1/(1+np.exp(-z))
        return ps
    def objective(x):
        return point_bound(region,supports,tilt,probabilities(x),terminal)
    initial=[log(supports[i]/(256-supports[i])) for i in free]
    fit=minimize(objective,initial,method='L-BFGS-B',bounds=[(-12.,18.)]*len(free),
                 options={'maxiter':60,'ftol':1e-11})
    x=initial if objective(initial)<fit.fun else fit.x
    return objective(x),probabilities(x)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,choices=(3,4),default=3)
    parser.add_argument('--grid',action='store_true')
    parser.add_argument('--fine',action='store_true')
    parser.add_argument('--basis-lattice',action='store_true')
    args=parser.parse_args()
    spectrum=authenticated_caps()
    caps=support_caps(spectrum,g=4,dimensions=dimension_caps())
    if args.basis_lattice:
        from basis_lattice import improve_caps,self_test as basis_test
        basis_test()
        caps=improve_caps(caps,spectrum)
    counts=[sum(row[u] for row in caps) for u in range(257)]
    data=local_data(args.groups)
    ctx.prec=192
    self_test(data)
    if args.grid:
        points=list(combinations_with_replacement((38,40,57,67,72,76,80,96,128,176,216,256),args.groups))
    else:
        points=[(u,)*args.groups for u in (38,57,72,76,80,96,100,128,176,216,256)]
        points += [(38,)+(256,)*(args.groups-1),(38,72)+(256,)*(args.groups-2)]
    best={u:0. for u in points}
    witnesses={u:None for u in points}
    tilts=TILTS
    if args.fine:
        from rank_three_flags import TILTS as tilts
    for text in tilts:
        _,region=build(data,text)
        for supports in points:
            value,_=optimize(region,supports,float(text))
            if value<best[supports]:
                best[supports]=value
                witnesses[supports]=text
        print('BINARY64 occupancy',args.groups,'tilt',text,flush=True)
    terms=[]
    for supports,value in best.items():
        score=(value+sum(log(counts[u]) for u in supports)+log(comb(2048,args.groups)))/log(2)
        terms.append((score,supports,value/log(2),witnesses[supports]))
    print('Worst sampled support vectors (log2 contribution, supports, log2 probability, tilt):',flush=True)
    for row in sorted(terms,reverse=True)[:20]:
        print(row,flush=True)
    print('Point/grid diagnostic only; no complete support coverage or outward union.')


if __name__=='__main__':
    main()

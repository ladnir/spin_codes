"""Selected-point, binary64 inner screen for the NEW independent-row route.

Reads the shape-resolved parent pilot cache only after its dependency hash
matches. Does not modify the shared-route driver, cache, or certificates.
No outward replay or full support/occupancy coverage is claimed.
"""
import argparse
from fractions import Fraction as Q
from math import comb, log
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from support import weighted_cdf_upper

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from categorical_screen import code_hash
from categorical_model import normalization_test
from bch_joint_support import authenticated_caps
from occupancy_sensitivity import float_placement
from occupancy_screen import matrix_for_probabilities
from occupancy_memory import TERMINAL
from two_group_screen import log_power_moment, log_binomial_mass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', required=True)
    parser.add_argument('--groups', type=int, default=64)
    parser.add_argument('--penalty', default='.75')
    parser.add_argument('--supports', type=int, nargs='+', default=[80,96,128,144,160,176,192,224,256])
    args = parser.parse_args()
    rho = Q(args.penalty)
    if not 0 < rho <= 1 or not 1 <= args.groups <= 2048 or any(not 1 <= u <= 256 for u in args.supports):
        parser.error('invalid penalty, group count, or supports')
    with np.load(args.cache, allow_pickle=False) as saved:
        if str(saved['code_hash']) != code_hash():
            raise ValueError('parent per-shape operator cache dependency mismatch')
        tilt = float(str(saved['tilt']))
        data = [(saved[f's{j}'],saved[f't{j}'],saved[f'c{j}']) for j in range(33)]
    normalization_test(data)
    # A factor rho for each full packet belongs to the inner moment;
    # reciprocal weighting belongs to the exact averaged outer CDF.
    ops = np.array([(matrices * float(rho)**counts[:,3,None,None]).max(axis=0)
                    for counts,matrices,_ in data])
    assert np.isfinite(ops).all() and (ops >= 0).all()
    print('NEW independent-row binary64 pilot; q',args.groups,'tilt',tilt,'penalty',rho,flush=True)
    region = float_placement(ops,args.groups)
    spectrum = authenticated_caps()
    counts = weighted_cdf_upper(spectrum,1 << 128,full_weight=1/rho)
    print('support,log2_weighted_count_cap,log2_point_upper,p',flush=True)
    for u in args.supports:
        def objective(p):
            moment = log_power_moment(matrix_for_probabilities(region,[p]*args.groups),256,TERMINAL)
            return moment+tilt*209715-args.groups*log_binomial_mass(256,u,p)
        if u == 256:
            p = 1.
        else:
            fit = minimize_scalar(lambda z:objective(1/(1+np.exp(-z))),
                                  bounds=(-12,18),method='bounded')
            p = 1/(1+np.exp(-fit.x))
        logcount = log(counts[u].numerator)-log(counts[u].denominator)
        score = (objective(p)+args.groups*logcount+log(comb(2048,args.groups)))/log(2)
        print(u,logcount/log(2),score,p,sep=',',flush=True)
    print('Nine-coordinate pilot only. Selected homogeneous supports, not a complete cover or certificate.')


if __name__ == '__main__':
    main()

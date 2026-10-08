"""Floating higher-distance screen for the checked width16/state8 map.

Only the deterministic GRS outer and output threshold change. The local
rank cover is the preserved superregular-return bound. Positive results
are class screens, not a whole-code certificate; failed bounds are not
counterexamples. No native implementation or benchmark runs here.
"""
import argparse
from fractions import Fraction
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry
from superregular_return_screen import bound as reference_bound
from state8_density_screen import bound as density_reference_bound


def output_threshold(n, distance):
    distance = Fraction(distance)
    if not 0 < distance < 1 or n < 1:
        raise ValueError('positive length and distance in (0,1) required')
    return (distance.numerator*n+distance.denominator-1)//distance.denominator-1


def bound(active, point, k_log=20, length=256, dimension=128,
          distance=Fraction(3, 10)):
    g = geometry(k_log, length, dimension)
    logy = -math.exp(point[0])*active/g['n']
    adjustment = -(output_threshold(g['n'], distance)-g['threshold'])*logy/math.log(2)
    return reference_bound(active, point, k_log, length, dimension, 8, True, 16)+adjustment


def search(active, k_log=20, length=256, dimension=128,
           distance=Fraction(3, 10), iterations=35):
    if not 0 < dimension < length <= 256:
        raise ValueError('invalid GF256 GRS geometry')
    g = geometry(k_log, length, dimension)
    if not 1 <= active <= g['rows'] or iterations < 1:
        raise ValueError('invalid active count or optimizer budget')
    options = []
    for point in ((5., math.log(255)), (5., math.log(16))):
        value = minimize(lambda p: bound(active, p, k_log, length, dimension, distance),
                         point, method='L-BFGS-B',
                         bounds=[(-2, 10), (0, math.log(255))],
                         options={'maxiter': iterations, 'ftol': 1e-10})
        options.append(value)
    best = min(options, key=lambda v: v.fun if np.isfinite(v.fun) else math.inf)
    return dict(outer=[length, dimension], rate=str(Fraction(dimension, length)),
                k_log=k_log, rows=g['rows'], n=g['n'], active=active,
                distance=str(Fraction(distance)), output_threshold=output_threshold(g['n'], distance),
                margin_bits=-float(best.fun), logtilt=float(best.x[0]), logtau=float(best.x[1]),
                optimizer_success=bool(best.success), iterations=int(best.nit),
                rigorous=False, whole_code=False, scope='one sampled active-row class',
                width=16, state=8, local_cover='preserved coarse superregular return')


def density_bound(point, weight, k_log=20, length=256, dimension=128,
                  distance=Fraction(3, 10)):
    g = geometry(k_log, length, dimension)
    adjustment = -(output_threshold(g['n'], distance)-g['threshold'])*point[1]/math.log(2)
    return density_reference_bound(point, weight, k_log, length, dimension, 16, True)+adjustment


def density_search(alpha, k_log=20, length=256, dimension=128,
                   distance=Fraction(3, 10), iterations=55, start=None):
    if not 0 < dimension < length <= 256 or not 0 < alpha < 1:
        raise ValueError('invalid geometry or input density')
    g = geometry(k_log, length, dimension)
    weight = max(1, round(alpha*g['n']))
    starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
    if start is not None:
        starts.insert(0, start)
    options = [minimize(lambda p: density_bound(p, weight, k_log, length, dimension, distance),
                        p, method='L-BFGS-B',
                        bounds=[(-30, 12), (-25, 0), (0, math.log(255)), (-20, 20)],
                        options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 20})
               for p in starts]
    best = min(options, key=lambda v: v.fun if np.isfinite(v.fun) else math.inf)
    return dict(outer=[length, dimension], rate=str(Fraction(dimension, length)),
                k_log=k_log, input_weight=weight, input_density=weight/g['n'],
                distance=str(Fraction(distance)), margin_bits=-float(best.fun),
                log_parameters=best.x.tolist(), optimizer_success=bool(best.success),
                iterations=int(best.nit), rigorous=False, whole_code=False,
                scope='one sampled inverse-input-weight class', width=16, state=8,
                local_cover='joint state8 selected-output-zero/return rank')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', type=int, nargs=2, default=(256, 128))
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--distance', type=Fraction, default=Fraction(3, 10))
    parser.add_argument('--active', type=int, nargs='+', default=[1, 4, 16, 64])
    parser.add_argument('--iterations', type=int, default=35)
    parser.add_argument('--density', type=float, nargs='+')
    args = parser.parse_args()
    if args.density:
        previous = None
        for alpha in args.density:
            result = density_search(alpha, args.k_log, *args.outer, args.distance,
                                    args.iterations, previous)
            previous = result['log_parameters']
            print(json.dumps(result), flush=True)
    else:
        for a in args.active:
            print(json.dumps(search(a, args.k_log, *args.outer, args.distance, args.iterations)), flush=True)

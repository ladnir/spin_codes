"""Bounded FLOAT proposals for a diverse frontier, never a coverage claim."""
from fractions import Fraction as F
from math import log2

import numpy as np
from scipy.optimize import minimize

from hd_shifted_swap_screen import CONSTRUCTION, bound
from hd_shifted_unoriented_screen import unsigned_histograms


FRONTIER = ((1, 0, 0), (2, 0, 0), (2, 1, 0),
            (16, 0, 0), (16, 4, 4), (16, 8, 7),
            (128, 0, 0), (128, 50, 50), (128, 64, 63),
            (512, 0, 0), (512, 200, 200), (512, 256, 255),
            (992, 450, 420), (1024, 400, 400), (1024, 512, 511),
            (2048, 0, 0), (2048, 800, 800), (2048, 1024, 1023),
            (4096, 0, 0), (4096, 1024, 1024), (4096, 2048, 2047),
            (5119, 1023, 1023), (5119, 2559, 2559),
            (1024, 512, 0), (2048, 1024, 0), (4096, 2048, 0),
            (4096, 0, 2048), (1024, 400, 380))


def exact_witness(result):
    """Propose rational inputs; directed replay must use these actual values."""
    return dict(z=str(F(result['mark'])), alpha=str(F(result['alpha'])),
                markers=[[str(F(2**value).limit_denominator(10**12)) for value in row]
                         for row in result['marker_log_odds']],
                construction=CONSTRUCTION, pair_swaps=True,
                orientation_aggregation='independent-fair-swaps')


def proposal(profile, *, max_evaluations=180, grid=64, start=None):
    if max_evaluations < 8:
        raise ValueError('At least eight objective evaluations required')
    profile = tuple(profile)
    histograms = unsigned_histograms(4096, *profile)
    tied = histograms[0] == histograms[1]
    live = np.array([n > 0 for hist in histograms for n in hist[1:]])
    if tied:
        live[2:] = False
    base = np.array([[log2(max(n, 1)/max(hist[0], 1)) if n else -40.
                      for n in hist[1:]] for hist in histograms]).ravel()
    best, evaluations = None, 0

    def objective(theta):
        nonlocal best, evaluations
        evaluations += 1
        ell, alpha, *markers = theta
        if not 1e-6 <= ell <= 16 or not .001 <= alpha <= 1 or any(abs(x) > 24 for x in markers):
            return 1e15
        odds = base.copy()
        odds[live] = markers
        if tied:
            odds[2:] = odds[:2]
        try:
            with np.errstate(over='raise', invalid='raise'):
                result = bound(profile, 2**-ell, alpha, odds.reshape((2, 2)), method='global')
            score = result['log2_upper_screen']
        except (ArithmeticError, OverflowError, FloatingPointError):
            return 1e15
        if np.isfinite(score) and (best is None or score < best['log2_upper_screen']):
            best = result
        return score if np.isfinite(score) else 1e15

    seeds = []
    for alpha in (.1, .25, .5, 1.):
        marker = base+np.tile((.5, .9), 2)*(alpha/.1)
        seeds.append([max(1e-4, profile[0]/1000), alpha, *marker[live]])
    if start is not None:
        seeds.insert(0, [-log2(start['mark']), start['alpha'],
                         *np.asarray(start['marker_log_odds']).ravel()[live]])
    scores = [objective(seed) for seed in seeds]
    seed = np.asarray(seeds[int(np.argmin(scores))])
    simplex = np.tile(seed, (len(seed)+1, 1))
    for axis, value in enumerate(seed):
        simplex[axis+1, axis] += (.25 if axis >= 2 else
                                max(1e-4, .06*value) if axis == 0 else
                                -.03 if value > .95 else .02)
    minimize(objective, seed, method='Nelder-Mead',
             options=dict(maxfev=max_evaluations-evaluations, initial_simplex=simplex,
                          xatol=1e-6, fatol=1e-3))
    if best is None:
        raise ArithmeticError('No finite swap witness')
    result = dict(profile=list(profile), objective_evaluations=evaluations,
                  evaluation_budget=max_evaluations, positive=best, exact_witness=exact_witness(best),
                  directed=False, whole_code_certificate=False, construction=CONSTRUCTION)
    if grid:
        try:
            with np.errstate(over='raise', invalid='raise'):
                result['alias'] = bound(profile, best['mark'], best['alpha'], best['marker_log_odds'],
                                        method='alias', grid=grid)
        except (ArithmeticError, FloatingPointError) as error:
            result['alias_error'] = str(error)
    return result


if __name__ == '__main__':
    import argparse, json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=int, nargs=3, action='append')
    parser.add_argument('--max-evaluations', type=int, default=180)
    parser.add_argument('--grid', type=int, default=64)
    args = parser.parse_args()
    for profile in args.profile or FRONTIER:
        print(json.dumps(proposal(profile, max_evaluations=args.max_evaluations, grid=args.grid),
                         allow_nan=False), flush=True)

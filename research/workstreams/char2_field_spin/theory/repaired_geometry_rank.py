"""Bounded FLOATING comparison of GRS33/22 and GRS64/44 rank screens.

Uses existing transfer bounds unchanged. The return-deletion comparison is
explicitly invalid as a proof and serves only to diagnose the bottleneck.
No output from this script is a directed or whole-code certificate.
"""

import argparse
from contextlib import contextmanager
import json
import math

from scipy.optimize import minimize

import cauchy_regional_rank_screen as regional
from cauchy_firstmoment_screen import geometry
from regional_refresh_screen import objective as density_objective


STARTS = {4: (5.32479700865, 4.96215263668),
          8: (5.13883585678, 2.88863774319),
          16: (5.07867455120, 2.33096175675),
          32: (5.019956, 2.11492535),
          64: (4.989729476, 2.030286964),
          128: (4.989729476, 2.030286964)}


@contextmanager
def invalid_delete_returns():
    original = regional.local_rank_transfer
    def without_returns(*args, **kwargs):
        value = original(*args, **kwargs)
        value[1][0] = 0
        return value
    regional.local_rank_transfer = without_returns
    try:
        yield
    finally:
        regional.local_rank_transfer = original


def sparse_screen(active, length, dimension, maxiter=35):
    g = geometry(20, length, dimension)
    def objective(point):
        logy = -math.exp(float(point[0]))*active/g['n']
        tau = min(255., math.exp(float(point[1])))
        return regional.bound(active, math.exp(logy), tau, 20, True, length, dimension)
    start = STARTS[active]
    initial = objective(start)
    result = minimize(objective, start, method='L-BFGS-B',
                      bounds=[(-2, 13), (0, math.log(255))],
                      options={'maxiter': maxiter, 'maxls': 12, 'ftol': 1e-10})
    best = result.x if result.fun < initial else start
    final = objective(best)
    with invalid_delete_returns():
        deleted = objective(best)
    return dict(outer=[length, dimension], k_log=20, active=active, rows=g['rows'],
                n=g['n'], threshold=g['threshold'], initial_margin_bits=-initial,
                margin_bits=-final, logtilt=float(best[0]), logtau=float(best[1]),
                invalid_return_deleted_margin_bits=-deleted,
                iterations=int(result.nit), evaluations=int(result.nfev),
                optimizer_success=bool(result.success), optimizer_message=str(result.message),
                rigorous=False, whole_code=False, starts=1,
                method='unchanged exact-regional 2-state joint-return rank bound')


def density_screen(alpha, length, dimension, start=None, maxiter=80):
    # Existing diagonal-refresh support cover; its spectral objective omits
    # finite route-conditioning and boundary factors. No interval assertion.
    if start is None:
        start = [math.log(32627/9127), math.log(44/291),
                 math.log(54371/3437), math.log(2722/4197)]
    result = minimize(density_objective, start,
                      args=(alpha, 1, length, length-dimension, 'diagonal'),
                      method='L-BFGS-B',
                      bounds=[(-30, 12), (-25, 0), (0, math.log(255)), (-20, 20)],
                      options={'maxiter': maxiter, 'maxls': 12, 'ftol': 1e-11})
    return dict(outer=[length, dimension], alpha=alpha,
                log2_upper_per_N=float(result.fun), log_parameters=result.x.tolist(),
                iterations=int(result.nit), optimizer_success=bool(result.success),
                rigorous=False, whole_code=False, finite_factors_included=False)


def checks():
    """Check geometry and diagnostic isolation; not a new transfer proof."""
    original = regional.local_rank_transfer
    cases = 0
    for length, dimension, expected in ((33, 22, (47680, 1573440, 157343)),
                                       (64, 44, (23872, 1527808, 152780))):
        g = geometry(20, length, dimension)
        assert (g['rows'], g['n'], g['threshold']) == expected
        logy = -math.exp(STARTS[4][0])*4/g['n']
        args = (4, math.exp(logy), math.exp(STARTS[4][1]), 20, True, length, dimension)
        before = regional.bound(*args)
        with invalid_delete_returns():
            deleted = regional.bound(*args)
        assert regional.local_rank_transfer is original
        after = regional.bound(*args)
        assert before == after and math.isfinite(before) and deleted <= before
        cases += 3
    try:
        with invalid_delete_returns():
            raise RuntimeError('test restoration on exceptional exit')
    except RuntimeError:
        pass
    assert regional.local_rank_transfer is original
    return cases+1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active', type=int, nargs='+', default=[4, 8, 16, 32, 64, 128])
    parser.add_argument('--geometries', type=int, nargs='+', choices=(33, 64), default=[33, 64])
    parser.add_argument('--maxiter', type=int, default=35)
    parser.add_argument('--dense', action='store_true')
    parser.add_argument('--dense-only', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.check:
        print(json.dumps({'geometry_and_diagnostic_checks': checks(), 'certificate': False}))
        raise SystemExit(0)
    for length in args.geometries:
        dimension = 22 if length == 33 else 44
        if not args.dense_only:
            for active in args.active:
                print(json.dumps(sparse_screen(active, length, dimension, args.maxiter)), flush=True)
        if args.dense or args.dense_only:
            start = None
            for alpha in (.12, .11, .10, .08):
                value = density_screen(alpha, length, dimension, start)
                start = value['log_parameters']
                print(json.dumps(value), flush=True)

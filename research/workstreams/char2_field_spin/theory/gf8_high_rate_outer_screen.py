"""Float-only outer-length exploration with the repaired GF256 8/4 inner.

Fix one regional route, independent nonzero input/state diagonals, and the
audited cubic maps. The outer is an [b,k] MDS code over GF256, b<=256.
Changing b changes the number of regions, not the packet geometry.
Density and active-row outputs are individual bound evaluations, not full
certificates or lower bounds on the construction's failure probability.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from gf8_packet_sparse_screen import log_coefficients
from gf8_packet_width_screen import packet
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


RETURN_MINIMA = (8, 6, 5, 4, 3, 2, 2, 1)


def profile(length, dimension, k_log=20):
    if not 0 < dimension < length <= 256 or not 6 <= k_log <= 30:
        raise ValueError('require 0<dimension<length<=256 and 6<=k_log<=30')
    result = geometry(k_log, length, dimension)
    result.update(outer=[length, dimension], q=256, width=8, state=4,
                  requested_k=2**k_log, padded_k=result['k'], k_log=k_log,
                  requested_rate=str(F(2**k_log, result['n'])),
                  nominal_rate=str(F(dimension, length)), distance='1/10',
                  return_minima=list(RETURN_MINIMA), rigorous=False, whole_code=False)
    return result


def bulk_source_log(alpha, logtau, logtilt, length, dimension):
    """Optimistic asymptotic diagnostic, not a valid finite source upper bound.

    Replaces log(exp(empty)+exp(occupied)) by max(empty,occupied).
    The omitted sum costs at most ln(2) per outer row. This diagnoses the
    relaxed two-product source envelope; it does not replace its proof.
    """
    empty = length*math.log1p(-alpha)
    occupied = -(length-dimension)*logtau + length*np.logaddexp(
        math.log1p(-alpha), math.log(alpha)+logtau+logtilt)
    return float(max(empty, occupied))


def middle_bound(point, weight, length=256, dimension=192, k_log=20,
                 components=False, source='feasible'):
    g = profile(length, dimension, k_log)
    if not 0 < weight < g['n'] or source not in ('feasible', 'relaxed', 'bulk'):
        raise ValueError('invalid weight or source diagnostic')
    logx, logy, logtau, logtilt = map(float, point)
    alpha = weight/g['n']
    row = (bulk_source_log(alpha, logtau, logtilt, length, dimension) if source == 'bulk'
           else source_row_log(alpha, logtau, logtilt, length, dimension, source == 'feasible'))
    logs = packet(logx, logy, logtau, 8, 4, RETURN_MINIMA)
    scale = float(logs.max())
    steps = g['n']//8
    terms = dict(route=length*modal_log_penalty(g['rows']), source=g['rows']*row,
                 inner=steps*scale+log_terminal_power(np.exp(logs-scale), steps),
                 selection=-weight*(logx+logtilt), output=-g['threshold']*logy,
                 rank_threshold=-logtau)
    value = sum(terms.values())/math.log(2)
    return (value, {key: x/math.log(2) for key, x in terms.items()}) if components else value


def middle_search(density, length=256, dimension=192, k_log=20,
                  iterations=60, start=None, source='feasible'):
    if not 0 < density < 1:
        raise ValueError('density must be strictly between zero and one')
    g = profile(length, dimension, k_log)
    weight = max(1, round(density*g['n']))
    starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(lambda p: middle_bound(p, weight, length, dimension, k_log, source=source),
                        point, method='L-BFGS-B',
                        bounds=[(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)],
                        options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25})
               for point in starts]
    best = min(choices, key=lambda r: r.fun if np.isfinite(r.fun) else math.inf)
    value, terms = middle_bound(best.x, weight, length, dimension, k_log, True, source)
    g.update(input_weight=weight, input_density=weight/g['n'], margin_bits=-value,
             log_parameters=best.x.tolist(), log2_components=terms,
             optimizer_success=bool(best.success), iterations=int(best.nit), source=source,
             scope=('optimistic source-envelope diagnostic, not an upper bound' if source == 'bulk'
                    else 'one sampled inverse-input-weight class, floating arithmetic'))
    return g


def sparse_bound(active, point, length=256, dimension=192, k_log=20, components=False):
    g = profile(length, dimension, k_log)
    if not 1 <= active <= g['rows']:
        raise ValueError('invalid active-row class')
    logtilt, logtau = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    tau = min(255., math.exp(logtau))
    coefficients = log_coefficients(active, math.exp(logy), tau, g['rows'], 8, 4, RETURN_MINIMA)
    regional = regional_from_coefficients(coefficients, active, g['rows'], math.log(tau))
    terms = dict(source=log_choose(g['rows'], active)-((length-dimension)*active+1)*math.log(tau),
                 output=-g['threshold']*logy, regional=log_terminal(regional, length))
    value = sum(terms.values())/math.log(2)
    return (value, {key: x/math.log(2) for key, x in terms.items()}) if components else value


def sparse_search(active, length=256, dimension=192, k_log=20, iterations=40, start=None):
    starts = [(math.log(640), math.log(255)), (math.log(640), math.log(8))]
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(lambda p: sparse_bound(active, p, length, dimension, k_log),
                        point, method='L-BFGS-B', bounds=[(-2, 10), (0, math.log(255))],
                        options={'maxiter': iterations, 'ftol': 1e-10, 'maxls': 16})
               for point in starts]
    best = min(choices, key=lambda r: r.fun if np.isfinite(r.fun) else math.inf)
    value, terms = sparse_bound(active, best.x, length, dimension, k_log, True)
    g = profile(length, dimension, k_log)
    g.update(active=active, margin_bits=-value, log_parameters=best.x.tolist(), log2_components=terms,
             optimizer_success=bool(best.success), iterations=int(best.nit),
             scope='one active-row class, floating arithmetic')
    return g


def expected_window_dimension(length, dimension, rows, occupancy=F(1, 10)):
    """Exact E[D] for one fixed equal-width window in every region.

    Uniform independent regional permutations give Binomial(b,occupancy)
    hits for each source row. Linearity of expectation does not require
    independence between source rows. This is not an adaptive-window test.
    """
    if not 0 < dimension < length or rows <= 0 or not 0 <= occupancy <= 1:
        raise ValueError('invalid fixed-window geometry')
    occupancy = F(occupancy)
    parity = length-dimension
    return rows*sum((h-parity)*math.comb(length, h)*occupancy**h*(1-occupancy)**(length-h)
                    for h in range(parity+1, length+1))


def window_screen(length, dimension, k_log=20):
    g = profile(length, dimension, k_log)
    width = (g['threshold']//length//8)*8
    expected = expected_window_dimension(length, dimension, g['rows'], F(width, g['rows']))
    threshold = 4*length+1
    g.update(window_width=width, total_window_size=width*length,
             expected_dimension=float(expected), actual_obstruction_threshold=threshold,
             fixed_window_markov_upper=float(min(F(1), expected/threshold)),
             scope='one fixed aligned window per region; no adaptive search or code certificate')
    return g


def window_witness(length, dimension, k_log=20, seed=1):
    """Test a realized route, without interpreting seeds as probability.

    The D>4b criterion proves failure for this realized route, independently
    of the nonzero diagonal values. Each row shortens to max(0,h-parity)
    dimensions; the b window-exit states impose at most 4b equations.
    """
    result = window_screen(length, dimension, k_log)
    rng = np.random.default_rng(seed)
    hits = np.zeros(result['rows'], dtype=np.int16)
    for _ in range(length):
        hits[rng.permutation(result['rows'])[:result['window_width']]] += 1
    dimension_sum = int(np.maximum(hits-(length-dimension), 0).sum())
    result.update(seed=seed, realized_shortened_dimension=dimension_sum,
                  exit_state_constraints=4*length,
                  kernel_dimension_lower=max(0, dimension_sum-4*length),
                  actual_route_obstruction=dimension_sum > 4*length,
                  scope='one realized seeded route; deterministic dimension obstruction only')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', type=int, nargs=2, action='append')
    parser.add_argument('--mode', choices=('middle', 'sparse', 'window'), default='middle')
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--density', type=float, nargs='+', default=[.03, .06, .10, .14, .20, .24])
    parser.add_argument('--active', type=int, nargs='+', default=[1, 8, 32, 64, 128])
    parser.add_argument('--source', choices=('feasible', 'relaxed', 'bulk'), default='feasible')
    parser.add_argument('--iterations', type=int, default=60)
    parser.add_argument('--window-seed', type=int,
                        help='also test one realized route; not a failure-probability estimate')
    args = parser.parse_args()
    for length, dimension in args.outer or [(256, 176), (256, 184), (256, 192), (256, 200)]:
        previous = None
        if args.mode == 'window':
            result = (window_screen(length, dimension, args.k_log) if args.window_seed is None else
                      window_witness(length, dimension, args.k_log, args.window_seed))
            print(json.dumps(result), flush=True)
            continue
        for value in args.density if args.mode == 'middle' else args.active:
            if args.mode == 'middle':
                result = middle_search(value, length, dimension, args.k_log,
                                       args.iterations, previous, args.source)
            else:
                result = sparse_search(value, length, dimension, args.k_log, args.iterations, previous)
            previous = result['log_parameters']
            print(json.dumps(result), flush=True)

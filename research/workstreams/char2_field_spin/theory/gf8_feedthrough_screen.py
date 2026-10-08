"""Positive quotient-rank cover for cubic8/4 with nonidentity feedthrough.

Y=(I+mu AC)X+eta AR, Q'=R+CX, with fresh independent nonzero input
and state diagonals. Selected output equations are eliminated against the
actual feedthrough columns before any fresh input diagonal is charged.
The remaining state-only matrix is bounded by its exact column matroid.
Numerical screens are not full distance certificates.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import log_choose, log_terminal_power
from gf8_high_rate_outer_screen import profile
from gf8_packet_maps import MUL, INV, cubic_matrices, matmul, rank
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from live_kernel_rank import live_kernel_bound
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


BIRTH_MINIMA = (8, 6, 5, 5, 3, 2, 2, 1)
RETURN_MINIMA = (8, 6, 6, 4, 4, 2, 2, 1)
SUBSETS = tuple(tuple(i for i in range(8) if mask >> i & 1) for mask in range(256))
STATE_SUBSETS = tuple(tuple(i for i in range(4) if mask >> i & 1) for mask in range(16))


def eliminate_inputs(inputs, states):
    """Return input rank and state-only equations after exact row operations."""
    if not inputs:
        return 0, ()
    columns = len(inputs[0])
    work = [list(left)+list(right) for left, right in zip(inputs, states)]
    pivot_row = 0
    for column in range(columns):
        pivot = next((r for r in range(pivot_row, len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = INV[work[pivot_row][column]]
        work[pivot_row] = [MUL[inverse][v] for v in work[pivot_row]]
        for r in range(pivot_row+1, len(work)):
            coefficient = work[r][column]
            if coefficient:
                work[r] = [a ^ MUL[coefficient][b] for a, b in zip(work[r], work[pivot_row])]
        pivot_row += 1
    return pivot_row, tuple(tuple(row[columns:]) for row in work[pivot_row:])


def state_profile(matrix):
    """(rank G, spark(G)-1), where full column rank has no live kernel.

    Every nonzero kernel vector uses a dependent column set. A minimal
    such set has size spark, so independent nonzero diagonal pivots bound
    each fixed-vector zero probability by 255^{-(spark-1)}.
    """
    r = rank(matrix)
    if r == 4:
        return 4, 4
    for size in range(1, r+2):
        for support in STATE_SUBSETS[1:]:
            if len(support) == size and rank([[row[j] for j in support] for row in matrix]) < size:
                return r, size-1
    raise AssertionError('a rank-deficient four-column matrix has a circuit')


def full_support_kernel(matrix, columns):
    """Exact feasibility of a full-support kernel vector when q>columns."""
    r = rank(matrix)
    if r == columns:
        return False
    return all(rank([[v for j, v in enumerate(row) if j != omitted] for row in matrix]) == r
               for omitted in range(columns))


def translate_mask(mask, translation):
    return sum(1 << (i ^ translation) for i in SUBSETS[mask])


@lru_cache(maxsize=1)
def rank_profiles():
    """Exhaust all support/output masks; use the proved GF4 translation symmetry.

    State-only circuit bounds take the worst of all eight map phases.
    Coefficients are sums over raw supports, not per-weight pointwise bounds.
    """
    phases = [cubic_matrices(p) for p in range(8)]
    ac = matmul(*phases[0])
    if any(matmul(*pair) != ac for pair in phases):
        raise AssertionError('feedthrough must be phase independent')
    feedthrough = tuple(tuple(int(i == j) ^ MUL[4][ac[i][j]] for j in range(8)) for i in range(8))
    continuation = [Counter() for _ in range(9)]
    birth = [Counter() for _ in range(9)]
    seen = set()
    for mask in range(256):
        if mask in seen:
            continue
        orbit = {translate_mask(mask, a) for a in range(4)}
        seen.update(orbit)
        inputs = SUBSETS[mask]
        weight = len(inputs)
        for zeros in SUBSETS:
            selected = tuple(tuple(feedthrough[z][i] for i in inputs) for z in zeros)
            state_columns = tuple(tuple(MUL[2][a[z][j]] for a, _ in phases for j in range(4))
                                  for z in zeros)
            input_rank, residual = eliminate_inputs(selected, state_columns)
            state_cases = [state_profile(tuple(tuple(row[4*p:4*p+4]) for row in residual)) for p in range(8)]
            state_rank = state_cases[0][0]
            if any(r != state_rank for r, _ in state_cases):
                raise AssertionError('phase change should preserve state-only rank')
            circuit_charge = min(charge for _, charge in state_cases)
            continuation[weight][len(zeros), input_rank, state_rank, circuit_charge] += len(orbit)
            if weight and full_support_kernel(selected, weight):
                # If selected output equations force C X=0, this is not a birth.
                c = phases[0][1]
                joined = selected + tuple(tuple(row[i] for i in inputs) for row in c)
                if rank(joined) > input_rank:
                    birth[weight][len(zeros), input_rank] += len(orbit)
    if seen != set(range(256)):
        raise AssertionError('missing raw supports')
    return tuple(continuation), tuple(birth)


@lru_cache(maxsize=256, typed=True)
def root_potential(tau, state_rank, circuit_charge, q=256):
    """Uniform live-state moment bound for the state-only matrix profile."""
    if state_rank == 4:
        return tau*0
    if state_rank == 0:
        return tau**0
    exact = not isinstance(tau, float)
    if not 1 <= tau <= q-1:
        raise ValueError('rank tilt must not exceed the nonzero scalar alphabet')
    probability = F(1, (q-1)**circuit_charge) if exact else (q-1.)**(-circuit_charge)
    answers = []
    for incoming in range(1, 5):
        maximum = min(incoming, 4-state_rank)
        mean = min((q**incoming-1)*probability, q**maximum-1)
        if exact:
            mean = F(mean)
        last_x, last_y = 0, tau*0
        for h in range(1, maximum+1):
            x, y = q**h-1, tau**(h-incoming)
            if mean <= x:
                answers.append(last_y+(y-last_y)*(mean-last_x)/(x-last_x))
                break
            last_x, last_y = x, y
    return max(answers)


@lru_cache(maxsize=256, typed=True)
def root_cap_weights(tau, cap, state_rank, circuit_charge, q=256):
    """Simultaneous upper-tail measure for surviving state-form rank."""
    tails = []
    exact = not isinstance(tau, float)
    for threshold in range(1, 5):
        candidates = []
        for incoming in range(1, cap+1):
            maximum = min(incoming, 4-state_rank)
            if threshold > maximum:
                candidates.append(tau*0)
                continue
            probability = F(1, (q-1)**circuit_charge) if exact else (q-1.)**(-circuit_charge)
            mean = min((q**incoming-1)*probability, q**maximum-1)
            if exact:
                mean = F(mean)
            xs = [q**h-1 for h in range(maximum+1)]
            ys = [tau**(h-incoming) if h >= threshold else tau*0 for h in range(maximum+1)]
            best = max(y for x, y in zip(xs, ys) if x <= mean)
            for left in range(maximum+1):
                if xs[left] > mean:
                    break
                for right in range(left+1, maximum+1):
                    if xs[right] >= mean:
                        ratio = (mean-xs[left])/(xs[right]-xs[left])
                        best = max(best, ys[left]+ratio*(ys[right]-ys[left]))
            candidates.append(best)
        tails.append(max(candidates))
    for index in range(2, -1, -1):
        tails[index] = max(tails[index], tails[index+1])
    return tuple(tails[h]-(tails[h+1] if h+1 < 4 else 0) for h in range(4))


def live_kernel_factor(tau, rank_charge, dimension_cap, live_coordinates, q=256):
    """Joint-diagonal rank potential from full-live kernel-vector counting.

    The common form space is fixed before ALL joint-matrix diagonals.
    This factor must not be multiplied by a charge using those same samples.
    """
    if not 1 <= tau <= q-1 or not 1 <= dimension_cap <= live_coordinates < q:
        raise ValueError('invalid joint live-kernel interface')
    return min(tau**(-rank_charge), live_kernel_bound(
        dimension_cap, live_coordinates, rank_charge, tau, q=q))


def polynomial_product(left, right, limit):
    if left.shape[1] == 2:
        return log_polynomial_product(left, right, limit)
    size = min(limit+1, len(left)+len(right)-1)
    states = left.shape[1]
    out = np.full((size, states, states), -math.inf)
    for degree, matrix in enumerate(left):
        count = min(len(right), size-degree)
        contribution = np.logaddexp.reduce(matrix[None, :, :, None]+right[:count, None, :, :], axis=2)
        out[degree:degree+count] = np.logaddexp(out[degree:degree+count], contribution)
    return out


def weighted_matrices(y, tau, rank_caps=False, live_kernel=True):
    if not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError('invalid marks')
    continuation, birth = rank_profiles()
    zero = y*0+tau*0
    if rank_caps:
        result = [[[zero for _ in range(5)] for _ in range(5)] for _ in range(9)]
        result[0][0][0] = y**0
        for w in range(9):
            count = math.comb(8, w)
            for (z, r, h, charge), multiplicity in continuation[w].items():
                factor = multiplicity*y**(8-z)*(1-y)**z/tau**r
                for cap in range(1, 5):
                    for surviving, mass in enumerate(root_cap_weights(tau, cap, h, charge), 1):
                        result[w][cap][min(4, surviving+w-r)] += factor*mass
            for (z, r), multiplicity in birth[w].items():
                result[w][0][min(4, w-r)] += multiplicity*y**(8-z)*(1-y)**z/tau**r
            if w:
                for cap in range(1, 5):
                    returns = zero
                    for k in range(min(w-1, cap-1)+1):
                        charge = (live_kernel_factor(tau, 4, cap-k+w, 4-k+w) if live_kernel else tau**-4)
                        returns += math.comb(4, k)/tau**k*charge
                    result[w][cap][0] = count*y**RETURN_MINIMA[w-1]*returns
                if w > 4:
                    charge = live_kernel_factor(tau, 4, w, w) if live_kernel else tau**-4
                    result[w][0][0] = count*y**w*charge
        return result
    result = [[[zero, zero], [zero, zero]] for _ in range(9)]
    result[0][0][0] = y**0
    for w in range(9):
        count = math.comb(8, w)
        result[w][1][1] = sum(m*y**(8-z)*(1-y)**z/tau**r*root_potential(tau, h, c)
                              for (z, r, h, c), m in continuation[w].items())
        if w:
            marginal_birth = sum(m*y**(8-z)*(1-y)**z/tau**r for (z, r), m in birth[w].items())
            result[w][0][1] = min(count*y**BIRTH_MINIMA[w-1], marginal_birth)
            # The unchanged state-return equation still has numerical rank4.
            flats = sum(math.comb(4, k)/tau**k for k in range(min(w-1, 3)+1))
            result[w][1][0] = count*flats*y**RETURN_MINIMA[w-1]/tau**4
            if w > 4:
                # C X=0 makes the feedthrough correction zero, so Y=X exactly.
                result[w][0][0] = count*y**w/tau**4
    return result


def packet(logx, logy, logtau, rank_caps=False, live_kernel=True, matrices=None):
    x, y, tau = math.exp(logx), math.exp(logy), min(255., math.exp(logtau))
    local = np.asarray((matrices or weighted_matrices)(y, tau, rank_caps, live_kernel), dtype=float)
    result = np.einsum('w,wij->ij', np.asarray([x**w for w in range(9)]), local)
    with np.errstate(divide='ignore'):
        return np.log(result)


def middle_bound(point, weight, dimension=96, k_log=20, rank_caps=False, live_kernel=True, matrices=None):
    g = profile(128, dimension, k_log)
    logx, logy, logtau, logtilt = point
    logs = packet(logx, logy, logtau, rank_caps, live_kernel, matrices)
    shift = float(logs.max())
    steps = g['n']//8
    return (128*modal_log_penalty(g['rows'])
            +g['rows']*source_row_log(weight/g['n'], logtau, logtilt, 128, dimension, True)
            +steps*shift+log_terminal_power(np.exp(logs-shift), steps)
            -weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def sparse_bound(active, point, dimension=96, k_log=20, rank_caps=False, live_kernel=True, matrices=None):
    g = profile(128, dimension, k_log)
    if len(point) not in (2, 3):
        raise ValueError('expected output/rank marks and optional minimum-weight mark')
    logtilt, logtau = point[:2]
    loglambda = point[2] if len(point) == 3 else 0.
    if loglambda < 0:
        raise ValueError('minimum-weight mark requires lambda>=1')
    logy = -math.exp(logtilt)*active/g['n']
    tau = min(255., math.exp(logtau))
    local = np.asarray((matrices or weighted_matrices)(math.exp(logy), tau, rank_caps, live_kernel)[:min(8, active)+1], dtype=float)
    with np.errstate(divide='ignore'):
        factor = np.log(local)
    states = local.shape[1]
    identity = np.full((states, states), -math.inf)
    np.fill_diagonal(identity, 0.)
    coefficients = np.array([identity])
    steps = g['rows']//8
    while steps:
        if steps & 1:
            coefficients = polynomial_product(coefficients, factor, active)
        steps >>= 1
        if steps:
            factor = polynomial_product(factor, factor, active)
    regional = regional_from_coefficients(coefficients, active, g['rows'], math.log(tau)+loglambda)
    shift = float(regional.max())
    terminal = 128*shift+log_terminal_power(np.exp(regional-shift), 128)
    return (log_choose(g['rows'], active)-((128-dimension)*active+1)*math.log(tau)
            -(129-dimension)*active*loglambda-g['threshold']*logy+terminal)/math.log(2)


def search(value, sparse=False, dimension=96, k_log=20, iterations=65, start=None,
           rank_caps=False, live_kernel=True, matrices=None, source_tilt=False):
    g = profile(128, dimension, k_log)
    if sparse:
        target = int(value)
        objective = lambda p: sparse_bound(target, p, dimension, k_log, rank_caps, live_kernel, matrices)
        starts = [(math.log(640), math.log(255)), (math.log(640), math.log(8))]
        bounds = [(-2, 10), (0, math.log(255))]
        if source_tilt:
            starts = [point+(0.,) for point in starts]+[(math.log(640), math.log(8), math.log(2))]
            bounds.append((0, 8))
    else:
        target = max(1, round(value*g['n']))
        objective = lambda p: middle_bound(p, target, dimension, k_log, rank_caps, live_kernel, matrices)
        starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
        bounds = [(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)]
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(objective, initial, method='L-BFGS-B', bounds=bounds,
                        options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25}) for initial in starts]
    best = min(choices, key=lambda c: c.fun if np.isfinite(c.fun) else math.inf)
    g.update(mode='active_rows' if sparse else 'input_weight', target=target,
             margin_bits=-float(best.fun), log_parameters=best.x.tolist(),
             optimizer_success=bool(best.success), iterations=int(best.nit), mu=4, eta=2,
             birth_minima=list(BIRTH_MINIMA), return_minima=list(RETURN_MINIMA),
             construction='gf256-cubic8/4-feedthrough-mu4-eta2-regional-nonzero-diagonals',
             rank_caps=rank_caps, live_kernel_return=live_kernel and rank_caps,
             minimum_source_weight_mark=source_tilt and sparse,
             scope='one weight class of the nonidentity-feedthrough positive cover; floating arithmetic')
    return g


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension', type=int, default=96)
    parser.add_argument('--density', type=float, nargs='*', default=[.03, .10, .14])
    parser.add_argument('--sparse', type=int, nargs='*', default=[64, 128])
    parser.add_argument('--iterations', type=int, default=65)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--rank-caps', action='store_true')
    parser.add_argument('--no-live-kernel', action='store_true')
    args = parser.parse_args()
    continuation, birth = rank_profiles()
    print(json.dumps({'rank_profiles_ready': True, 'continuation_classes': [len(x) for x in continuation],
                      'birth_classes': [len(x) for x in birth], 'rigorous': False}), flush=True)
    for sparse, values in ((True, args.sparse), (False, args.density)):
        previous = None
        for value in values:
            result = search(value, sparse, args.dimension, args.k_log, args.iterations, previous,
                            args.rank_caps, not args.no_live_kernel)
            previous = result['log_parameters']
            print(json.dumps(result), flush=True)

"""Compose proved local refinements; floating screening ONLY.

The independent directed sparse ladder does not call this module. A negative
screen must be replayed with directed local moments before being certified.
"""

from fractions import Fraction
from functools import lru_cache
from math import comb, exp2, log2

from hd_randomized_outer_state_screen import add, total
from hd_randomized_outer_state_profile import local_profiles, power as power3
from hd_randomized_outer_flat_moment import flat_moments
from hd_randomized_outer_joint_quiet import joint_moment_logs


@lru_cache(None)
def occupancy_counts():
    counts = {(0, 0):1}
    for unused in range(8):
        following = {}
        for (j, u), value in counts.items():
            for a in range(9):
                key = j+a, u+(a != 0)
                following[key] = following.get(key, 0)+value*comb(8, a)
        counts = following
    return tuple(tuple((u, log2(value)-log2(comb(64, j)))
                       for (jj, u), value in counts.items() if jj == j)
                 for j in range(65))


@lru_cache(None)
def refined_profiles(log_z, *, joint=True, flat=True, returns=True):
    q, d = 65536, 65536**4-1
    epsilon, c = 1/q, 1-1/q
    mark = Fraction(2**float(log_z))
    base = local_profiles(log_z)
    result = [[list(row) for row in packet] for packet in base]
    quiet = joint_moment_logs(mark)[0] if joint else None
    density = flat_moments(mark) if flat else None
    t = add(-16, log2(c)+log_z)
    for j in range(65):
        if joint:
            result[j][0][0] = min(result[j][0][0], quiet[j])
        if flat:
            uniform = density[j]
            live = min(base[j][2][2]+16, log2(d)+uniform)
            gamma = total(weight+min(j*t, max(0, 64-8*u-3)*log_z)
                          for u, weight in occupancy_counts()[j])
            gamma = min(gamma, log2(d)+uniform)
            # A cap for the nonzero target columns also bounds the zero
            # column after restricting h!=0. Keep the older stronger cap too.
            bulk = add(log2(c)+uniform, -16+gamma)
            result[j][1][1] = min(result[j][1][1], bulk)
            if j:
                result[j][1][0] = min(result[j][1][0], bulk-log2(d))
            result[j][2][1] = min(result[j][2][1], log2(c)+live)
            result[j][2][2] = min(result[j][2][2], -16+live)
            if j:
                result[j][2][0] = min(result[j][2][0], 48-log2(d)+live)
        if returns and j:
            rank = min(j, 4)
            factor = log2(epsilon+c*(q**rank-1)/d)
            untouched = total(weight+max(0, 64-8*u-3)*log_z
                              for u, weight in occupancy_counts()[j])
            result[j][2][0] = min(result[j][2][0], untouched-16*rank+factor)
            if joint:
                result[j][2][0] = min(result[j][2][0], quiet[j]+factor)
    return tuple(tuple(tuple(row) for row in packet) for packet in result)


def marker_bound(rows, log_z, alpha, marker, *, joint=True, flat=True, returns=True):
    operators = refined_profiles(log_z, joint=joint, flat=flat, returns=returns)
    weights = [log2(comb(64, j))+j*log2(marker)+(64-j)*log2(1-marker) for j in range(65)]
    matrix = tuple(tuple(total(weights[j]+alpha*operators[j][a][b] for j in range(65))
                         for b in range(3)) for a in range(3))
    run = power3(matrix, 32768)
    choose = log2(comb(8192, rows))
    normalizer = choose+rows*log2(marker)+(8192-rows)*log2(1-marker)
    return (choose+alpha*(rows*(4096-128*log2(65535))-log2(65535)-838860*log_z)
            +total(run[0])-256*normalizer)


def screen(rows):
    best = (float('inf'), None)
    for log_z in (-1/512, -1/64, -.125, -.5, -1, -2, -4, -8, -12, -14, -16):
        for alpha in (1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1):
            for factor in (.5, 1, 2, 4):
                marker = rows/8192*factor
                if not 0 < marker < 1:
                    continue
                value = marker_bound(rows, log_z, alpha, marker)
                if value < best[0]:
                    best = value, dict(log_z=log_z, alpha=alpha, marker=marker)
    return dict(active_rows=rows, diagnostic_log2=best[0], witness=best[1],
                directed=False, whole_code_certificate=False)


def optimized_screen(rows):
    """Same bounded tilt/exponent grid, with a continuous marker search."""
    from scipy.optimize import minimize_scalar
    best = (float('inf'), None)
    for log_z in (-1/512, -1/64, -.125, -.5, -1, -2, -4, -8, -12, -14, -16):
        for alpha in (1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1):
            def objective(log_odds):
                return marker_bound(rows, log_z, alpha, 1/(1+exp2(-log_odds)))
            result = minimize_scalar(objective, bounds=(-24, 24), method='bounded',
                                     options={'xatol':1e-7})
            if result.success and result.fun < best[0]:
                best = result.fun, dict(log_z=log_z, alpha=alpha,
                                       marker=1/(1+exp2(-result.x)))
    return dict(active_rows=rows, diagnostic_log2=best[0], witness=best[1],
                optimized_marker=True, directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', nargs='+', type=int, default=[16, 128, 512, 1024, 2048, 4096, 6144])
    parser.add_argument('--optimize-marker', action='store_true')
    arguments = parser.parse_args()
    for rows in arguments.rows:
        print(json.dumps((optimized_screen if arguments.optimize_marker else screen)(rows)), flush=True)

"""Exploratory half-rate systematic Toeplitz outer screen, not a certificate.

Each block maps x in GF16^k to (x,T*x), for an independent uniform k-by-k
Toeplitz T. For every nonzero fixed x, T*x is uniform in GF16^k. The
per-mask expected multiplicity is at most 16^-k * 15^w. The same positive
envelope holds for an independent fully random [2k,k] outer block.
No timing or native encoder is included.
"""

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp

from gf16_sizes import local_model
from reduced_distance_certificate import log_polynomial_product


def screen(k_log, active, state, dimension, z):
    rows = (1 << k_log) // dimension
    outer = 2 * dimension
    length = outer * rows
    threshold = (length + 9) // 10 - 1
    if (1 << k_log) % dimension or rows % 16:
        raise ValueError("geometry must divide exactly")
    with np.errstate(divide="ignore"):
        factors = np.log(local_model(z, state)[:active + 1])
    for w in range(len(factors)):
        factors[w] += math.log(math.comb(16, w))
    regions = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    power = rows // 16
    while power:
        if power & 1:
            regions = log_polynomial_product(regions, factors, active)
        power >>= 1
        if power:
            factors = log_polynomial_product(factors, factors, active)
    for h in range(active + 1):
        regions[h] -= math.lgamma(rows + 1) - math.lgamma(h + 1) - math.lgamma(rows - h + 1)
    source = np.array([math.log(math.comb(active, h)) + h * math.log(15)
                       - active * .5 * math.log(16) for h in range(active + 1)])
    matrix = logsumexp(regions + source[:, None, None], axis=0)
    vector = np.array([0., -math.inf])
    for _ in range(outer):
        vector = logsumexp(vector[:, None] + matrix, axis=0)
    moment = float(logsumexp(vector))
    choose = math.lgamma(rows + 1) - math.lgamma(active + 1) - math.lgamma(rows - active + 1)
    return -(choose - math.log(15) - threshold * math.log(z) + moment) / math.log(2)


def search(k_log, active, state, dimension):
    length = 2 << k_log
    def objective(logtilt):
        z = math.exp(-math.exp(logtilt) * active / length)
        return -screen(k_log, active, state, dimension, z)
    optimum = minimize_scalar(objective, bounds=(1, 10), method="bounded", options={"xatol": 2e-5})
    return {"k_log": k_log, "active_rows": active, "state": state, "outer": [2 * dimension, dimension],
            "margin_bits": -float(optimum.fun),
            "z": math.exp(-math.exp(optimum.x) * active / length), "rigorous": False}


def middle_screen(k_log, active, state, dimension, x, z):
    from coupled_chernoff import log_terminal_power
    rows = (1 << k_log) // dimension
    outer = 2 * dimension
    threshold = ((2 << k_log) + 9) // 10 - 1
    local = local_model(z, state)
    matrix = sum(math.comb(16, w) * (15*x)**w * (1+x)**(16-w) * local[w]
                 for w in range(17))
    choose = math.lgamma(rows + 1) - math.lgamma(active + 1) - math.lgamma(rows - active + 1)
    value = (-(outer-1)*choose - dimension*active*math.log(16) - math.log(15)
             -threshold*math.log(z) - outer*active*math.log(x)
             + log_terminal_power(matrix, rows*outer//16))
    return -value / math.log(2)


def middle_search(k_log, active, state, dimension):
    rows = (1 << k_log) // dimension
    def objective(point):
        logx, logz = point
        if not -35 < logx < 10 or not -20 < logz < 0:
            return 1e100
        return -middle_screen(k_log, active, state, dimension, math.exp(logx), math.exp(logz))
    best = None
    for initial_z in [.5, .9, .99, .999]:
        initial_x = active / (max(1, rows-active)*16)
        optimum = minimize(objective, [math.log(initial_x), math.log(initial_z)], method="Nelder-Mead",
                           options={"maxiter": 300, "xatol": 1e-9, "fatol": 1e-7})
        candidate = {"k_log": k_log, "active_rows": active, "state": state,
                     "outer": [2*dimension, dimension], "margin_bits": -float(optimum.fun),
                     "x": math.exp(optimum.x[0]), "z": math.exp(optimum.x[1]), "rigorous": False}
        if best is None or candidate["margin_bits"] > best["margin_bits"]:
            best = candidate
    return best


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--active", nargs="+", type=int, default=[1, 4, 16, 64, 128])
    parser.add_argument("--state", nargs="+", type=int, default=[4, 6, 8])
    parser.add_argument("--dimension", nargs="+", type=int, default=[16])
    parser.add_argument("--k-log", type=int, default=20)
    parser.add_argument("--middle", action="store_true")
    args = parser.parse_args()
    for dimension in args.dimension:
        for state in args.state:
            for active in args.active:
                method = middle_search if args.middle else search
                print(json.dumps(method(args.k_log, active, state, dimension)), flush=True)

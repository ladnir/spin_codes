"""Exploratory dense-source exponent for BAA outer packet-support tails.

This is a binary64 optimization screen, not a certified asymptotic bound.
It uses a positive 2x2 packet transfer and Chernoff coefficient bounds.
"""

import argparse
from math import exp, log, sqrt

import numpy as np
from scipy.optimize import minimize, minimize_scalar


def entropy(p: float) -> float:
    if p <= 0 or p >= 1:
        return 0.0
    return -p * log(p) - (1 - p) * log(1 - p)


def golay_exponent(alpha: float) -> float:
    def value(log_t: float) -> float:
        t = exp(log_t)
        poly = 1 + 759*t**8 + 2576*t**12 + 759*t**16 + t**24
        return log(poly) / 24 - alpha * log_t

    if alpha == 0 or alpha == 1:
        return 0.0
    fit = minimize_scalar(value, bounds=(-8, 8), method="bounded",
                          options={"xatol": 1e-11})
    return fit.fun


def accumulator_exponent(alpha: float, beta: float) -> float:
    if alpha / 2 > beta or beta > 1 - alpha / 2:
        return float("-inf")
    return (beta * entropy(alpha / (2 * beta))
            + (1 - beta) * entropy(alpha / (2 * (1 - beta)))
            - entropy(alpha))


def packet_transfer_radius(x: float, u: float) -> float:
    x2 = x*x
    x4 = x2*x2
    a = 1 + u * (6*x2 + x4)
    b = 4*u*x*(1+x2)
    c = x*(1+3*u) + 4*u*x*x2
    d = u*(1+6*x2+x4)
    return (a + d + sqrt((a-d)*(a-d) + 4*b*c)) / 2


def final_packet_tail_exponent(beta: float, rho: float) -> tuple[float, float, float]:
    """Upper exponent for V/R<=rho conditioned on c/b=beta."""
    def value(point: np.ndarray) -> float:
        log_x, log_u = point
        x, u = exp(log_x), exp(log_u)
        radius = packet_transfer_radius(x, u)
        return (log(radius) / 4 - beta*log_x
                - (rho/4)*log_u - entropy(beta))

    fit = minimize(value, np.array([0.0, -1.0]),
                   method="L-BFGS-B", bounds=[(-8, 8), (-20, 0)],
                   options={"ftol": 1e-12, "maxiter": 1000})
    return fit.fun, fit.x[0], fit.x[1]


def screen(rho: float, grid: int, min_alpha: float) -> None:
    if not 0 < rho < 1:
        raise ValueError("rho must be strictly between 0 and 1")
    if not 0 < min_alpha < 1:
        raise ValueError("minimum alpha must be strictly between 0 and 1")
    assert abs(packet_transfer_radius(1.0, 1.0) - 16.0) < 1e-12
    samples = []
    alphas = np.unique(np.concatenate((
        np.geomspace(min_alpha, 0.1, grid),
        np.linspace(0.1, 1 - 1e-6, grid),
    )))
    for alpha in alphas:
        g = golay_exponent(float(alpha))
        betas = np.linspace(alpha/2, 1-alpha/2, grid)
        for beta in betas:
            beta = float(beta)
            p = accumulator_exponent(float(alpha), beta)
            if p == float("-inf"):
                continue
            f, log_x, log_u = final_packet_tail_exponent(beta, rho)
            samples.append((g+p+f, float(alpha), beta, f, log_x, log_u))
    samples.sort(reverse=True)
    print({"rho": rho, "min_alpha": min_alpha, "grid": grid,
           "worst_sampled_points": samples[:8],
           "scope": "binary64 grid/optimization diagnostic; not a uniform proof"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho", type=float, default=0.16)
    parser.add_argument("--grid", type=int, default=12)
    parser.add_argument("--min-alpha", type=float, default=0.008)
    args = parser.parse_args()
    screen(args.rho, args.grid, args.min_alpha)

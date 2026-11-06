"""Binary64 saddle screen for the selected-outer Turbo Q1 exponent.

Splits Golay/first-accumulator, final packet occupancy, and inner
transfer. The output is a diagnostic grid, not a proof certificate.
"""

import argparse
from math import exp, log

import numpy as np
from scipy.optimize import minimize

from packet_baa_q1_screen import LIVE_BIT_DENSITY, STATE_NONZERO
from packet_tail_exponent import (
    accumulator_exponent, entropy, golay_exponent, packet_transfer_radius,
)


def source_exponent(beta: float, samples: int) -> tuple[float, float]:
    upper_alpha = min(1.0, 2*beta, 2*(1-beta))
    values = []
    for alpha in np.linspace(0.0, upper_alpha, samples):
        alpha = float(alpha)
        values.append((golay_exponent(alpha)
                       + accumulator_exponent(alpha, beta), alpha))
    return max(values)


def packet_coefficient_exponent(beta: float, rho: float) -> tuple[float, float, float]:
    def value(point: np.ndarray) -> float:
        log_x, log_u = point
        radius = packet_transfer_radius(exp(log_x), exp(log_u))
        return (log(radius)/4 - beta*log_x
                -rho*log_u/4 - entropy(beta))

    candidates = []
    for initial in ([0.0, 0.0], [0.0, -1.0], [0.0, 1.0]):
        fit = minimize(value, np.asarray(initial), method="L-BFGS-B",
                       bounds=[(-30, 30), (-30, 30)])
        candidates.append((fit.fun, fit.x[0], fit.x[1]))
    return min(candidates)


def inner_radius(y: float, gamma: float) -> float:
    E = exp(-gamma)
    J = -np.expm1(-gamma)/gamma
    a = 1.0
    b = y*J
    c = y*J/STATE_NONZERO
    d = E*(1+y*(1-1/STATE_NONZERO))
    return (a+d+((a-d)*(a-d)+4*b*c)**0.5)/2


def inner_exponent(rho: float, delta: float) -> tuple[float, float, float]:
    def value(point: np.ndarray) -> float:
        log_y, log_gamma = point
        gamma = exp(log_gamma)
        return (log(inner_radius(exp(log_y), gamma))
                -rho*log_y -entropy(rho)
                +gamma*delta/LIVE_BIT_DENSITY)/4

    candidates = []
    for initial in ([-1.0, 0.0], [0.0, 0.0], [1.0, 0.0]):
        fit = minimize(value, np.asarray(initial), method="L-BFGS-B",
                       bounds=[(-20, 20), (-8, 4)])
        candidates.append((fit.fun, fit.x[0], exp(fit.x[1])))
    return min(candidates)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho-floor", type=float, default=0.19)
    parser.add_argument("--grid", type=int, default=16)
    parser.add_argument("--source-grid", type=int, default=48)
    parser.add_argument("--delta", type=float, default=0.11)
    parser.add_argument("--schedule", type=float, default=39/4)
    args = parser.parse_args()
    if not (0 < args.rho_floor < 1 and args.grid >= 2):
        raise ValueError("invalid packet floor or grid")

    samples = []
    betas = np.linspace(0.0001, 0.9999, args.grid)
    rhos = np.linspace(args.rho_floor, 1, args.grid)
    inner = {float(rho): inner_exponent(float(rho), args.delta)
             for rho in rhos}
    for beta in betas:
        beta = float(beta)
        source, alpha = source_exponent(beta, args.source_grid)
        for rho in rhos:
            rho = float(rho)
            outer, log_x, log_u = packet_coefficient_exponent(beta, rho)
            inner_value, log_y, gamma = inner[rho]
            total = (log(2)/args.schedule + source + outer + inner_value)
            samples.append((total, alpha, beta, rho, source, outer,
                            inner_value, gamma, log_x, log_u, log_y))
    samples.sort(reverse=True)
    print({"rho_floor": args.rho_floor, "grid": args.grid,
           "source_grid": args.source_grid,
           "worst_samples": samples[:10],
           "scope": "binary64 grid saddle screen; not a certificate"})


if __name__ == "__main__":
    main()

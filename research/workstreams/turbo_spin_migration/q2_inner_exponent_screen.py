"""Continuum two-block Turbo inner exponent diagnostic, not a certificate."""

import argparse
from math import exp, log

import numpy as np
from scipy.optimize import minimize

from packet_baa_q1_screen import LIVE_BIT_DENSITY, STATE_NONZERO
from packet_tail_exponent import entropy


def kernels(gamma: float):
    E = exp(-gamma)
    J = -np.expm1(-gamma)/gamma
    H = 2*(gamma-1+E)/(gamma*gamma)
    G = 2*(1-(1+gamma)*E)/(gamma*gamma)
    r = 1/STATE_NONZERO
    k0 = np.array([[1.0, 0.0], [0.0, E]])
    k1 = np.array([[0.0, J], [r*J, (1-r)*E]])
    k2 = np.array([[r*H, (1-r)*G],
                   [r*(1-r)*G, r*G+(1-r)**2*E]])
    return k0, k1, k2


def perron(matrix: np.ndarray) -> float:
    a, b = matrix[0]
    c, d = matrix[1]
    return (a+d+((a-d)**2+4*b*c)**0.5)/2


def screen(rho1: float, rho2: float, delta: float, schedule: float):
    def value(point: np.ndarray) -> float:
        log_u, log_v, log_gamma = point
        u, v, gamma = exp(log_u), exp(log_v), exp(log_gamma)
        k0, k1, k2 = kernels(gamma)
        radius = perron(k0+(u+v)*k1+u*v*k2)
        inner = (log(radius)-rho1*log_u-rho2*log_v
                 -entropy(rho1)-entropy(rho2)
                 +gamma*delta/LIVE_BIT_DENSITY)/4
        return 2*log(2)/schedule + inner

    fits = []
    for start in ([0, 0, 0], [1, 1, 0], [-1, -1, 0]):
        fit = minimize(value, np.asarray(start, dtype=float),
                       method="L-BFGS-B",
                       bounds=[(-16, 16), (-16, 16), (-8, 4)])
        fits.append(fit)
    best = min(fits, key=lambda fit: fit.fun)
    return {"rho1": rho1, "rho2": rho2,
            "schedule_plus_inner_nats_per_bit": best.fun,
            "u": exp(best.x[0]), "v": exp(best.x[1]),
            "gamma": exp(best.x[2])}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho1", type=float, default=0.19)
    parser.add_argument("--rho2", type=float, default=0.19)
    parser.add_argument("--delta", type=float, default=0.11)
    parser.add_argument("--schedule", type=float, default=39/4)
    args = parser.parse_args()
    print(screen(args.rho1, args.rho2, args.delta, args.schedule))

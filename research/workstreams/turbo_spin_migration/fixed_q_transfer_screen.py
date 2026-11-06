"""Binary64 equal-packet-fraction screen for fixed occupation Q>=3.

This checks whether the exact continuum transfer has room before trying a
uniform certificate. It does not handle unequal packet counts or growing Q.
"""

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize


M = (1 << 22)-1
P0 = (1 << 21)/M
SCHEDULE = 39/4


def entropy(x):
    return 0 if x in (0, 1) else -x*math.log(x)-(1-x)*math.log1p(-x)


def kernels(q, gamma):
    d = np.diag([0., -gamma])
    p = np.array([[0., 1.], [1/M, 1-1/M]])
    matrix = np.zeros((2*(q+1), 2*(q+1)))
    for i in range(q+1):
        matrix[2*i:2*i+2, 2*i:2*i+2] = d
        if i < q:
            matrix[2*i:2*i+2, 2*i+2:2*i+4] = p
    transfer = expm(matrix)
    return [math.factorial(a)*transfer[:2, 2*a:2*a+2]
            for a in range(q+1)]


def screen(q, rho, outer):
    def objective(point):
        log_u, log_gamma = point
        u, gamma = math.exp(log_u), math.exp(log_gamma)
        k = kernels(q, gamma)
        t = sum(math.comb(q, a)*u**a*k[a] for a in range(q+1))
        eig = np.linalg.eigvals(t)
        radius = max(eig.real)
        return (q*(outer+math.log(2)/SCHEDULE)
                +(math.log(radius)-q*rho*log_u-q*entropy(rho)
                  +gamma*.11/P0)/4)

    fits = []
    for start in ((0, 0), (1, 1), (2, 2)):
        fits.append(minimize(objective, start, method="L-BFGS-B",
                             bounds=[(-20, 20), (-10, 6)]))
    fit = min(fits, key=lambda item: item.fun)
    return {"q": q, "rho": rho, "outer": outer,
            "total_exponent": fit.fun, "per_row": fit.fun/q,
            "u": math.exp(fit.x[0]), "gamma": math.exp(fit.x[1])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", type=int, nargs="+", default=[3, 4, 8])
    parser.add_argument("--state-bits", type=int, default=22)
    args = parser.parse_args()
    global M, P0
    M = (1 << args.state_bits)-1
    P0 = (1 << (args.state_bits-1))/M
    bands = json.loads(Path(__file__).with_name("packet_envelope_32.json").read_text())["bands"]
    for q in args.q:
        for index in (0, 8, 16, 24, 31):
            band = bands[index]
            rho = (float(Fraction(band["rho_lo"]))+float(Fraction(band["rho_hi"])))/2
            print(screen(q, rho, band["upper_nats_per_bit"]), flush=True)


if __name__ == "__main__":
    main()

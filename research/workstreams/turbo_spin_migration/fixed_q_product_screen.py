"""Binary64 diagnostic for the packet-route fixed-Q product norm."""

import json
import math
from fractions import Fraction
from pathlib import Path
from scipy.optimize import minimize_scalar


def entropy(x):
    return 0.0 if x in (0.0, 1.0) else -x*math.log(x)-(1-x)*math.log1p(-x)


def main():
    data = json.loads(Path(__file__).with_name("packet_envelope_32.json").read_text())
    sigma, tau, v = 127/250, 133/125, 3/1600
    r = 1/((1 << 22)-1)
    p0 = (1 << 21)/((1 << 22)-1)
    g3 = -minimize_scalar(
        lambda y: tau*y-(4/3)*math.log(1-y+y/sigma),
        bounds=(0, 1), method="bounded").fun
    outcomes = []
    for i, row in enumerate(data["bands"]):
        lo, hi = (float(Fraction(row[name])) for name in ("rho_lo", "rho_hi"))

        def objective(log_u):
            u = math.exp(log_u)
            m = max(1+u*sigma*v, u*r/v+sigma*(1+u))
            return (row["upper_nats_per_bit"]
                    +max((-entropy(x)-x*log_u+math.log(m))/4
                         for x in (lo, hi))
                    +g3/4+0.11*tau/(4*p0)+math.log(2)/(39/4))

        fit = minimize_scalar(objective, bounds=(-30, 30), method="bounded")
        outcomes.append((fit.fun, i, fit.x))
    print({"g3": g3, "worst_bands": sorted(outcomes, reverse=True)[:10],
           "scope": "binary64 diagnostic, not a certificate"})


if __name__ == "__main__":
    main()

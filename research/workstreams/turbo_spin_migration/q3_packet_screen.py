"""Binary64 three-position packet-transfer screen with unequal supports."""

import itertools
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


def kernels(gamma):
    e = math.exp(-gamma)
    i0 = -math.expm1(-gamma)/gamma
    i1 = (1-(1+gamma)*e)/gamma**2
    i2 = (2-(2+2*gamma+gamma**2)*e)/gamma**3
    weights = [3*(i0-2*i1+i2), 6*(i1-i2), 3*i2]
    r = 1/M
    p = np.array([[0., 1.], [r, 1-r]])
    k = np.zeros((2, 2))
    for path in itertools.product((0, 1), repeat=4):
        prob = math.prod(p[path[i], path[i+1]] for i in range(3))
        count = sum(path)
        factor = (1. if count == 0 else e if count == 4
                  else weights[count-1])
        k[path[0], path[-1]] += prob*factor
    # The lower kernels use the already checked two-packet formulas.
    j = i0
    h = 2*(i0-i1)
    g = 2*i1
    return (np.diag([1., e]),
            np.array([[0., j], [r*j, (1-r)*e]]),
            np.array([[r*h, (1-r)*g],
                      [r*(1-r)*g, r*g+(1-r)**2*e]]), k)


def screen(indices, bands):
    rows = [bands[index] for index in indices]
    rho = [(float(Fraction(row["rho_lo"]))+float(Fraction(row["rho_hi"])))/2
           for row in rows]
    outer = [row["upper_nats_per_bit"] for row in rows]

    def objective(point):
        u = [math.exp(x) for x in point[:3]]
        gamma = math.exp(point[3])
        k = kernels(gamma)
        t = k[0]+sum(u)*k[1]
        t += sum(u[i]*u[j] for i in range(3) for j in range(i+1, 3))*k[2]
        t += math.prod(u)*k[3]
        a, b, c, d = t.ravel()
        radius = (a+d+math.sqrt((a-d)**2+4*b*c))/2
        return (sum(outer)+3*math.log(2)/SCHEDULE
                +(math.log(radius)-sum(rho[i]*point[i] for i in range(3))
                  -sum(entropy(x) for x in rho)+gamma*.11/P0)/4)

    fit = minimize(objective, [0., 0., 0., 1.], method="L-BFGS-B",
                   bounds=[(-16, 16)]*3+[(-8, 5)])
    return fit.fun, tuple(indices), tuple(rho), tuple(fit.x)


def main():
    bands = json.loads(Path(__file__).with_name("packet_envelope_32.json").read_text())["bands"]
    points = [0, 4, 8, 12, 16, 20, 24, 28, 31]
    outcomes = [screen(indices, bands)
                for indices in itertools.combinations_with_replacement(points, 3)]
    outcomes.sort(reverse=True)
    print({"triples_checked": len(outcomes), "worst": outcomes[:12],
           "scope": "binary64, sampled bands at their midpoints; not certificate"})


if __name__ == "__main__":
    main()

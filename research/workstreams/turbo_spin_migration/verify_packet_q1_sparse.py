"""Outward check of the sparse-source Turbo Q1 exponent at 19% support."""

import math
from fractions import Fraction

from mpmath import iv


iv.dps = 80
ETA = Fraction(1, 50000)
RHO = Fraction(19, 100)
DELTA = Fraction(11, 100)
SCHEDULE = Fraction(39, 4)
P0 = Fraction(1 << 21, (1 << 22) - 1)
M = (1 << 22) - 1
Y = Fraction(6)
GAMMA = Fraction(2)


def as_interval(value: Fraction):
    return iv.mpf(value.numerator)/value.denominator


def upper(value) -> float:
    return math.nextafter(float(value.b), math.inf)


def entropy(value: Fraction):
    if value == 0 or value == 1:
        return iv.mpf(0)
    x = as_interval(value)
    return -x*iv.log(x)-(1-x)*iv.log(1-x)


def main() -> None:
    eta, y, gamma = map(as_interval, (ETA, Y, GAMMA))
    E = iv.exp(-gamma)
    J = (1-E)/gamma
    a, b = iv.mpf(1), y*J
    c = y*J/M
    d = E*(1+y*(1-iv.mpf(1)/M))
    radius = (a+d+iv.sqrt((a-d)**2+4*b*c))/2

    # C(b/24,j)4095^j <= (4095 e b/(24j))^j. Its normalized
    # exponent increases through j/b<=eta.
    source = eta*iv.log(4095*iv.e/(24*eta))
    values = []
    for rho in (RHO, Fraction(1)):
        r = as_interval(rho)
        value = (iv.log(iv.mpf(2))/as_interval(SCHEDULE) + source
                 + (iv.log(radius)-r*iv.log(y)-entropy(rho)
                    +gamma*as_interval(DELTA)/as_interval(P0))/4)
        values.append(upper(value))
    # The rho-dependent expression is convex; its maximum on [RHO,1]
    # occurs at an endpoint.
    assert max(values) < 0
    print({"rho_floor": str(RHO), "eta": str(ETA),
           "witness_y": str(Y), "witness_gamma": str(GAMMA),
           "endpoint_uppers": values,
           "max_upper": max(values), "result": "closed"})


if __name__ == "__main__":
    main()

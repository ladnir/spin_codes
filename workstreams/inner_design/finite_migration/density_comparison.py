"""Outward Poisson-binomial/binomial density ratio with inactive rows.

Nonzero Bernoulli probabilities belong to [pmin,pmax] or equal one.
The mean of these Q probabilities belongs to [nlo,nhi]. See
SHORT_LENGTH_BOUNDS.md for the tilting and Fourier argument.
"""
from fractions import Fraction as F
import math

from flint import arb


def number(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def variance_term(p, x):
    return p * (1 - p) * x / (1 - p + p * x)**2


def variance_floor(ps, mean, xlo, xhi):
    """Two-support LP for one tilted Bernoulli variance at fixed mean."""
    points = [(p, min(variance_term(p, xlo), variance_term(p, xhi)))
              for p in sorted(set(ps))] + [(F(1), F(0))]
    candidates = [v for p, v in points if p == mean]
    for index, (p, v) in enumerate(points):
        for r, w in points[index + 1:]:
            if p <= mean <= r:
                candidates.append(((r - mean) * v + (mean - p) * w) / (r - p))
    assert candidates
    return min(candidates)


def factor(rows, lo, hi, nlo, nhi, ps):
    """Upper ratio for every integer Q and mean in the supplied rectangle."""
    assert type(rows) is int and 1 <= lo <= hi <= rows
    pmin, pmax = min(ps), max(ps)
    assert all(isinstance(v, F) for v in (nlo, nhi, *ps))
    assert 0 < pmin <= pmax < 1 and pmin <= nlo <= nhi < 1
    theta_lo, theta_hi = F(lo, rows) * nlo, F(hi, rows) * nhi
    assert 0 < theta_lo <= theta_hi < 1
    best = arb(1)
    for j in range(1, hi + 1):
        # j > Q has probability zero. j=0 and j=rows have ratio <=1
        # directly by concavity of log(1-p) and log(p).
        if j == rows:
            continue
        y = F(j, rows)
        xs = (y * (1 - theta_hi) / (theta_hi * (1 - y)),
              y * (1 - theta_lo) / (theta_lo * (1 - y)))
        # The lower convex envelope ends at (1,0), with all earlier values
        # positive, so it decreases with mean. Q>=lo and nu<=nhi suffice.
        variance = lo * variance_floor(ps, nhi, *xs)
        assert variance > 0
        # At fixed tilt, Jensen's pgf upper bound increases with Q.
        # For Q=hi, the pgf ratio is largest at theta nearest j/rows.
        theta = min(theta_hi, max(theta_lo, y))
        argument = 1 - y + F(rows, hi) * (y - theta)
        assert argument > 0
        log_ratio = (hi * number(argument).log()
                     + (rows - hi) * number(1 - y).log()
                     - rows * number(1 - theta).log())
        # Fourier inversion: max mass <= exp(-V) I_0(V). Use an
        # explicit product rather than relying on a scaled-function convention.
        v = number(variance)
        max_mass = (-v).exp() * v.bessel_i(0)
        denominator = number(math.comb(rows, j)) * number(y)**j * number(1 - y)**(rows - j)
        ratio = (log_ratio.exp() * max_mass / denominator).upper()
        best = max(best, ratio)
    return best.upper()

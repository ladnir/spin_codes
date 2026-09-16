"""Routing density with an upper bound on the fraction of all-one rows."""
from fractions import Fraction as F
from functools import lru_cache
import math

from flint import arb
import density_comparison as base


@lru_cache(maxsize=128)
def variance_vertices(ps, nlo, nhi, maximum):
    """Vertices of the mean slab with all-one mass in [0, maximum]."""
    assert 0 <= maximum < 1 and min(ps) <= nlo <= nhi < 1
    candidates = set()
    def add(alpha, p, r, mass):
        if 0 <= alpha <= maximum and 0 <= mass <= 1-alpha:
            mean = alpha + mass*p + (1-alpha-mass)*r
            if nlo <= mean <= nhi:
                candidates.add((p, r, mass, 1-alpha-mass))
    # On either alpha face, a mean endpoint uses at most two ordinary
    # probabilities. Interior means can attain an ordinary vertex.
    for alpha in (F(0), maximum):
        for p in ps:
            add(alpha, p, p, 1-alpha)
        for i, p in enumerate(ps):
            for r in ps[i+1:]:
                for nu in (nlo, nhi):
                    add(alpha, p, r, (nu-alpha-(1-alpha)*r)/(p-r))
    # Vertices away from the alpha faces use one ordinary probability and
    # the all-one probability, with a tight mean endpoint.
    for p in ps:
        for nu in (nlo, nhi):
            alpha = (nu-p)/(1-p)
            add(alpha, p, p, 1-alpha)
    assert candidates, 'Infeasible constrained mean'
    return tuple(sorted(candidates))


def factor(rows, lo, hi, nlo, nhi, ps, maximum):
    assert type(rows) is int and 1 <= lo <= hi <= rows
    assert all(type(v) is F for v in (nlo, nhi, maximum, *ps))
    assert 0 < min(ps) <= max(ps) < 1
    vertices = variance_vertices(tuple(sorted(set(ps))), nlo, nhi, maximum)
    theta_lo, theta_hi = F(lo, rows)*nlo, F(hi, rows)*nhi
    best = arb(1)
    for j in range(1, hi+1):
        if j == rows:
            continue
        y = F(j, rows)
        xs = (y*(1-theta_hi)/(theta_hi*(1-y)), y*(1-theta_lo)/(theta_lo*(1-y)))
        values = {p: min(base.variance_term(p, x) for x in xs) for p in ps}
        variance = lo*min(a*values[p]+b*values[r] for p, r, a, b in vertices)
        assert variance > 0
        theta = min(theta_hi, max(theta_lo, y))
        argument = 1-y+F(rows, hi)*(y-theta)
        assert argument > 0
        log_ratio = (hi*base.number(argument).log()+(rows-hi)*base.number(1-y).log()
                     -rows*base.number(1-theta).log())
        v = base.number(variance)
        mass = (-v).exp()*v.bessel_i(0)
        denominator = base.number(math.comb(rows, j))*base.number(y)**j*base.number(1-y)**(rows-j)
        best = max(best, (log_ratio.exp()*mass/denominator).upper())
    return best.upper()

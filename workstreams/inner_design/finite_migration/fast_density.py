"""Skip variance programs when a unit point-mass bound already suffices."""
from fractions import Fraction as F
import math

from flint import arb
import constant_density as frozen


def factor(rows,lo,hi,nlo,nhi,ps,maximum):
    assert type(rows) is int and 1 <= lo <= hi <= rows
    assert all(type(v) is F for v in (nlo,nhi,maximum,*ps))
    assert 0 < min(ps) <= max(ps) < 1
    vertices = frozen.variance_vertices(tuple(sorted(set(ps))),nlo,nhi,maximum)
    theta_lo,theta_hi = F(lo,rows)*nlo,F(hi,rows)*nhi
    best = arb(1)
    for j in range(1,hi+1):
        if j == rows:
            continue
        y = F(j,rows)
        theta = min(theta_hi,max(theta_lo,y))
        argument = 1-y+F(rows,hi)*(y-theta)
        assert argument > 0
        log_ratio = (hi*frozen.base.number(argument).log()+(rows-hi)*frozen.base.number(1-y).log()
                     -rows*frozen.base.number(1-theta).log())
        denominator = (frozen.base.number(math.comb(rows,j))*frozen.base.number(y)**j
                       *frozen.base.number(1-y)**(rows-j))
        prefactor = (log_ratio.exp()/denominator).upper()
        # Every tilted point mass is <=1. No variance calculation for this
        # j can raise the current maximum when this upper bound is smaller.
        if prefactor <= best:
            continue
        xs = (y*(1-theta_hi)/(theta_hi*(1-y)),y*(1-theta_lo)/(theta_lo*(1-y)))
        values = {p:min(frozen.base.variance_term(p,x) for x in xs) for p in ps}
        variance = lo*min(a*values[p]+b*values[r] for p,r,a,b in vertices)
        assert variance > 0
        v = frozen.base.number(variance)
        mass = (-v).exp()*v.bessel_i(0)
        best = max(best,(prefactor*mass).upper())
    return best.upper()

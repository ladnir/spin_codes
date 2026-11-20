"""Entrywise convex sequence dominating all region matrices.

For a convex sequence f, a Poisson-binomial sum with Q trials and mean Qp
has E f(J) <= E f(Bin(Q,p)). Pairwise averaging probabilities proves this:
at fixed pair sum, the expectation increases with their product times the
nonnegative second difference of f. Iterating averages yields equal p.
"""
from flint import arb


def convexify(region):
    result=[tuple(v.upper() for v in row) for row in region]
    for j in range(len(result)-3,-1,-1):
        result[j]=tuple(max(v,(2*a-b).upper()) for v,a,b in zip(result[j],result[j+1],result[j+2]))
    return result


def assert_convex(region):
    # All entries are exact endpoints. Use exact dyadic rationals for checks.
    from fractions import Fraction as F
    def rational(value):
        m,e=value.man_exp()
        return F(int(m))*F(2)**int(e)
    points=[tuple(rational(v) for v in row) for row in region]
    assert all(a-2*b+c>=0 for x,y,z in zip(points,points[1:],points[2:]) for a,b,c in zip(x,y,z))

"""Directed fixed-point bounds on log2 of a positive rational.

Only integer arithmetic establishes bounds. No large integer powers, float
logs, decimal rounding assumptions, or external interval package are used.
"""
from functools import lru_cache
from fractions import Fraction

BITS=96
Q=1<<BITS
TERMS=64


def ceildiv(n,d):return -((-n)//d)


def ln_normalized(n,d):
    """Q-scaled bounds for ln(n/d), where 1<=n/d<=2."""
    assert d<=n<=2*d
    if n==d:return 0,0
    # ln(x)=2*atanh((x-1)/(x+1)). Bound each positive series term.
    zl=Q*(n-d)//(n+d);zu=ceildiv(Q*(n-d),n+d)
    z2l=zl*zl//Q;z2u=ceildiv(zu*zu,Q)
    tl,tu=zl,zu;sl=su=0
    for j in range(TERMS):
        sl+=tl//(2*j+1);su+=ceildiv(tu,2*j+1)
        tl=tl*z2l//Q;tu=ceildiv(tu*z2u,Q)
    # Actual z<=1/3<1/2. The tail of 2*sum z^(2j+1)/(2j+1)
    # after TERMS terms is <=4/(3*(2*TERMS+1)*4^TERMS).
    tail=ceildiv(4*Q,3*(2*TERMS+1)*(1<<(2*TERMS)))
    return 2*sl,2*su+tail


LN2_LO,LN2_HI=ln_normalized(2,1)
assert 0<LN2_LO<=LN2_HI<Q


@lru_cache(maxsize=8192)
def log2_bounds(value,scale=1<<24):
    value=Fraction(value);n,d=value.numerator,value.denominator
    assert n>0 and d>0 and isinstance(scale,int) and 0<scale<=Q
    exponent=n.bit_length()-d.bit_length()
    below=(n<d<<exponent) if exponent>=0 else (n<<(-exponent)<d)
    if below:exponent-=1
    if exponent>=0:d<<=exponent
    else:n<<=-exponent
    assert d<=n<2*d
    lo,hi=ln_normalized(n,d)
    return exponent*scale+scale*lo//LN2_HI,exponent*scale+ceildiv(scale*hi,LN2_LO)


def log2_upper_scaled(value,scale=1<<24):return log2_bounds(value,scale)[1]

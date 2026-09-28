"""Outward Holder envelopes for persistent packet-law types.

A finite family contains exact pairs (h_t, theta_t), where h_t>=0 and
theta_t is a probability law on packet weights 0..4. The type persists
across m regions. For each packet weight b, define

    S_b = sum_t h_t * theta_t(b)**m,    nu_b >= S_b**(1/m).

Holder's inequality gives, for every trajectory b_1,...,b_m,

    sum_t h_t * product_i theta_t(b_i) <= product_i nu_(b_i).

The helper rounds each root upward to an exact dyadic, then defines the
exact rationals z=sum_b nu_b and phi_b=nu_b/z. Thus the right side is
z**m * product_i phi_(b_i). This dominates a persistent-type mixture; it
does not replace persistent types by independently resampled types.

All inputs and final comparisons are exact. Arb proposes root enclosures;
integer-power checks establish the dyadic endpoint and its predecessor.
An empty or zero-mass family returns z=0 and phi concentrated at weight0.
The caller may discard that family because its dominating measure is zero.
"""
from fractions import Fraction as Q

from flint import arb, ctx


def _exact(value):
    # In particular, do not turn a binary64 proposal into a certificate
    # input by silently treating its represented value as an intended bound.
    if not isinstance(value,(int,Q)):
        raise TypeError('exact integer or Fraction input required')
    return Q(value)


def _parameters(m, precision=None):
    if not isinstance(m,int) or m<1:
        raise ValueError('positive integer trajectory length required')
    if precision is not None and (not isinstance(precision,int) or precision<32):
        raise ValueError('integer precision of at least 32 bits required')


def _power_of_two(exponent):
    return Q(1 << exponent) if exponent>=0 else Q(1,1 << -exponent)


def _fraction(point):
    value=point.fmpq()
    return Q(int(value.p),int(value.q))


def dyadic_root_bound(value, m, *, precision=192):
    """Return (upper, quantum), with upper the least valid grid endpoint.

    The positive quantum is a power of two on a precision-bit relative
    grid. Exact checks ensure upper**m>=value and (upper-quantum)**m<value.
    For value zero return (0,0). The grid itself is selected from an Arb
    enclosure; its particular scale is not part of the validity premise.
    """
    value=_exact(value)
    _parameters(m,precision)
    if value<0:
        raise ValueError('nonnegative root argument required')
    if not value:
        return Q(0),Q(0)
    previous=ctx.prec
    ctx.prec=precision
    try:
        root=(arb(value.numerator)/value.denominator).root(m)
        upper=root.upper()
        if not root.is_finite() or not upper>0:
            raise ArithmeticError('failed to obtain a finite positive root enclosure')
        mantissa,exponent=upper.man_exp()
        binary_exponent=int(mantissa).bit_length()-1+int(exponent)
        quantum=_power_of_two(binary_exponent-precision+1)
        scaled_upper=_fraction(upper)/quantum
        candidate=-(-scaled_upper.numerator//scaled_upper.denominator)
        target=value/quantum**m
        def valid(integer):
            return integer**m*target.denominator>=target.numerator
        if not valid(candidate):
            raise ArithmeticError('proposed root endpoint failed the exact upper-bound check')
        while candidate>0 and valid(candidate-1):
            candidate-=1
        result=candidate*quantum
        assert result**m>=value and (result-quantum)**m<value
        return result,quantum
    finally:
        ctx.prec=previous


def coordinate_sums(types, *, m=256):
    """Return the five exact S_b values; reject inexact input coefficients."""
    _parameters(m)
    sums=[Q(0)]*5
    for h,theta in types:
        h=_exact(h)
        theta=tuple(_exact(p) for p in theta)
        if h<0 or len(theta)!=5 or any(p<0 for p in theta) or sum(theta)!=1:
            raise ValueError('nonnegative h and five probabilities summing to one required')
        if h:
            for b,p in enumerate(theta):
                sums[b]+=h*p**m
    return tuple(sums)


def holder_envelope(types, *, m=256, precision=192):
    """Return exact (z,phi) dominating every length-m packet trajectory.

    types is an iterable of (h,theta) pairs with integer/Fraction entries.
    phi has five nonnegative entries summing exactly to one. z can be zero;
    in that case phi is the canonical point mass on the zero packet.
    """
    _parameters(m,precision)
    sums=coordinate_sums(types,m=m)
    nu=tuple(dyadic_root_bound(value,m,precision=precision)[0] for value in sums)
    z=sum(nu,Q(0))
    phi=tuple(value/z for value in nu) if z else (Q(1),Q(0),Q(0),Q(0),Q(0))
    assert z>=0 and all(p>=0 for p in phi) and sum(phi)==1
    # Constant trajectories are necessary constraints and independently
    # check the normalization of every returned coordinate.
    assert all((z*phi[b])**m>=sums[b] for b in range(5))
    return z,phi

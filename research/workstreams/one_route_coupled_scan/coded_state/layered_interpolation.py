"""Rational interpolation scale tau=2^m, allowing noninteger proof weights.

Integers retain the historical convention m in [0,32]. A LogWeight stores
tau itself exactly; its exponent is log2(tau), not a rounded real number.
"""
from dataclasses import dataclass
from fractions import Fraction
import math
from rational_log_bounds import log2_bounds


@dataclass(frozen=True)
class LogWeight:
    tau: Fraction

    def __post_init__(self):
        object.__setattr__(self,'tau',Fraction(self.tau))
        assert 1<=self.tau<=1<<32

    def __float__(self):return math.log2(self.tau)


def exact_tau(weight):
    if isinstance(weight,LogWeight):return weight.tau
    assert isinstance(weight,int) and 0<=weight<=32
    return Fraction(1<<weight)


def exponent_bounds(weight,scale):
    lo,hi=log2_bounds(exact_tau(weight),scale)
    # The route D/2 cap uses an integer half-coefficient.
    return lo,2*((hi+1)//2)

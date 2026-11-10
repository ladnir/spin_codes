"""Audited producer formulas with compact UNTRUSTED Map22 count candidates.

This deliberately uses no capped birth, inverse Walsh transform, or raw
state arrays. Kernel authentication of the compact counts remains required.
"""
from fractions import Fraction as Q
from math import comb
from flint import arb, arb_mat

LEVELS = [16, 24, 28, 32, 36, 40, 48, 64]
S = 1 << 22


def aq(x):
    x = Q(x)
    return arb(x.numerator) / x.denominator


def up(x):
    if not x.is_finite():
        raise ArithmeticError('finite enclosure required')
    v = x.upper()
    if v < 0:
        raise ArithmeticError('negative upper bound for nonnegative quantity')
    return v


def coefficients(profile, factors):
    values = [factors[0][0] * 0 + 1]
    for count, (inactive, active) in zip(profile, factors):
        for _ in range(count):
            new = [values[0] * 0 for _ in range(len(values) + 1)]
            for j, value in enumerate(values):
                new[j] += inactive * value
                new[j+1] += active * value
            values = new
    return values


class CoarsePhysical:
    def __init__(self, candidate):
        self.profiles = candidate['profiles']
        self.counts = candidate['multiplicities']
        self.births = {b['occupancy']: b for b in candidate['small_births']}
        if (sum(self.counts) != S or len(self.profiles) != 25 or
                any(sum(h) != 16 for h in self.profiles)):
            raise ValueError('complete candidate profile partition required')
        self.weights = [sum(r*n for r,n in enumerate(h)) for h in self.profiles]
        self.live = [k for k,w in enumerate(self.weights) if w]
        self.classes = [[k for k,w in enumerate(self.weights) if w == v] for v in LEVELS]
        feedback = [Q(0) for _ in range(17)]
        for h, count in zip(self.profiles, self.counts):
            polynomial = coefficients(h, [(Q(1), Q(1) if r == 0 else Q(-1, 15)) for r in range(5)])
            for j, v in enumerate(polynomial):
                feedback[j] += count * v / (S * comb(16, j))
        if any(not 0 <= x <= 1 for x in feedback) or feedback[0] != 1:
            raise ArithmeticError('exact feedback candidate invalid')
        self.feedback_zero = feedback

    def family(self, tilt):
        z = (-aq(tilt)).exp()
        divisors = [comb(16, j) for j in range(17)]
        moments, seconds = [], []
        for weight_z, target in [(z, moments), (z*z, seconds)]:
            factors = [(weight_z**r, ((1+weight_z)**4-weight_z**r)/15) for r in range(5)]
            for h in self.profiles:
                target.append([up(v/divisors[j]) for j,v in enumerate(coefficients(h, factors))])
        factors = [(arb(1), ((1+z)**(4-r)*(1-z)**r-1)/15) for r in range(5)]
        zero = [arb(0) for _ in range(17)]
        for h, count in zip(self.profiles, self.counts):
            for j,v in enumerate(coefficients(h, factors)):
                zero[j] += count*v/(S*divisors[j])
        g = ((1+z)**4-1)/15
        family = []
        for j in range(17):
            rows = [[arb(0) for _ in range(11)] for _ in range(11)]
            if j == 0:
                rows[0][0] = 1
            elif j in self.births:
                b = self.births[j]
                for i, polynomial in enumerate(b['coefficients']):
                    v = sum((n*z**w for w,n in enumerate(polynomial)), arb(0))/b['denominator']
                    rows[0][0 if i == 0 else i+2] = up(v)
            else:
                rows[0][0] = up(zero[j])
                rows[0][1] = up(g**j)
            nonzero = aq(1-self.feedback_zero[j])
            def refresh(i, first, second):
                rows[i][0] = arb(0) if j == 0 else up(min(first, up((second*nonzero).sqrt()))/(S-1))
                rows[i][2] = first
            refresh(1, max(moments[k][j] for k in self.live), max(seconds[k][j] for k in self.live))
            refresh(2, up(sum((self.counts[k]*moments[k][j] for k in self.live), arb(0))/(S-1)),
                up(sum((self.counts[k]*seconds[k][j] for k in self.live), arb(0))/(S-1)))
            for i,indices in enumerate(self.classes, 3):
                refresh(i, max(moments[k][j] for k in indices), max(seconds[k][j] for k in indices))
            family.append(arb_mat(rows))
        return family


def binary_target_floor(bits, precision=256):
    """Exact dyadic lower bound on 2**(-bits), using only integer powers."""
    bits = Q(bits)
    power = precision*bits.denominator-bits.numerator
    if power < 0:
        raise ValueError('precision must exceed requested bit margin')
    target = 1 << power
    lo, hi = 0, 1 << ((power+bits.denominator-1)//bits.denominator+1)
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid**bits.denominator <= target:
            lo = mid
        else:
            hi = mid
    return Q(lo, 1 << precision)

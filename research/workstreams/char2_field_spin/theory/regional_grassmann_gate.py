"""Exact fresh-orientation root kernel for the separate refresh construction.

This is a proof reference, not a native encoder or a distance certificate.
The incoming common quotient and all source forms are fixed before a fresh
uniform GL(s,q) orientation is sampled.  Its image is a uniform r-subspace.
Do not reuse this distribution after conditioning on this orientation's
zero-output tests.  See REGIONAL_GRASSMANN_GATE.md for the interface.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from math import isqrt, log2


def _integer(value, name, lower=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < lower:
        raise ValueError(f"{name} must be an integer at least {lower}")
    return value


def _field_order(q):
    _integer(q, "q", 2)
    p = next((d for d in range(2, isqrt(q) + 1) if q % d == 0), q)
    remainder = q
    while remainder % p == 0:
        remainder //= p
    if remainder != 1:
        raise ValueError("q must be a prime power")


@lru_cache(maxsize=None)
def gaussian_binomial(n, k, q):
    """The number of k-dimensional subspaces of F_q^n, exactly."""
    _field_order(q)
    _integer(n, "n")
    if not isinstance(k, int) or isinstance(k, bool):
        raise ValueError("k must be an integer")
    if k < 0 or k > n:
        return 0
    k = min(k, n-k)
    numerator = denominator = 1
    for i in range(k):
        numerator *= q**(n-i)-1
        denominator *= q**(k-i)-1
    value, remainder = divmod(numerator, denominator)
    assert remainder == 0
    return value


@lru_cache(maxsize=None)
def root_rank_distribution(q, s, r, z):
    """Return P[h=0..s] for h=dim(U intersect K), codim(K)=z.

    U is uniform among r-dimensional subspaces of F_q^s.  All probabilities
    are Fractions; the mass at h=0 is retained in this distribution.
    """
    _field_order(q)
    _integer(s, "s", 1)
    _integer(r, "r")
    _integer(z, "z")
    if r > s or z > s:
        raise ValueError("require 0 <= r,z <= s")
    result = [Fraction(0)]*(s+1)
    denominator = gaussian_binomial(s, r, q)
    for h in range(max(0, r-z), min(r, s-z)+1):
        count = (gaussian_binomial(s-z, h, q)
                 * gaussian_binomial(z, r-h, q)
                 * q**((s-z-h)*(r-h)))
        result[h] = Fraction(count, denominator)
    assert sum(result) == 1
    return tuple(result)


def zero_phase_action(q, s, r, z, f0, d0, tau, continuation, *, fresh=True):
    """Exact root average followed by a conservative diagonal action.

    This interface requires an incoming *nonzero witness* state, 1<=r<=s.
    z counts zero-input/zero-output positions.  f0/d0 count fresh/dependent
    inputs forced to zero output.  continuation[h] bounds all subsequent
    operations uniformly conditional on their history.  Return value is an
    upper expected alive tau^(-equation-rank-increase) times continuation.

    fresh=False takes the prior worst-case root rank, for comparison only.
    Freshness is not resampled after the root equations.  The maximum over
    post-diagonal ranks deliberately retains arbitrary hidden dependencies.
    """
    _field_order(q)
    _integer(s, "s", 1)
    for value, name in ((z, "z"), (f0, "f0"), (d0, "d0")):
        _integer(value, name)
    _integer(r, "r", 1)
    tau = Fraction(tau)
    values = tuple(map(Fraction, continuation))
    if r > s or len(values) != s+1 or min(values) < 0 or not 1 <= tau <= q-1:
        raise ValueError("invalid rank, continuation, or tau")
    probabilities = root_rank_distribution(q, s, r, min(z, s))
    contributions = []
    for u in range(1, s+1):
        if not probabilities[u]:
            continue
        payoff = max(values[max(1, u-d0):u+1])
        candidate = payoff/tau**(r-u+f0+d0)
        contributions.append(probabilities[u]*candidate if fresh else candidate)
    return sum(contributions, Fraction(0)) if fresh else max(contributions, default=Fraction(0))


def alive_root_contraction(q, s, r, z, tau):
    return zero_phase_action(q, s, r, z, 0, 0, tau, [1]*(s+1))


def cost(refresh_symbols):
    """Conservative processed-byte arithmetic volumes, not runtime ratios.

    Charge a full additional four GFNI512 products at each refresh.  This
    ignores possible savings from replacing a pre-existing phase boundary.
    The once-per-region setup cost is already negligible in this comparison.
    """
    if refresh_symbols not in (16, 64):
        raise ValueError("supported diagnostic refresh lengths are 16 and 64")
    baseline = Fraction(181, 88)
    extra = Fraction(4, refresh_symbols)*Fraction(32, 22)
    return baseline, extra, baseline+extra, extra/baseline


def check():
    for q in (2, 4, 16, 256):
        for r in range(5):
            for z in range(5):
                probability = root_rank_distribution(q, 4, r, z)
                assert sum(probability) == 1 and min(probability) >= 0
                if z == 0:
                    assert probability[r] == 1
                if z == 4:
                    assert probability[0] == 1
    ratios = []
    for r in range(1, 5):
        for z in range(1, 4):
            contraction = alive_root_contraction(256, 4, r, z, 255)
            ratio = contraction*255**z
            assert ratio <= Fraction(257, 256)
            ratios.append((ratio, r, z))
        assert alive_root_contraction(256, 4, r, 4, 255) == 0
    largest = max(ratios)
    # There is genuine positive slack; do not replace the exact kernel by
    # tau^-z without a multiplicative correction.
    assert alive_root_contraction(256, 4, 2, 1, 255) > Fraction(1, 255)
    return largest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        ratio, r, z = check()
        print(f"exact kernel checks passed; max E*255^z={float(ratio):.12f}, r={r}, z={z}")
    for r in range(1, 5):
        print("r", r, [(z, round(-log2(float(alive_root_contraction(256, 4, r, z, 255))), 9))
                       for z in range(1, 4)])
    for period in (64, 16):
        base, extra, total, overhead = cost(period)
        print(f"refresh {period}: {base}+{extra}={total}; arithmetic volume +{100*float(overhead):.6f}%")

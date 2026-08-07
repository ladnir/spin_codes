"""Small standard-library outward interval layer for base-two logarithms.

``Decimal.ln`` is correctly rounded by Python's decimal implementation.  We
evaluate it with guard digits, enlarge by a much larger explicit pad, and use
directed rounding for all subsequent arithmetic.  The resulting intervals are
intended for certificate comparisons with many bits of slack, not tight
numerical analysis.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from functools import lru_cache


PRECISION = 100
GUARD_PRECISION = 130
LN_PAD = Decimal("1e-110")


@dataclass(frozen=True)
class Interval:
    lo: Decimal
    hi: Decimal

    @staticmethod
    def exact(value: int | Decimal) -> "Interval":
        item = Decimal(value)
        return Interval(item, item)

    def __add__(self, other: "Interval") -> "Interval":
        with localcontext() as ctx:
            ctx.prec = PRECISION
            ctx.rounding = ROUND_FLOOR
            lo = self.lo + other.lo
            ctx.rounding = ROUND_CEILING
            hi = self.hi + other.hi
        return Interval(lo, hi)

    def __sub__(self, other: "Interval") -> "Interval":
        with localcontext() as ctx:
            ctx.prec = PRECISION
            ctx.rounding = ROUND_FLOOR
            lo = self.lo - other.hi
            ctx.rounding = ROUND_CEILING
            hi = self.hi - other.lo
        return Interval(lo, hi)

    def __mul__(self, other: "Interval") -> "Interval":
        with localcontext() as ctx:
            ctx.prec = PRECISION
            ctx.rounding = ROUND_FLOOR
            lo = min(
                self.lo * other.lo,
                self.lo * other.hi,
                self.hi * other.lo,
                self.hi * other.hi,
            )
            # Recompute under upward rounding rather than reusing products.
            ctx.rounding = ROUND_CEILING
            hi = max(
                self.lo * other.lo,
                self.lo * other.hi,
                self.hi * other.lo,
                self.hi * other.hi,
            )
        return Interval(lo, hi)

    def __truediv__(self, other: "Interval") -> "Interval":
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval divisor contains zero")
        with localcontext() as ctx:
            ctx.prec = PRECISION
            ctx.rounding = ROUND_FLOOR
            lo = min(
                self.lo / other.lo,
                self.lo / other.hi,
                self.hi / other.lo,
                self.hi / other.hi,
            )
            ctx.rounding = ROUND_CEILING
            hi = max(
                self.lo / other.lo,
                self.lo / other.hi,
                self.hi / other.lo,
                self.hi / other.hi,
            )
        return Interval(lo, hi)

    def times_int(self, value: int) -> "Interval":
        return self * Interval.exact(value)


def _ln_decimal_bounds(value: Decimal) -> Interval:
    if value <= 0:
        raise ValueError("logarithm requires a positive value")
    with localcontext() as ctx:
        ctx.prec = GUARD_PRECISION
        middle = value.ln()
        # The pad is over ten orders of magnitude larger than the maximum
        # rounding unit at this guard precision for every value used here.
        lo = middle - LN_PAD
        hi = middle + LN_PAD
    return Interval(lo, hi)


@lru_cache(maxsize=None)
def ln_int(value: int) -> Interval:
    return _ln_decimal_bounds(Decimal(value))


@lru_cache(maxsize=None)
def ln_fraction(value: Fraction) -> Interval:
    return ln_int(value.numerator) - ln_int(value.denominator)


LN2 = ln_int(2)


def log2_int(value: int) -> Interval:
    return ln_int(value) / LN2


def log2_fraction(value: Fraction) -> Interval:
    return ln_fraction(value) / LN2


# Classical Archimedean rational bounds 333/106 < pi < 355/113 give
# 333/53 < 2*pi < 710/113.  Robbins' factorial bounds below therefore need no
# floating pi constant.
LN_TWO_PI = Interval(
    ln_fraction(Fraction(333, 53)).lo,
    ln_fraction(Fraction(710, 113)).hi,
)


@lru_cache(maxsize=None)
def ln_factorial(n: int) -> Interval:
    """Robbins interval for ln(n!), with exact small factorial fallback."""

    if n < 0:
        raise ValueError("factorial argument must be nonnegative")
    if n <= 256:
        return ln_int(math.factorial(n))
    log_n = ln_int(n)
    half = Interval.exact(Decimal("0.5"))
    center = log_n * (Interval.exact(n) + half) - Interval.exact(n) + LN_TWO_PI * half
    lower_correction = Interval.exact(1) / Interval.exact(12 * n + 1)
    upper_correction = Interval.exact(1) / Interval.exact(12 * n)
    return Interval((center + lower_correction).lo, (center + upper_correction).hi)


@lru_cache(maxsize=None)
def log2_binom(n: int, k: int) -> Interval:
    if k < 0 or k > n:
        raise ValueError("binomial index outside range")
    k = min(k, n - k)
    if k == 0:
        return Interval.exact(0)
    natural = ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
    return natural / LN2


def self_check() -> None:
    """Check that small exact logarithms lie inside the outward intervals."""

    for n in range(2, 180, 7):
        for k in {1, n // 3, n // 2, n - 1}:
            interval = log2_binom(n, k)
            direct = log2_int(math.comb(n, k))
            if interval.hi < direct.lo or direct.hi < interval.lo:
                raise SystemExit("outward log2 binomial self-check failed")

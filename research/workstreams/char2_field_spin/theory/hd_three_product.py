"""Three-product two-state post-inner: exact local law, not a certificate.

Independently sample nonzero d1,d2,f at each position. Use
  a=s1+d1*x; y=a+f*s2; b=s2+d2*x; (s1,s2)=(a,b).
Thus c=d1 in the earlier four-product ensemble. The distribution changes.
"""
from fractions import Fraction as F


def support_law(old, active, q):
    if old not in range(4) or active not in (0, 1) or q < 4:
        raise ValueError("invalid support or field size")
    p = F(1, q-1)
    if not active:
        if old == 0:
            return {(0, 0): F(1)}
        if old == 3:
            return {(3, 0): p, (3, 1): 1-p}
        return {(old, 1): F(1)}
    result = {}
    for first in (0, 1):
        a = (p if first == 0 else 1-p) if old & 2 else F(first)
        for second in (0, 1):
            b = (p if second == 0 else 1-p) if old & 1 else F(second)
            if not a*b:
                continue
            # With old s2=0, y is exactly the new s1.
            zero = (p if first else F(0)) if old & 1 else F(1-first)
            for emitted, probability in ((0, zero), (1, 1-zero)):
                if probability:
                    result[2*first+second, emitted] = a*b*probability
    return result


def collapsed_majorant(v, q):
    """Entrywise zero/live envelope also valid for the four-product inner.

    Individual live support classes are each bounded by the same row;
    nonnegative continuation therefore preserves this majorization.
    """
    v = F(v)
    if not 0 < v <= 1 or q < 4:
        raise ValueError("need a positive output mark at most one")
    p = F(1, q-1)
    b = v+(1-v)*p
    return (((F(1), F(0)), (F(0), b)),
            ((F(0), v), (v*p*p, b)))


def exact_probability(n, weight, cutoff, q):
    """Average over all input supports, fixed arbitrary nonzero values."""
    from math import comb
    if not 0 <= weight <= n or not 0 <= cutoff <= n:
        raise ValueError("invalid finite class")
    states = {(0, 0, 0): F(1)}
    for _ in range(n):
        following = {}
        for (used, old, emitted), mass in states.items():
            for active in (0, 1):
                if used+active > weight:
                    continue
                for (new, out), probability in support_law(old, active, q).items():
                    if emitted+out > cutoff:
                        continue
                    key = used+active, new, emitted+out
                    following[key] = following.get(key, F(0))+mass*probability
        states = following
    return sum((mass for (used, _, _), mass in states.items() if used == weight), F(0))/comb(n, weight)


if __name__ == '__main__':
    for old in range(4):
        print(old, {str(key): str(value) for key, value in support_law(old, 1, 65536).items()})

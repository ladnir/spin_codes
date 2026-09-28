"""Verified polynomial bound for even binary codes of length n and d>=38.

SciPy proposes coefficients only. Every accepted bound is checked using exact
rational arithmetic at every permitted nonzero distance. No files are written.
"""

from fractions import Fraction as F
from math import comb, log2
from scipy.optimize import linprog


def krawtchouk(n, w):
    values = [1, n - 2*w]
    for j in range(1, n):
        numerator = (n-2*w)*values[-1] - (n-j+1)*values[-2]
        value, remainder = divmod(numerator, j+1)
        assert remainder == 0
        values.append(value)
    return values[:n+1]


def bound(n, d=38):
    if n < d:
        return F(1), ()
    allowed = list(range(d + d % 2, n + 1, 2))
    binomial = [comb(n, j) for j in range(n+1)]
    ratios = [[F(v, binomial[j]) for j, v in enumerate(krawtchouk(n, w))]
              for w in allowed]
    result = linprog([1.] * n,
                     A_ub=[[float(v) for v in row[1:]] for row in ratios],
                     b_ub=[-1.] * len(allowed), bounds=(0, None), method='highs')
    if not result.success:
        raise ArithmeticError(f'LP failed at n={n}: {result.message}')
    # x_j = coefficient_j * binomial(n,j); impose coefficient_0=1.
    x = [F(0)] + [F(float(v)).limit_denominator(10**12) for v in result.x]
    assert all(v >= 0 for v in x)
    values = [sum((v * r for v, r in zip(x, row)), F(0)) for row in ratios]
    if not all(v < 0 for v in values):
        raise ArithmeticError(f'No strictly negative polynomial tail at n={n}')
    scale = max([F(1)] + [-1/v for v in values])
    x = [v * scale for v in x]
    assert all(1 + sum((v*r for v,r in zip(x,row)), F(0)) <= 0 for row in ratios)
    # F(0) = 1+sum_j x_j; constant Krawtchouk coefficient is exactly one.
    upper = 1 + sum(x, F(0))
    return upper, tuple(x)


def support_minima(max_rank=16, d=38):
    minima = [None] * (max_rank + 1)
    remaining = set(range(1, max_rank + 1))
    for n in range(d, 257):
        upper, _ = bound(n, d)
        # At a length not excluded by this bound, fix a lower bound on the
        # minimum possible support. This is not an existence statement.
        for h in sorted(remaining):
            if upper >= 1 << h:
                minima[h] = n
                remaining.remove(h)
        if not remaining:
            return minima
    raise AssertionError('support scan exhausted')


def dimension_caps(length=256, dimension=128, d=38, last_lp=104):
    """Dimension upper bounds for shortened subcodes on each fixed support.

    Shortening on u-v coordinates loses at most u-v dimensions. Thus a
    dimension bound k_v at smaller length v gives k_u <= k_v + u-v.
    A failed numerical proposal contributes no bound.
    """
    caps = []
    accepted = failed = 0
    for n in range(length + 1):
        if n < d:
            caps.append(0)
            continue
        cap = min(dimension, n, caps[-1] + 1)
        # The Hamming bound: disjoint balls of radius floor((d-1)/2).
        volume = sum(comb(n,j) for j in range((d-1)//2+1))
        while (1 << cap) * volume > 1 << n:
            cap -= 1
        # Griesmer for linear codes, also retained when LP proposals fail.
        while sum((d+(1<<j)-1) >> j for j in range(cap)) > n:
            cap -= 1
        if n <= last_lp:
            try:
                upper, _ = bound(n,d)
            except ArithmeticError:
                failed += 1
            else:
                accepted += 1
                while 1 << cap > upper:
                    cap -= 1
        caps.append(cap)
    print(f'Shortening dimension caps: {accepted} exact LP witnesses accepted, {failed} proposals rejected', flush=True)
    return caps


def self_test():
    for n in range(1, 15):
        for w in range(n+1):
            values = krawtchouk(n,w)
            for j,v in enumerate(values):
                direct = sum((-1)**i * comb(w,i) * comb(n-w,j-i)
                             for i in range(max(0,j-n+w),min(w,j)+1))
                assert v == direct
    # Familiar tiny even codes: repetition and even-parity code exist.
    for n in (4,6,8):
        assert bound(n,n)[0] >= 2
        assert bound(n,2)[0] >= 1 << (n-1)
    # RM(1,3), the [8,4,4] extended Hamming code.
    assert bound(8,4)[0] >= 16
    print('Krawtchouk recurrence and small even-code checks passed', flush=True)


if __name__ == '__main__':
    self_test()
    for n in (76, 80, 87, 96, 104):
        upper, _ = bound(n)
        print(n, 'log2 exact size upper', log2(upper.numerator)-log2(upper.denominator), flush=True)
    print('Rank support lower bounds:', support_minima())
    caps = dimension_caps()
    print('Selected shortened dimensions:', [(u,caps[u]) for u in (80,87,91,96,104,128,160,192,256)])

"""Uniform permitted-pair term via a sharper positive-occupancy kernel.

This is one term of the positive source expansion, not a distance certificate.
The two-state upper control applies to both post-inner distributions.
"""
from fractions import Fraction as F
from math import comb

from hd_component_rank import arb_fraction, gap_moment_upper, power, run_kernel_upper


def first_return_sum(rows, base):
    """Sum over ALL h of the single-return coefficient in kernel F_h[1,1].

    The omitted common factor is v^2*p^2. Summing, not maximizing,
    gives a positive uniform upper bound on each occupancy.
    """
    result = base*0
    for gap in range(rows-1):
        coefficient = F(2*(rows+1)*(rows-gap-1), (gap+1)*(gap+2)*(gap+3))
        scalar = arb_fraction(coefficient) if type(base).__module__ == 'flint.types.arb' else coefficient
        result = result*base+scalar
    return result


def uniform_positive(rows, q, mark, minimum=1):
    """Directed 2x2 entrywise upper for every occupancy minimum..rows."""
    from flint import arb
    if not 0 < mark < 1 or q < 4:
        raise ValueError('need field size at least four and mark strictly between zero and one')
    mark = arb_fraction(mark)
    p = arb(1)/(q-1)
    base = mark+(1-mark)*p
    dense = base**rows
    lam = rows*p*p
    if not lam/3 < 1:
        raise ValueError('two-return tail needs a smaller parameter range')
    tail = lam*lam/(2*(1-lam/3))
    correction = mark*mark*p*p*first_return_sum(rows, base)+arb(tail.upper())
    eta = (1+p*p)**rows-1
    if not 1 <= minimum <= rows or rows < 2:
        raise ValueError('need at least two positions and a valid minimum')
    birth = mark*base**(rows-1)*gap_moment_upper(rows, minimum, 1, base)
    # One return spans the first through last active input. Coupling a
    # uniform h-subset to a uniform two-subset shows its range is larger.
    range_two = 2*((rows-1)-rows*base+dense)/(rows*(rows-1)*(1-base)**2)
    returned = mark*mark*p*p*range_two+p*p*((1+p*p)**(rows-1)-1)
    return [[returned, birth+eta],
            [p*p*(birth+eta), dense+correction]]


def pair_matrix(rows, q, mark):
    from flint import arb
    p, v = arb(1)/(q-1), arb_fraction(mark)
    base = v+(1-v)*p
    empty = [[arb(1), arb(0)], [arb(0), base**rows]]
    positive = uniform_positive(rows, q, mark)
    two_live = uniform_positive(rows, q, mark, 2)
    arbitrary = [[arb(1), positive[0][1]], positive[1][:]]
    full = run_kernel_upper(rows, rows, q, mark)
    result = [[arb(0) for _ in range(4)] for _ in range(4)]
    Q = q-1
    for i in range(2):
        for j in range(2):
            result[i][j] = empty[i][j]+Q*positive[i][j]
            result[i][2+j] = Q*positive[i][j]+Q*Q*two_live[i][j]
            result[2+i][j] = arbitrary[i][j]+Q*positive[i][j]
            result[2+i][2+j] = Q*positive[i][j]+Q*Q*full[i][j]
    return result


def receipt(denominator=1580, precision=128):
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, regions, q = 32768, 64, 65536
        mark = F(denominator-1, denominator)
        moment = sum(power(pair_matrix(rows, q, mark), regions)[0])-1
        cutoff = (2*rows*regions)//5
        per_profile = moment*arb(q-1)**-65*arb_fraction(mark)**-cutoff
        total = comb(rows, 2)*per_profile
        return dict(mark=str(mark), precision=precision,
                    per_profile_log2_upper=str(per_profile.log()/arb(2).log()),
                    all_pairs_log2_upper=str(total.log()/arb(2).log()),
                    all_pairs_below_2_to_minus_40=bool(total < arb(2)**-40),
                    any_fixed_regionwise_row_order=True,
                    whole_distance_certificate=False)
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    print(receipt())

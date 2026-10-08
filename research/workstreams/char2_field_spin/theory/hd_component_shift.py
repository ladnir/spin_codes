"""Changed ensemble: independent cyclic group offsets before the prefix.

Proof diagnostics only. No timing or full-distance certificate is implied.
"""
from fractions import Fraction as F
from math import comb

from hd_component_pair import uniform_positive
from hd_component_rank import arb_fraction, gap_moment_upper, power, run_kernel_upper
from hd_component_sparse import minimum_kernel


def birth_sum(rows, q, mark):
    """Sum over h=1..R of the no-return birth term, by a positive identity."""
    from flint import arb
    v, p = arb_fraction(mark), arb(1)/(q-1)
    base = v+(1-v)*p
    result = arb(0)
    for gap in range(rows):
        result = result*base+arb(rows+1)/((gap+1)*(gap+2))
    return v*result


def birth_return_sum(rows, q, mark):
    """Directed sum over h of exactly-one-return birth mass."""
    from flint import arb
    v, p = arb_fraction(mark), arb(1)/(q-1)
    base = v+(1-v)*p
    result = arb(0)
    for gap in range(rows-2):
        coefficient = arb(6*(rows+1)*(rows-gap-2))/((gap+2)*(gap+3)*(gap+4))
        result = result*base+coefficient
    return v**3*p*p*result


def shift_matrix(rows, q, mark, active, group=4):
    """Uniform mean matrix for all fixed cut sets under group-aligned shifts."""
    from flint import arb
    if rows % group or not 1 <= active <= rows:
        raise ValueError('group must divide rows, and cut count must be valid')
    v, p = arb_fraction(mark), arb(1)/(q-1)
    base = v+(1-v)*p
    lam = rows*p*p
    if not lam/3 < 1:
        raise ValueError('higher return tail outside allowed range')
    # For a birth with >=2 returns, bound all mark factors by 1 and
    # C(h-t-1,t) by h^t/t!. This controls each h, not the sum over h.
    tail = lam*lam/(2*(1-lam/3))
    single_return_sum = birth_return_sum(rows, q, mark)
    shared = uniform_positive(rows, q, mark)
    positive = [None]+[minimum_kernel(rows, q, mark, m, shared)
                       for m in range(1, active+1)]
    for m in range(1, active+1):
        no_return = v*base**(rows-1)*gap_moment_upper(rows, m, 1, base)
        positive[m][0][1] = no_return+single_return_sum+tail
        positive[m][1][0] = p*p*positive[m][0][1]
    full = run_kernel_upper(rows, rows, q, mark)
    empty = [[arb(1), arb(0)], [arb(0), base**rows]]
    # With unequal carried flags, occupancy varies injectively while a
    # fixed circular gap crosses the boundary. There are active such gaps;
    # sampling R/group offsets costs at most group*active/R per occupancy.
    multiplicity = arb(group*active)/rows
    sum_birth = birth_sum(rows, q, mark)
    average = []
    for minimum in range(active+1):
        m = max(1, minimum)
        if minimum > 1:
            previous_birth = v*base**(rows-1)*gap_moment_upper(rows, minimum-1, 1, base)
            sum_birth -= previous_birth
        birth_bound = multiplicity*(sum_birth+single_return_sum+rows*tail)
        # Both bounds are valid; select an upper endpoint, never a midpoint.
        birth_bound = arb(min(birth_bound.upper(), positive[m][0][1].upper()))
        zero_return = positive[m][0][0]
        if minimum == 0:
            zero_return += multiplicity
        average.append([[zero_return, birth_bound],
                        [p*p*birth_bound, positive[m][1][1]]])
    result = [[arb(0) for _ in range(4)] for _ in range(4)]
    for old in (0, 1):
        for new in (0, 1):
            for live in range(active):
                if old == new == 0 and live == 0:
                    kernel = empty
                elif old == new == 1 and live == active-1:
                    kernel = full
                elif old != new:
                    kernel = average[live+new]
                else:
                    kernel = positive[live+new]
                weight = comb(active-1, live)*(q-1)**(live+new)
                for i in range(2):
                    for j in range(2):
                        result[2*old+i][2*new+j] += weight*kernel[i][j]
    return result


def receipt(active=3, denominator=1200, precision=128, group=4):
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, regions, q = 32768, 64, 65536
        mark = F(denominator-1, denominator)
        moment = sum(power(shift_matrix(rows, q, mark, active, group), regions)[0])-1
        total = (comb(rows, active)*moment*arb(q-1)**(-32*active-1)
                 *arb_fraction(mark)**-838860)
        return dict(active=active, mark=str(mark), group=group,
                    log2_upper=str(total.log()/arb(2).log()),
                    below_2_to_minus_40=bool(total < arb(2)**-40),
                    changed_ensemble=True, whole_distance_certificate=False)
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--active', type=int, default=3)
    parser.add_argument('--denominators', type=int, nargs='+', default=[800, 1000, 1200, 1400])
    parser.add_argument('--group', type=int, default=4)
    args = parser.parse_args()
    for denominator in args.denominators:
        print(receipt(args.active, denominator, group=args.group), flush=True)

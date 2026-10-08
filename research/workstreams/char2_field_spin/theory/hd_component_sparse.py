"""Minimum-component bounds for sparse permitted-row classes.

These are positive upper bounds; a failed bound is not a code obstruction.
"""
from fractions import Fraction as F
from math import comb

from hd_component_pair import uniform_positive
from hd_component_rank import arb_fraction, gap_moment_upper, power, run_kernel_upper


def minimum_kernel(rows, q, mark, minimum, shared=None):
    """Uniform upper for h>=minimum, retaining the minimum range on returns."""
    from flint import arb
    p, v = arb(1)/(q-1), arb_fraction(mark)
    base = v+(1-v)*p
    if shared is None:
        shared = uniform_positive(rows, q, mark)
    if not 1 <= minimum <= rows:
        raise ValueError('positive minimum occupancy required')
    eta = (1+p*p)**rows-1
    birth = v*base**(rows-1)*gap_moment_upper(rows, minimum, 1, base)
    # One return spans the first and last live input. Coupling a uniform
    # minimum-subset inside an h-subset proves stochastic range domination.
    range_h = max(2, minimum)
    returned = (v*v*p*p*base**(rows-2)
                * gap_moment_upper(rows, range_h, 2, base)
                + p*p*((1+p*p)**(rows-1)-1))
    return [[returned, birth+eta], [p*p*(birth+eta), shared[1][1]]]


def class_matrix(rows, q, mark, active):
    """Uniform over every ordered set of active cuts in every region."""
    from flint import arb
    if not 1 <= active <= rows:
        raise ValueError('invalid number of permitted rows')
    v, p = arb_fraction(mark), arb(1)/(q-1)
    base = v+(1-v)*p
    empty = [[arb(1), arb(0)], [arb(0), base**rows]]
    full = run_kernel_upper(rows, rows, q, mark)
    shared = uniform_positive(rows, q, mark)
    positive = [None]+[minimum_kernel(rows, q, mark, m, shared)
                       for m in range(1, active+1)]
    arbitrary = [[arb(1), positive[1][0][1]], positive[1][1][:]]
    result = [[arb(0) for _ in range(4)] for _ in range(4)]
    for old in (0, 1):
        for new in (0, 1):
            for live in range(active):
                if old == new == 0 and live == 0:
                    kernel = empty
                elif old == new == 1 and live == active-1:
                    kernel = full
                elif live+new == 0:
                    kernel = arbitrary
                else:
                    kernel = positive[live+new]
                weight = comb(active-1, live)*(q-1)**(live+new)
                for i in range(2):
                    for j in range(2):
                        result[2*old+i][2*new+j] += weight*kernel[i][j]
    return result


def receipt(active=3, denominator=1000, precision=128):
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, regions, q = 32768, 64, 65536
        mark = F(denominator-1, denominator)
        moment = sum(power(class_matrix(rows, q, mark, active), regions)[0])-1
        per_profile = moment*arb(q-1)**(-32*active-1)*arb_fraction(mark)**-838860
        total = comb(rows, active)*per_profile
        return dict(active=active, mark=str(mark), precision=precision,
                    log2_upper=str(total.log()/arb(2).log()),
                    below_2_to_minus_40=bool(total < arb(2)**-40),
                    whole_distance_certificate=False)
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--active', type=int, default=3)
    parser.add_argument('--denominators', type=int, nargs='+', default=[600, 800, 1000, 1200, 1400])
    args = parser.parse_args()
    for denominator in args.denominators:
        print(receipt(args.active, denominator), flush=True)

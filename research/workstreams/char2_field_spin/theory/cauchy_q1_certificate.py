"""Directed ONE-active-row certificate for the Cauchy/diagonal-state inner.

This is not a whole-code certificate.  Every arithmetic operation in the
positive transfer uses the existing 80-digit outward-rounded Interval type.
Frozen rational marks are scaled from y=9997/10000 at K20.  No optimizer or
floating arithmetic is used by certify().
"""

import argparse
from fractions import Fraction
import json
import math

from subfield_certificate import Interval, ONE, ZERO, integer_margin


def rows_for(k_log):
    return ((2**k_log+22*64-1)//(22*64))*64


WITNESSES = {16: Fraction(93553, 94000),
             18: Fraction(373553, 374000),
             20: Fraction(9997, 10000)}


def interval(value):
    value = Fraction(value)
    return Interval.rational(value.numerator, value.denominator)


def positive_power_sum(base, count):
    """Return base**count and sum(base**j,j=0..count-1), without subtraction."""
    if count < 0:
        raise ValueError('negative geometric length')
    result_power, result_sum = ONE, ZERO
    block_power, block_sum = base, ONE
    while count:
        if count & 1:
            result_sum += result_power*block_sum
            result_power *= block_power
        count >>= 1
        if count:
            block_sum += block_power*block_sum
            block_power *= block_power
    return result_power, result_sum


def marginal(mark, weight):
    """Directed M_w for w=0 or1; fresh state/input diagonals exposed once."""
    if weight not in (0, 1) or not 0 < mark <= 1:
        raise ValueError('one-row marginal requires w in {0,1}, 0<y<=1')
    y, omitted = interval(mark), interval(1-mark)
    return sum((math.comb(weight, i)*math.comb(16-weight, z)
                *y**(16-i-z)*omitted**(i+z)/255**(i+z)
                for i in range(weight+1) for z in range(4)), ZERO)


def regional_matrices(rows, mark):
    if rows < 16 or rows % 16:
        raise ValueError('regional length must be a positive multiple of16')
    packets = rows//16
    y, m0, m1 = interval(mark), marginal(mark, 0), marginal(mark, 1)
    power, geometric = positive_power_sum(m0, packets)
    average = geometric/packets
    empty = [[ONE, ZERO], [ZERO, power]]
    # Every return from a one-input packet has all16 output coordinates
    # nonzero.  The return fixes all four nonzero state diagonals.
    active = [[ZERO, y*average],
              [y**16*average/255**4, m1*m0**(packets-1)]]
    return empty, active


def projective_support_count(weight, parity=10):
    if weight <= parity:
        return 0
    count = sum((-1)**j*math.comb(weight-1, j)*256**(weight-parity-1-j)
                for j in range(weight-parity))
    if count <= 0:
        raise ArithmeticError('invalid positive MDS support multiplicity')
    return count


def q1_expression(rows, mark, *, length=32, parity=10, threshold=None):
    """Outward interval around the proved Q1 upper expression.

    Optional small length/parity values support independent Fraction tests;
    the public certificate always uses length32 and parity10.
    """
    if not 0 <= parity < length or length > 256:
        raise ValueError('invalid MDS shape')
    if threshold is None:
        threshold = (length*rows+9)//10-1
    if not 0 <= threshold < length*rows:
        raise ValueError('invalid output threshold')
    empty, active = regional_matrices(rows, mark)
    current = [[ONE, ZERO]]+[[ZERO, ZERO] for _ in range(length)]
    for region in range(length):
        following = [[ZERO, ZERO] for _ in range(length+1)]
        for weight in range(region+1):
            for present, matrix in ((0, empty), (1, active)):
                for end in range(2):
                    following[weight+present][end] += sum(
                        (current[weight][start]*matrix[start][end] for start in range(2)), ZERO)
        current = following
    total = sum(((current[weight][0]+current[weight][1])
                 *projective_support_count(weight, parity)
                 for weight in range(parity+1, length+1)), ZERO)
    return rows*total/interval(mark)**threshold


def certify(k_log):
    if k_log not in WITNESSES:
        raise ValueError('only frozen K16/K18/K20 witnesses are certified')
    rows, mark = rows_for(k_log), WITNESSES[k_log]
    n = 32*rows
    threshold = (n+9)//10-1
    bound = q1_expression(rows, mark, threshold=threshold)
    if not 0 < bound.lo <= bound.hi < 1:
        raise ArithmeticError('invalid Q1 probability interval')
    margin = integer_margin(bound.hi)
    if margin < 40:
        raise ArithmeticError('Q1 margin below40bits')
    return dict(construction='cauchy4-fresh-nonzero-state-diagonal-every16',
                outer=[32,22], rate='11/16', active_outer_rows=1,
                k_log=k_log, k=22*rows, rows=rows, n=n,
                output_threshold=threshold, y=str(mark),
                return_output_weight=16, return_probability_upper='1/255^4',
                arithmetic='80-digit outward-rounded Decimal intervals',
                failure_upper=str(bound.hi), class_margin_bits_lower=margin,
                rigorous=True, whole_code=False,
                scope='exactly one nonzero outer row; all other source rows zero')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, nargs='+', default=[16,18,20])
    args = parser.parse_args()
    for k_log in args.k_log:
        print(json.dumps(certify(k_log)), flush=True)

"""Directed finite certificate for a complete intermediate input-density class.

Proves the return-pivot support cover for every input weight in the
specified interval. It is NOT a whole-code certificate: lower densities
remain open and the independent dense theorem covers densities >= .24.
Floating optimization only supplies rational witnesses. Every bound and
the sum over all covered weights use outward-rounded interval arithmetic.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
import math

import return_cover
from subfield_certificate import Interval, ONE, ZERO, integer_margin, interval_matrix_product


@dataclass(frozen=True)
class Witness:
    x: Fraction
    y: Fraction
    tau: Fraction
    tilt: Fraction


# Frozen rational witnesses. The search may propose replacements, but the
# ordinary replay never trusts an optimizer or its floating objective.
WITNESSES = tuple((Fraction(density), Witness(*(Fraction(value) for value in values)))
    for density, values in (
        ('19/100', ('50334/4505','219/8372','180234/4601','1852/8277')),
        ('1/5', ('121894/8559','187/9341','397747/8352','337/1782')),
        ('21/100', ('136620/7387','145/9668','172437/2912','1551/9830')),
        ('11/50', ('62052/2525','19/1734','352908/4667','1117/8622')),
        ('23/100', ('215494/6409','50/6443','542415/5437','342/3281')),
        ('6/25', ('433371/9068','19/3597','1027630/7493','627/7675'))))


def interval(value):
    value = Fraction(value)
    return Interval.rational(value.numerator, value.denominator)


def choose_witness(density):
    selected = return_cover.screen(density)['logarithmic_parameters']
    values = [Fraction(math.exp(value)).limit_denominator(10000) for value in selected[:4]]
    return Witness(*values)


def transfer(witness):
    x, y, tau = (interval(value) for value in (witness.x, witness.y, witness.tau))
    if not witness.x > 0 or not 0 < witness.y <= 1 or not 1 <= witness.tau <= 255:
        raise ValueError('invalid Chernoff/rank witness')
    xp, yp, gp = [x**i for i in range(17)], [y**i for i in range(17)], [(ONE/tau)**i for i in range(21)]
    result = [[ONE, ZERO], [ZERO, ZERO]]
    for w in range(1, 17):
        value = math.comb(16, w)*xp[w]*yp[w]
        result[0][1] += value
        if w >= 5:
            result[0][0] += value*gp[4]
    for union in range(13, 17):
        for only_input in range(union+1):
            for only_output in range(union-only_input+1):
                both = union-only_input-only_output
                wi, wo = only_input+both, only_output+both
                multiplicity = (math.comb(16, union)*math.comb(union, only_input)
                                *math.comb(union-only_input, only_output))
                value = multiplicity*xp[wi]*yp[wo]*gp[only_input]
                result[1][1] += value
                if wi and wo and (both or wo != 1):
                    result[1][0] += value*gp[min(both, 4)]
    return result


def terminal_power(matrix, exponent):
    result = [[ONE, ZERO], [ZERO, ONE]]
    while exponent:
        if exponent & 1:
            result = interval_matrix_product(result, matrix)
        exponent >>= 1
        if exponent:
            matrix = interval_matrix_product(matrix, matrix)
    return result[0][0]+result[0][1]


def row_upper(lower, upper, witness, length=32, parity=10):
    """Uniform bound over a density interval; no monotonicity assumption hidden."""
    beta = witness.tau*witness.tilt
    endpoint = upper if beta >= 1 else lower
    return interval(1-lower)**length + interval(witness.tau)**(-parity)*interval(1+(beta-1)*endpoint)**length


def certify(k_log, first=1900, last=2400, grid=10000, witnesses=None):
    if not 0 < first < last <= grid:
        raise ValueError('invalid density range')
    length, dimension, period = 32, 22, 64
    rows = ((2**k_log+dimension*period-1)//(dimension*period))*period
    n, threshold = length*rows, (length*rows+9)//10-1
    if witnesses is None:
        witnesses = WITNESSES
    prefactors = {}
    for density, witness in witnesses:
        prefactors[density] = (interval(rows+1)**length
            *terminal_power(transfer(witness), n//16)
            /(interval(witness.tau)*interval(witness.y)**threshold))
    total, intervals, least_margin = ZERO, 0, None
    for index in range(first, last):
        lo = (index*n+grid-1)//grid
        hi = ((index+1)*n+grid-1)//grid-1
        if lo > hi:
            continue
        lower, upper = Fraction(lo, n), Fraction(hi, n)
        density, witness = min(witnesses, key=lambda item: abs(item[0]-(lower+upper)/2))
        base = witness.x*witness.tilt
        exponent = lo if base >= 1 else hi
        bound = ((hi-lo+1)*prefactors[density]*row_upper(lower, upper, witness)**rows
                 /interval(base)**exponent)
        if not 0 < bound.hi < 1:
            raise ArithmeticError(f'failed cell {index}/{grid}: upper={bound.hi}')
        margin = integer_margin(bound.hi)
        least_margin = margin if least_margin is None else min(least_margin, margin)
        total += bound
        intervals += 1
    margin = integer_margin(total.hi)
    if margin < 40:
        raise ArithmeticError(f'class margin {margin} is below target')
    return dict(k_log=k_log, k=dimension*rows, n=n, rows=rows,
                minimum_input_weight=(first*n+grid-1)//grid,
                maximum_input_weight=(last*n+grid-1)//grid-1,
                density_min=str(Fraction(first, grid)), density_max_exclusive=str(Fraction(last, grid)),
                output_threshold=threshold, intervals=intervals,
                failure_upper=str(total.hi), class_margin_bits_lower=margin,
                minimum_cell_margin_bits_lower=least_margin,
                rigorous=True, whole_code=False,
                scope='all words in this inverse-input-density interval only')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, nargs='+', default=(16, 18, 20))
    parser.add_argument('--first', type=int, default=1900)
    parser.add_argument('--search', action='store_true', help='propose new rational witnesses before directed replay')
    args = parser.parse_args()
    witnesses = ([(Fraction(index, 100), choose_witness(index/100)) for index in range(19, 25)]
                 if args.search else WITNESSES)
    print(json.dumps(dict(rational_witnesses=[dict(density=str(density),
                        **{name: str(getattr(witness, name)) for name in ('x','y','tau','tilt')})
                     for density, witness in witnesses])), flush=True)
    for k_log in args.k_log:
        print(json.dumps(certify(k_log, args.first, witnesses=witnesses)), flush=True)

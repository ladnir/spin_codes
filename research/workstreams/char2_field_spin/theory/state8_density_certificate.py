"""Directed K20 middle-input certificate for symmetric width16/state8.

Fixed GRS[256,176], one regional routing stage, independent nonzero input
and state diagonals. Covers every W in [.03N,.24N), not the whole code.
The ordinary replay uses frozen rational witnesses and no optimizer.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb

from diagonal_refresh_rank import rank_potential_bound
from return_cover_certificate import Witness, interval, row_upper, terminal_power
from state8_joint_return import minimum_return_weight, split_joint_return_terms
from subfield_certificate import HIGH, LOW, Interval, ONE, ZERO, integer_margin


# Large intermediate positive factors are cancelled only after bounding.
# Change exponent limits, not precision or directed rounding modes.
LOW.Emax = HIGH.Emax = 1_000_000_000
LOW.Emin = HIGH.Emin = -1_000_000_000

WITNESSES = tuple((F(density), Witness(*(F(value) for value in values)))
    for density, values in (
        ('3/100', ('2846/9467', '6127/8442', '44970/6877', '37191/9643')),
        ('1/25', ('2595/6262', '2382/3631', '22919/3255', '7212/2569')),
        ('3/50', ('697/1114', '591/1061', '78151/9313', '10281/5926')),
        ('2/25', ('2889/3529', '2489/5086', '9729/980', '11533/9615')),
        ('1/10', ('5685/5206', '3438/8263', '56244/4555', '142/165')),
        ('3/25', ('13637/9327', '3418/9991', '81515/5323', '863/1351')),
        ('13/100', ('9383/5067', '2805/9907', '147269/8194', '1690/3121')),
        ('7/50', ('22408/9357', '2169/9890', '61305/3017', '3797/8099')),
        ('3/20', ('9587/2233', '618/6083', '126828/5657', '3644/8839')),
        ('9/50', ('150391/9400', '106/5507', '318567/8246', '2500/9999')),
        ('21/100', ('44185/788', '47/9800', '610070/6199', '1207/9929')),
        ('6/25', ('1628077/9637', '8/5339', '255', '238/4139'))))


def interval_minimum(left, right):
    return Interval(min(left.lo, right.lo), min(left.hi, right.hi))


def transfer(witness):
    """Positive two-state packet cover, all arithmetic rounded outward."""
    if not (witness.x > 0 and 0 < witness.y <= 1 and
            1 <= witness.tau <= 255 and witness.tilt > 0):
        raise ValueError('invalid rational witness')
    x, y, tau = (interval(value) for value in (witness.x, witness.y, witness.tau))
    xp, yp = [x**j for j in range(17)], [y**j for j in range(17)]
    zp, gp = [(ONE-y)**j for j in range(17)], [(ONE/tau)**j for j in range(25)]
    gamma = [interval(max(rank_potential_bound(r, z, witness.tau, 256, 8)
                          for r in range(1, 9))) for z in range(8)]
    result = [[ONE, ZERO], [ZERO, ZERO]]
    for w in range(17):
        multiplicity = comb(16, w)*xp[w]
        if w:
            result[0][1] += multiplicity*yp[w]
            if w > 8:
                result[0][0] += multiplicity*yp[w]*gp[8]
            flats = interval(sum(F(comb(8, k), 1)/witness.tau**k
                                 for k in range(min(w-1, 7)+1)))
            pure = flats*yp[minimum_return_weight(w, 16)]*gp[8]
            joint = ZERO
            for i, z, count, charge in split_joint_return_terms(w, 16):
                selected_flats = interval(sum(F(comb(8, k), 1)/witness.tau**k
                                             for k in range(min(w-1, 7-z)+1)))
                joint += selected_flats*count*yp[16-i-z]*zp[i+z]*gp[charge]
            result[1][0] += multiplicity*interval_minimum(pure, joint)
        for i in range(w+1):
            for z in range(min(7, 16-w)+1):
                result[1][1] += (multiplicity*comb(w, i)*comb(16-w, z)
                                  *yp[16-i-z]*zp[i+z]*gp[i]*gamma[z])
    return result


def weight_geometry(k_log=20):
    length, dimension, period = 256, 176, 64
    rows = ((2**k_log+dimension*period-1)//(dimension*period))*period
    n = rows*length
    return dict(k_log=k_log, length=length, dimension=dimension, rows=rows,
                k=dimension*rows, n=n, threshold=(n+9)//10-1)


def integer_span(n, first=F(3, 100), last=F(6, 25)):
    if not 0 < first < last <= 1:
        raise ValueError('invalid half-open density interval')
    lo = (first.numerator*n+first.denominator-1)//first.denominator
    hi = (last.numerator*n+last.denominator-1)//last.denominator-1
    return lo, hi


def cell_upper(lo, hi, n, rows, witness, prefactor):
    if not 0 < lo <= hi < n:
        raise ValueError('invalid integer weight cell')
    source = row_upper(F(lo, n), F(hi, n), witness, length=256, parity=80)
    base = witness.x*witness.tilt
    endpoint = lo if base >= 1 else hi
    return (hi-lo+1)*prefactor*source**rows/interval(base)**endpoint


def retarget_witness(anchor_density, witness, density):
    """Keep the frozen source odds; choose a rational tilt for this cell."""
    if not 0 < anchor_density < 1 or not 0 < density < 1:
        raise ValueError('invalid witness density')
    odds = witness.tilt*anchor_density/(1-anchor_density)
    return Witness(witness.x, witness.y, witness.tau, odds*(1-density)/density)


def refined_cells(first, last, initial_span, evaluate, maximum_upper):
    """Deterministic, gap-free integer bisection; evaluate returns(key,bound)."""
    if first > last or initial_span < 1:
        raise ValueError('invalid initial cell partition')
    stack = [(lo, min(last, lo+initial_span-1))
             for lo in range(first, last+1, initial_span)][::-1]
    while stack:
        lo, hi = stack.pop()
        key, bound = evaluate(lo, hi)
        if not bound.hi.is_finite() or bound.hi <= 0:
            raise ArithmeticError('nonpositive or nonfinite directed bound')
        if bound.hi <= maximum_upper:
            yield lo, hi, key, bound
        elif lo == hi:
            raise ArithmeticError(f'uncertified single input weight {lo}: {bound.hi}')
        else:
            middle = (lo+hi)//2
            stack.extend(((middle+1, hi), (lo, middle)))


def certify(k_log=20, first=F(3, 100), last=F(6, 25), initial_grid=1000,
            cell_bits=256, witnesses=WITNESSES):
    if k_log != 20:
        raise ValueError('this checkpoint certifies K20 only')
    if initial_grid < 1 or cell_bits < 40 or not witnesses:
        raise ValueError('invalid replay parameters')
    g = weight_geometry(k_log)
    first_weight, last_weight = integer_span(g['n'], first, last)
    prefactors = tuple(interval(g['rows']+1)**256
                       *terminal_power(transfer(w), g['n']//16)
                       /(interval(w.tau)*interval(w.y)**g['threshold'])
                       for _, w in witnesses)
    def evaluate(lo, hi):
        midpoint = F(lo+hi, 2*g['n'])
        candidates = sorted(range(len(witnesses)), key=lambda i: abs(witnesses[i][0]-midpoint))[:3]
        options = [(i, cell_upper(lo, hi, g['n'], g['rows'],
                                  retarget_witness(*witnesses[i], midpoint), prefactors[i]))
                   for i in candidates]
        return min(options, key=lambda item: item[1].hi)
    maximum = interval(F(1, 2**cell_bits)).lo
    span = (g['n']+initial_grid-1)//initial_grid
    total, next_weight, cells, largest = ZERO, first_weight, 0, ZERO.hi
    uses, smallest_width, largest_width = [0]*len(witnesses), None, 0
    digest = sha256()
    for lo, hi, index, bound in refined_cells(first_weight, last_weight, span, evaluate, maximum):
        if lo != next_weight or hi < lo:
            raise ArithmeticError('density cells have a gap, overlap, or reversed endpoints')
        next_weight = hi+1
        total += bound
        largest = max(largest, bound.hi)
        cells += 1
        uses[index] += 1
        width = hi-lo+1
        smallest_width = width if smallest_width is None else min(smallest_width, width)
        largest_width = max(largest_width, width)
        digest.update(f'{lo}:{hi}:{index}\n'.encode('ascii'))
    if next_weight != last_weight+1:
        raise ArithmeticError('density cells do not reach the required endpoint')
    margin = integer_margin(total.hi)
    if margin < 40:
        raise ArithmeticError(f'complete class margin below40: {margin}')
    return dict(construction='symmetric-Cauchy width16/state8 fresh diagonals',
                outer=[256, 176], width=16, state=8, q=256,
                rate='11/16', k_log=k_log, k=g['k'], n=g['n'],
                rows=g['rows'], output_threshold=g['threshold'],
                minimum_input_weight=first_weight, maximum_input_weight=last_weight,
                covered_weights=last_weight-first_weight+1, density_min=str(first),
                density_max_exclusive=str(last), checked_cells=cells,
                smallest_cell_weights=smallest_width, largest_cell_weights=largest_width,
                partition_sha256=digest.hexdigest(), witness_cell_counts=uses,
                failure_upper=str(total.hi), class_margin_bits_lower=margin,
                minimum_cell_margin_bits_lower=integer_margin(largest),
                precision_digits=80, rigorous=True, whole_code=False,
                scope='complete inverse-input-weight interval only')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial-grid', type=int, default=1000)
    parser.add_argument('--cell-bits', type=int, default=256)
    args = parser.parse_args()
    print(json.dumps(certify(initial_grid=args.initial_grid, cell_bits=args.cell_bits)), flush=True)

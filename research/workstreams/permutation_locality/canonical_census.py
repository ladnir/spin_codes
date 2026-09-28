"""Exact local feedback census for the four-row/four-column block route.

Uniform column permutations and independent within-column row permutations
make each orbit depend only on the sorted four column weights. This is a
local cancellation check, not a full SPIN distance certificate.
"""
from collections import Counter
from fractions import Fraction
from math import comb,factorial,prod
import group_moment


def census():
    _,columns,_ = group_moment.maps()
    totals,zero = Counter(),Counter()
    for window in range(8):
        syndrome = [0]*65536
        for mask in range(65536):
            if mask:
                syndrome[mask] = syndrome[mask&(mask-1)] ^ columns[16*window+(mask&-mask).bit_length()-1]
            shape = tuple(sorted(((mask>>(4*j))&15).bit_count() for j in range(4)))
            totals[shape] += 1
            zero[shape] += syndrome[mask] == 0
    assert len(totals) == comb(8,4) == 70
    for shape,total in totals.items():
        arrangements = factorial(4)//prod(factorial(n) for n in Counter(shape).values())
        assert total == 8*arrangements*prod(comb(4,a) for a in shape)
    bad = {shape:count for shape,count in zero.items() if count and any(shape)}
    assert bad == {(0,2,3,3):1,(1,2,2,3):1}
    return totals,zero


def main():
    totals,zero = census()
    bad = {shape:count for shape,count in zero.items() if count and any(shape)}
    for shape,count in bad.items():
        print('column weights',shape,'zero choices',count,'total',totals[shape],
              'probability',Fraction(count,totals[shape]))
    print('Other nonempty shapes have no zero feedback. No full-distance claim.')


if __name__ == '__main__':
    main()

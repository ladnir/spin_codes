"""Exact coverage checks for disjoint boxes of integer row-type counts."""
import itertools
import math


def lattice_count(lower,upper,total):
    if len(lower) != len(upper) or len(lower) < 2:
        raise ValueError('Invalid count-box dimensions')
    if any(type(x) is not int for x in [*lower,*upper,total]):
        raise ValueError('Count-box coordinates must be integers')
    if total < 0 or any(x < 0 for x in lower):
        raise ValueError('Negative type counts')
    if any(a > b for a,b in zip(lower,upper)):
        return 0
    remaining = total-sum(lower)
    if remaining < 0 or sum(upper) < total:
        return 0
    widths = [b-a+1 for a,b in zip(lower,upper)]
    dimensions = len(widths)
    count = 0
    for size in range(dimensions+1):
        for subset in itertools.combinations(widths,size):
            free = remaining-sum(subset)
            if free >= 0:
                count += (-1)**size*math.comb(free+dimensions-1,dimensions-1)
    assert count >= 0
    return count


def check(boxes,total,minimum,categories):
    if not 1 <= minimum <= total or not boxes:
        raise ValueError('Invalid dense domain')
    root_lower,root_upper = [0]*categories,[total]*categories
    root_upper[0] -= minimum
    expected = lattice_count(root_lower,root_upper,total)
    count = 0
    accepted = []
    for box in boxes:
        lo,hi = box['lower'],box['upper']
        if len(lo) != categories or len(hi) != categories:
            raise ValueError('Inconsistent box dimensions')
        if any(a < c or b > d for a,b,c,d in zip(lo,hi,root_lower,root_upper)):
            raise ValueError('Box outside the dense domain')
        size = lattice_count(lo,hi,total)
        if size == 0:
            raise ValueError('Empty dense box')
        for other_lo,other_hi in accepted:
            intersection_lo = [max(a,b) for a,b in zip(lo,other_lo)]
            intersection_hi = [min(a,b) for a,b in zip(hi,other_hi)]
            if lattice_count(intersection_lo,intersection_hi,total):
                raise ValueError('Dense boxes overlap on an integer type')
        count += size
        accepted.append((lo,hi))
    if count != expected:
        raise ValueError('Dense boxes leave an uncovered integer type')
    return count

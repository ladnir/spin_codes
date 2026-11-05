"""Exact constant-time feasibility for unsigned alternating-pair profiles.

An indexed binary cycle has length 2*M. Its profile (r,a,b) records its
weight and its two alternating counts of 11 edges. For 0<r<2*M the exact
domain is

    max(0,r-M) <= a,b <= floor(r/2),  not (2*a == 2*b == r).

The two constant cycles have profiles (0,0,0) and (2*M,M,M).

Necessity follows from the two matching histograms. Being full under both
matchings would make a nonconstant cycle constant. For sufficiency when
0<r<M, the lane counts (ceil(r/2),floor(r/2)) satisfy the exact witness
inequalities in hd_shifted_general_count.feasible_witness. The same holds
at r=M unless a=b=0; that profile occupies one complete lane. Complement
the r<M construction for r>M. Thus no enumeration of lane counts is needed.

Box helpers intersect integer ranges with this domain. They count neither
supports nor probabilities and make no distance claim.
"""


def _validate_pairs(pairs):
    if not isinstance(pairs, int) or pairs < 1:
        raise ValueError('pairs must be a positive integer')


def unsigned_witness(pairs, weight, even_full, odd_full):
    """Return feasible lane counts (E,O), or None for an empty profile."""
    _validate_pairs(pairs)
    if (any(not isinstance(n, int) for n in (weight, even_full, odd_full))
            or not 0 <= weight <= 2*pairs):
        return None
    lower, upper = max(0, weight-pairs), weight//2
    if not lower <= even_full <= upper or not lower <= odd_full <= upper:
        return None
    if 0 < weight < 2*pairs and 2*even_full == 2*odd_full == weight:
        return None
    if weight == pairs and even_full == odd_full == 0:
        return pairs, 0
    return (weight+1)//2, weight//2


def unsigned_is_feasible(pairs, weight, even_full, odd_full):
    """Whether an indexed binary cycle realizes the unsigned profile."""
    return unsigned_witness(pairs, weight, even_full, odd_full) is not None


def _ranges(ranges):
    ranges = tuple(tuple(interval) for interval in ranges)
    if (len(ranges) != 3 or any(len(interval) != 2
            or any(not isinstance(n, int) for n in interval)
            or interval[0] > interval[1] for interval in ranges)):
        raise ValueError('Expected ordered inclusive integer ranges (r,a,b)')
    return ranges


def unsigned_box_witness(pairs, ranges, *, exclude_zero=False):
    """Return a feasible integer (r,a,b) in the box, or None.

    Ranges may extend outside physical counts. Intersect them with the
    domain first. Every weight in the resulting interval has valid
    matching histograms. If its first point is forbidden, the next weight
    is odd and therefore cannot be full under both matchings.
    """
    _validate_pairs(pairs)
    (rlo, rhi), (alo, ahi), (blo, bhi) = _ranges(ranges)
    alo, ahi = max(0, alo), min(pairs, ahi)
    blo, bhi = max(0, blo), min(pairs, bhi)
    if alo > ahi or blo > bhi:
        return None
    rlo = max(rlo, int(exclude_zero), 2*alo, 2*blo)
    rhi = min(rhi, 2*pairs, pairs+ahi, pairs+bhi)
    if rlo > rhi:
        return None
    a, b = max(alo, rlo-pairs), max(blo, rlo-pairs)
    if 0 < rlo < 2*pairs and 2*a == 2*b == rlo:
        if rlo == rhi:
            return None
        rlo += 1
        a, b = max(alo, rlo-pairs), max(blo, rlo-pairs)
    return rlo, a, b


def unsigned_box_is_empty(pairs, ranges, *, exclude_zero=False):
    """Exact emptiness, using only a fixed number of integer operations."""
    return unsigned_box_witness(pairs, ranges, exclude_zero=exclude_zero) is None


def unsigned_box_all_feasible(pairs, ranges, *, exclude_zero=False):
    """Whether every integer triple in the supplied box is feasible."""
    _validate_pairs(pairs)
    (rlo, rhi), (alo, ahi), (blo, bhi) = _ranges(ranges)
    if (not int(exclude_zero) <= rlo <= rhi <= 2*pairs
            or not 0 <= alo <= ahi <= pairs or not 0 <= blo <= bhi <= pairs
            or rlo < 2*max(ahi, bhi) or rhi > pairs+min(alo, blo)):
        return False
    # The only excluded interior profiles lie on (r,a,b)=(2*k,k,k).
    forbidden_lo = max(1, alo, blo, (rlo+1)//2)
    forbidden_hi = min(pairs-1, ahi, bhi, rhi//2)
    return forbidden_lo > forbidden_hi

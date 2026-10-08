"""Joint component-length control for sparse permitted-row classes.

The positive birth kernel is discretely convex in occupancy. Keeping it
inside the full component transfer reduces length choices to simplex
corners, rather than independently maximizing each matrix entry. A failed
class bound is not a code obstruction or a whole-distance certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb, factorial

from hd_component_pair import first_return_sum
from hd_component_rank import (arb_fraction, gap_moment_exact, gap_moment_upper,
                               run_kernel_exact, run_kernel_upper, component_moment)


def birth_exact(rows, h, q, mark):
    """No-return zero-to-live term, in a form independent of gap moments."""
    v, p = F(mark), F(1, q-1)
    base = v+(1-v)*p
    return v*sum(comb(t, h-1)*base**t for t in range(h-1, rows))/comb(rows, h)


def birth_mixture(rows, h, q, mark):
    """Positive mixture proving discrete convexity for h=1,...,rows."""
    v, p = F(mark), F(1, q-1)
    base = v+(1-v)*p
    return v*base**(rows-1)+v*(1-base)*sum(
        base**(t-1)*F(comb(t, h) if t >= h else 0, comb(rows, h))
        for t in range(1, rows))


def length_corners(rows, active, zero_first):
    """Vertices of one face of the integer component-length simplex."""
    if not 1 <= active <= rows:
        raise ValueError("invalid number of permitted cuts")
    base = (0,)+(1,)*active if zero_first else (1,)*(active+1)
    excess = rows-sum(base)
    if excess < 0:
        return ()
    indices = range(1, active+1) if zero_first else range(active+1)
    return tuple(dict.fromkeys(tuple(x+(excess if j == i else 0)
                                    for j, x in enumerate(base)) for i in indices))


def control_factory(rows, q, mark, directed=False):
    """Return the shared positive control indexed by (occupancy, minimum).

    Occupancy zero and full remain exact (or directed run-sum uppers).
    Other occupancies keep their exact birth term, minimum-span return
    bound, and the pair helper's uniform live-to-live correction.
    """
    if directed:
        from flint import arb
        scalar = arb_fraction
        gap = gap_moment_upper
        full = run_kernel_upper(rows, rows, q, mark)
    else:
        scalar = F
        gap = gap_moment_exact
        full = run_kernel_exact(rows, rows, q, mark)
    v, p = scalar(mark), scalar(F(1, q-1))
    zero, one = scalar(0), scalar(1)
    base = v+(1-v)*p
    dense, lam = base**rows, rows*p*p
    if not lam/3 < 1:
        raise ValueError("two-return tail requires R/(q-1)^2 < 3")
    eta = (1+p*p)**rows-1
    live_live = (dense+v*v*p*p*first_return_sum(rows, base)
                 + lam*lam/(2*(1-lam/3)))
    @lru_cache(None)
    def control(h, minimum=1):
        if h == 0:
            return [[one, zero], [zero, dense]]
        if h == rows:
            return full
        if not 1 <= minimum <= h < rows:
            raise ValueError("invalid positive occupancy/minimum")
        birth = v*base**(rows-1)*gap(rows, h, 1, base)
        span = max(2, minimum)
        returned = (v*v*p*p*base**(rows-2)*gap(rows, span, 2, base)
                    + p*p*((1+p*p)**(rows-1)-1))
        return [[returned, birth+eta], [p*p*(birth+eta), live_live]]
    return control


def run_kernel_interval(rows, low, high, q, mark, returns=3, directed=False):
    """Uniform run-sum upper for every occupancy low..high.

    For fixed selected-gap count, the gap moment decreases with occupancy.
    Each retained run coefficient increases with occupancy. Their separate
    endpoints therefore bound the whole term without an unconditional eta.
    """
    if not 0 <= low <= high <= rows:
        raise ValueError("invalid occupancy interval")
    if directed:
        from flint import arb
        scalar, gap = arb_fraction, gap_moment_upper
    else:
        scalar, gap = F, gap_moment_exact
    v, p = scalar(mark), scalar(F(1, q-1))
    zero, one = scalar(0), scalar(1)
    base = v+(1-v)*p
    if high == 0:
        return [[one, zero], [zero, base**rows]]
    positive_low = max(1, low)
    result = [[zero, zero], [zero, zero]]
    for t in range(min((high-1)//2, returns)+1):
        start = max(positive_low, 2*t+1)
        result[0][1] += (comb(high-t-1, t)*v**(2*t+1)*p**(2*t)
                         *base**(rows-2*t-1)*gap(rows, start, t+1, base))
    for t in range(min(high//2, returns)+1):
        start = max(positive_low, 2*t)
        result[1][1] += (comb(high-t, t)*v**(2*t)*p**(2*t)
                         *base**(rows-2*t)*gap(rows, start, t, base))
    for t in range(1, min(high//2, returns+1)+1):
        start = max(positive_low, 2*t)
        result[0][0] += (comb(high-t-1, t-1)*v**(2*t)*p**(2*t)
                         *base**(rows-2*t)*gap(rows, start, t+1, base))
    lam, next_t = high*p*p, returns+1
    if not lam/(next_t+1) < 1:
        raise ValueError("increase retained returns for this interval")
    tail = lam**next_t/factorial(next_t)/(1-lam/(next_t+1))
    result[0][1] += tail
    result[1][1] += tail
    result[0][0] += p*p*tail
    result[1][0] = p*p*result[0][1]
    if low == 0:
        # F0 has only diagonal entries; the live diagonal is already bounded.
        result[0][0] = arb(max(result[0][0].upper(), 1)) if directed else max(result[0][0], 1)
    return result


def grouped_triple_matrix(group_low, group_high, rows, q, mark, kernel):
    """Entrywise upper for the four lane triples within one positioned group.

    The physical group index lies in the inclusive supplied interval.
    Birth/return ranges retain complementary exterior occupancies; only
    positions within this single bin, not arbitrary positions, are merged.
    """
    from itertools import combinations
    result = None
    for lanes in combinations(range(4), 3):
        small = (4*group_low+lanes[0], lanes[1]-lanes[0], lanes[2]-lanes[1],
                 rows-4*group_high-lanes[2])
        large = (4*group_high+lanes[0], small[1], small[2],
                 rows-4*group_low-lanes[2])
        zero = kernel(0, 0)[0][0]*0
        matrix = [[zero for _ in range(4)] for _ in range(4)]
        for old, new in product((0, 1), repeat=2):
            for inside in product((0, 1), repeat=2):
                flags = (old, *inside, new)
                if old == new:
                    # The exterior sum is fixed: do not add independent
                    # endpoints of two complementary length intervals.
                    low = high = old*(rows-lanes[2]+lanes[0])+sum(
                        flag*length for flag, length in zip(inside, small[1:3]))
                else:
                    low = sum(flag*length for flag, length in zip(flags, small))
                    high = sum(flag*length for flag, length in zip(flags, large))
                block = kernel(low, high)
                weight = (q-1)**(sum(inside)+new)
                for i in range(2):
                    for j in range(2):
                        matrix[2*old+i][2*new+j] += weight*block[i][j]
        if result is None:
            result = matrix
        else:
            from flint import arb
            result = [[arb(max(result[i][j].upper(), matrix[i][j].upper()))
                       for j in range(4)] for i in range(4)]
    return result


def group_bins(groups, subdivision=4):
    """Fine near both scan ends, geometric in the interior; singleton zero."""
    if subdivision < 1:
        raise ValueError("positive subdivision required")
    edges = {0, groups, *range(1, min(subdivision, groups))}
    step = 1
    while step < groups:
        for multiple in range(subdivision, 2*subdivision):
            edge = multiple*step
            if edge < groups:
                edges.add(edge)
                edges.add(groups-edge)
        step *= 2
    ordered = sorted(edges)
    return tuple((left, right-1) for left, right in zip(ordered, ordered[1:]))


def grouped_triple_receipt(denominator=1000, precision=128, subdivision=4):
    """Directed sum for all four lane triples inside every L4 source group."""
    from flint import arb, ctx
    from hd_two_state_regional_strides import STRIDES
    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, q, regions = 32768, 65536, 64
        groups, mark = rows//4, F(denominator-1, denominator)
        @lru_cache(None)
        def kernel(low, high):
            return run_kernel_interval(rows, low, high, q, mark, directed=True)
        bins = group_bins(groups, subdivision)
        matrices = [grouped_triple_matrix(low, high, rows, q, mark, kernel)
                    for low, high in bins]
        group_bin = [0]*groups
        for index, (low, high) in enumerate(bins):
            group_bin[low:high+1] = [index]*(high-low+1)
        inverses = [pow(stride, -1, groups) for stride in STRIDES]
        # Group zero is fixed by every odd multiplier. Preserve each lane
        # triple there: merging the j1=0 and j1=1 faces would invent cycles.
        zero_moments = [component_moment([lanes]*regions, rows, q,
                                         lambda h: kernel(h, h))
                        for lanes in combinations(range(4), 3)]
        aggregate = sum(zero_moments)
        largest = arb(0)
        worst_group = None
        for group in range(1, groups):
            vector = [arb(1), arb(0), arb(0), arb(0)]
            for inverse in inverses:
                matrix = matrices[group_bin[inverse*group % groups]]
                vector = [sum(vector[i]*matrix[i][j] for i in range(4)) for j in range(4)]
            moment = sum(vector)-1
            aggregate += 4*moment
            if moment.upper() > largest.upper():
                largest, worst_group = moment, group
        scale = arb(q-1)**-97*arb_fraction(mark)**-838860
        total = aggregate*scale
        return dict(mark=str(mark), precision=precision, bins=len(bins), subdivision=subdivision,
                    occupancy_intervals=kernel.cache_info().currsize,
                    covered_permitted_triples=4*groups,
                    log2_upper=str(total.log()/arb(2).log()),
                    worst_nonzero_group=worst_group,
                    worst_nonzero_profile_log2_upper=str((largest*scale).log()/arb(2).log()),
                    group_zero_profiles=[dict(cuts=lanes,
                        log2_upper=str((moment*scale).log()/arb(2).log()))
                        for lanes, moment in zip(combinations(range(4), 3), zero_moments)],
                    below_2_to_minus_40=bool(total < arb(2)**-40),
                    all_permitted_triples=False, whole_distance_certificate=False)
    finally:
        ctx.prec = previous


def geometry_matrix(lengths, q, control):
    """Full four-state component/post transfer with joint lengths retained."""
    active = len(lengths)-1
    if lengths[0] < 0 or any(x <= 0 for x in lengths[1:]):
        raise ValueError("only the first component may have zero length")
    zero = control(0)[0][0]*0
    result = [[zero for _ in range(4)] for _ in range(4)]
    for old, new in product((0, 1), repeat=2):
        for interior in product((0, 1), repeat=active-1):
            flags = (old, *interior, new)
            h = sum(flag*length for flag, length in zip(flags, lengths))
            minimum = sum(flag for flag, length in zip(flags, lengths) if length > 0)
            kernel = control(h, minimum)
            weight = (q-1)**(sum(interior)+new)
            for i in range(2):
                for j in range(2):
                    result[2*old+i][2*new+j] += weight*kernel[i][j]
    return result


def corner_matrices(rows, active, q, mark, directed=False):
    control = control_factory(rows, q, mark, directed)
    corners = (*length_corners(rows, active, True),
               *length_corners(rows, active, False))
    return [(lengths, geometry_matrix(lengths, q, control)) for lengths in corners]


def rowwise_upper(matrices, regions, directed=False):
    """Backward positive upper; each row retains its joint vertex choice."""
    value = [matrices[0][0][0]*0+1]*4
    for _ in range(regions):
        candidates = [[sum(matrix[i][j]*value[j] for j in range(4))
                       for i in range(4)] for matrix in matrices]
        if directed:
            from flint import arb
            value = [arb(max(candidate[i].upper() for candidate in candidates))
                     for i in range(4)]
        else:
            value = [max(candidate[i] for candidate in candidates) for i in range(4)]
    return value[0]-1


def receipt(active=3, denominator=1000, precision=128):
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, q, regions = 32768, 65536, 64
        mark = F(denominator-1, denominator)
        corners = corner_matrices(rows, active, q, mark, directed=True)
        moment = rowwise_upper([matrix for _, matrix in corners], regions, directed=True)
        total = (comb(rows, active)*moment*arb(q-1)**(-32*active-1)
                 * arb_fraction(mark)**-838860)
        return dict(active=active, mark=str(mark), precision=precision,
                    corners=len(corners), log2_upper=str(total.log()/arb(2).log()),
                    below_2_to_minus_40=bool(total < arb(2)**-40),
                    any_fixed_regionwise_row_order=True,
                    whole_distance_certificate=False)
    finally:
        ctx.prec = previous


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--active", type=int, default=3)
    parser.add_argument("--denominators", type=int, nargs="+", default=[800, 1000, 1200])
    parser.add_argument("--grouped-triples", action="store_true")
    parser.add_argument("--bin-subdivision", type=int, default=4)
    args = parser.parse_args()
    for denominator in args.denominators:
        print(grouped_triple_receipt(denominator, subdivision=args.bin_subdivision) if args.grouped_triples
              else receipt(args.active, denominator), flush=True)

"""All-valid run-aligned boxes for the alternating-pair support domain.

Coordinates are (r,c,d), where c=r-a-b and d=a-b. The inverse is
a=(r-c+d)/2, b=(r-c-d)/2. Integer profiles require r-c+d even.
This helper deliberately rejects boxes crossing a histogram boundary.
The nonconstant domain c>=1 excludes both constant supports.
"""

from fractions import Fraction as F
from itertools import product
from math import prod


def validate_run_box(ranges, pairs=4096):
    ranges = tuple(tuple(row) for row in ranges)
    if (not isinstance(pairs, int) or pairs < 1 or len(ranges) != 3
            or any(len(row) != 2 or any(not isinstance(n, int) for n in row)
                   or row[0] > row[1] for row in ranges)
            or ranges[1][0] < 1):
        raise ValueError('Three integer intervals, positive pairs, and c>=1 required')
    for r, c, d in product(*(tuple(dict.fromkeys(row)) for row in ranges)):
        # Twice the histogram entries; no rounding or parity assumption is
        # appropriate for continuous vertices of the convex envelope.
        if min(2*pairs-r-c+d, 2*(c-d), r-c+d,
               2*pairs-r-c-d, 2*(c+d), r-c-d) < 0:
            raise ValueError('Every run-box corner must have nonnegative histograms')
    return ranges


def run_box_vertices(ranges, pairs=4096):
    """Return at most eight exact (r,a,b) vertices, retaining half-integers."""
    ranges = validate_run_box(ranges, pairs)
    return tuple((F(r), F(r-c+d, 2), F(r-c-d, 2))
                 for r, c, d in product(*(tuple(dict.fromkeys(row)) for row in ranges)))


def run_box_profile_count(ranges, pairs=4096):
    """Exact number of feasible integer (r,a,b) profiles, not supports."""
    ranges = validate_run_box(ranges, pairs)
    lengths = tuple(hi-lo+1 for lo, hi in ranges)
    # The character (-1)^(r-c+d) factors over the three intervals. For an
    # odd-length interval its character sum is the parity of its low end.
    differences = tuple((1 if lo % 2 == 0 else -1) if length % 2 else 0
                        for (lo, _), length in zip(ranges, lengths))
    return (prod(lengths)+prod(differences))//2


def run_box_support_count_upper(ranges, x, y, pairs=4096):
    """Bound all indexed supports at once; do not multiply by profile count."""
    from hd_shifted_run_count import transfer_value

    ranges = validate_run_box(ranges, pairs)
    value = transfer_value(pairs, x, y)  # also validates exact positive markers
    if not run_box_profile_count(ranges, pairs):
        return 0
    x, y = F(x), F(y)
    minimum = x**ranges[0][0 if x >= 1 else 1]*y**ranges[1][0 if y >= 1 else 1]
    upper = value/minimum
    return upper.numerator//upper.denominator


def run_box_receipt(witness, ranges, x, y):
    """Directed run-box envelope using an existing CachedBound instance.

    Only its regional numerators and count-free convex scalar are reused.
    The run transfer bounds the total support count. The parity-aware
    profile count is ledger metadata, not an extra probability factor.
    """
    from flint import arb, ctx
    from hd_shifted_alias_certificate import rational
    from hd_shifted_coverage import rational_ceiling
    from hd_shifted_run_count import transfer_value
    from hd_shifted_unoriented_certificate import matrix_product
    from hd_sumrank_paired_certificate import scalar_power

    ranges = validate_run_box(ranges)
    profiles = run_box_profile_count(ranges)
    if not profiles:
        raise ValueError('Run box contains no integer profile')
    vertices = run_box_vertices(ranges)
    supports = run_box_support_count_upper(ranges, x, y)
    with ctx.workprec(witness.precision):
        if witness.grid is None or witness.uniform_aliases:
            regions = witness.regions
        else:
            regions = []
            for parity, atlas in enumerate(witness.regions):
                singles = [r-2*p[parity+1] for p in vertices for r in (p[0],)]
                full = [p[parity+1] for p in vertices]
                floor = lambda v: v.numerator//v.denominator
                ceil = lambda v: -((-v.numerator)//v.denominator)
                targets = ((floor(min(singles)), ceil(max(singles))),
                           (floor(min(full)), ceil(max(full))))
                regions.append(atlas.numerator(targets))
        run = scalar_power(matrix_product(*regions), 128)
        operator = (run[0]+run[1]).log()/arb(2).log()
        scalar = max(witness.scalar(point, count_matching=None).upper() for point in vertices)
        log_count = arb(supports).log()/arb(2).log()
        separate = operator+scalar+log_count
        log_x = rational(F(x)).log()/arb(2).log()
        log_y = rational(F(y)).log()/arb(2).log()
        log_transfer = rational(transfer_value(4096, x, y)).log()/arb(2).log()
        # Sum_U 2^scalar(U) <= G(x,y) * max_box 2^(scalar-r log x-c log y).
        # The run-marker term is affine, so the same eight vertices suffice.
        tilted_scalar = max((witness.scalar(point, count_matching=None)
                             -rational(point[0])*log_x
                             -rational(point[0]-point[1]-point[2])*log_y).upper()
                            for point in vertices)
        coupled = operator+log_transfer+tilted_scalar
        bound = min(separate.upper(), coupled.upper())
        ceiling = rational_ceiling(bound)
        separate_ceiling = rational_ceiling(separate)
        coupled_ceiling = rational_ceiling(coupled)
    return dict(profile_coordinates=['r', 'c', 'd'], ranges=ranges,
                integer_profile_triples=profiles, convex_vertices=len(vertices),
                count_method='coupled_run_transfer_or_separate', count_markers=[str(F(x)), str(F(y))],
                log2_support_count_enclosure=str(log_count), log2_operator_enclosure=str(operator),
                log2_scalar_enclosure=str(scalar), log2_upper_rational=str(ceiling),
                separate_upper_rational=str(separate_ceiling),
                coupled_upper_rational=str(coupled_ceiling),
                log2_run_transfer_enclosure=str(log_transfer),
                log2_tilted_scalar_enclosure=str(tilted_scalar),
                log2_upper_enclosure=str(bound), directed=True,
                covers_all_even_odd_splits=True, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    import json
    from hd_shifted_cover_backend import CachedBound

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r-radius', type=int, default=32)
    parser.add_argument('--c-radius', type=int, default=4)
    parser.add_argument('--d-radius', type=int, default=4)
    parser.add_argument('--grid', type=int, default=128)
    parser.add_argument('--precision', type=int, default=256)
    args = parser.parse_args()
    marker = (F(26675286813512279, 301536634975809212), F(179919979, 892083004))
    witness = CachedBound(F(2172973163350011, 4503599627370496),
                          F(3695833225349633, 36028797018963968), (marker, marker),
                          grid=args.grid, precision=args.precision)
    box = ((1024-args.r_radius, 1024+args.r_radius),
           (224-args.c_radius, 224+args.c_radius), (-args.d_radius, args.d_radius))
    print(json.dumps(dict(witness=witness.witness(), receipt=run_box_receipt(witness, box, F(25,31), F(7,775)))))

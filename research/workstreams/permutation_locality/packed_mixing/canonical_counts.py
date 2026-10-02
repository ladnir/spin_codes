"""Shortened-code counting for canonical packed blocks and exact CDF transport.

No claim of a whole-code certificate. The CLI authenticates the retained BCH
premises, builds exactly checked shortened-dimension bounds, and reports only
the resulting one-group expected support CDF. Numeric log2 displays are
diagnostics; stored caps and convolutions are exact integers/rationals.
"""
import argparse
from fractions import Fraction as Q
from math import comb, log2, prod
from pathlib import Path
import sys

from local_models import convolve, full_block, integer, checked_pmf


def rank_total(dimension, rows, rank):
    if rank > min(dimension, rows):
        return 0
    numerator = prod(((1 << dimension)-(1 << i))*((1 << rows)-(1 << i)) for i in range(rank))
    denominator = prod((1 << rank)-(1 << i) for i in range(rank))
    result, remainder = divmod(numerator, denominator)
    if remainder:
        raise ArithmeticError('rank count must be integral')
    return result


def canonical_rank_caps(dimensions, width=8, rows=4, original_rank_cdfs=None):
    """Count tuples contained in at most H specified-width canonical blocks.

    dimensions[u] must bound the shortened-code dimension on EVERY u-subset.
    Every tuple with <=H active blocks lies in some union of H blocks. Counting
    all such supersets overcounts, so C(groups,H)*rank_total(d(width*H),rows,r)
    is a valid rank-r cumulative cap. Optional original CDFs are authenticated
    externally; only valid caps may be provided.
    """
    dimensions = tuple(dimensions)
    length = len(dimensions)-1
    integer(width, 1, length, 'width')
    integer(rows, 1, 16, 'rows')
    if length % width or dimensions[0] != 0:
        raise ValueError('valid dimensions and integral canonical block count required')
    for u, dimension in enumerate(dimensions):
        integer(dimension, 0, u, 'shortened dimension')
    if any(a > b for a, b in zip(dimensions, dimensions[1:])):
        raise ValueError('shortened-dimension caps must be nondecreasing')
    groups, dimension = length//width, dimensions[-1]
    if original_rank_cdfs is not None:
        if len(original_rank_cdfs) != rows or any(len(row) != length+1 for row in original_rank_cdfs):
            raise ValueError('one full original CDF per rank required')
        for row in original_rank_cdfs:
            if row[0] != 0 or any(type(x) is not int or x < 0 for x in row) or any(a > b for a, b in zip(row, row[1:])):
                raise ValueError('nonnegative integer cumulative rank caps required')
    result = []
    for r in range(1, rows+1):
        total = rank_total(dimension, rows, r)
        row = []
        for h in range(groups+1):
            bound = min(total, comb(groups, h)*rank_total(dimensions[width*h], rows, r))
            if original_rank_cdfs is not None:
                bound = min(bound, original_rank_cdfs[r-1][width*h])
            row.append(bound)
        for h in range(groups-1, -1, -1):
            row[h] = min(row[h], row[h+1])
        result.append(tuple(row))
    return tuple(result)


def total_caps(rank_caps):
    return tuple(map(sum, zip(*rank_caps)))


def transport_cdf(caps, local):
    """Expected output-support CDF after independent rank-free local maps.

    local is the support PMF of a nonempty block and must assign no zero mass.
    For each output cutoff, the kernel Pr[sum of H iid local supports <=cutoff]
    decreases with H. Abel summation therefore permits input CDF differences
    INSIDE THIS TRANSFORM only. The returned object is another CDF cap, not
    pointwise shell domination; its differences are not shell bounds either.
    """
    caps = tuple(map(Q, caps))
    local = checked_pmf(local)
    if len(caps) < 2 or caps[0] != 0 or local[0] != 0:
        raise ValueError('nonzero-message CDF and nonzero local support required')
    if any(a < 0 or a > b for a, b in zip(caps, caps[1:])) or caps[-1] < 0:
        raise ValueError('nonnegative nondecreasing cumulative caps required')
    size = (len(caps)-1)*(len(local)-1)
    comparison = [Q(0)]*(size+1)
    law = (Q(1),)
    for h in range(1, len(caps)):
        law = convolve(law, local)
        difference = caps[h]-caps[h-1]
        if difference:
            for w, p in enumerate(law):
                comparison[w] += difference*p
    result, running = [], Q(0)
    for value in comparison:
        running += value
        result.append(running)
    if running != caps[-1]:
        raise ArithmeticError('transport must preserve comparison total mass')
    return tuple(result)


def display_log(value):
    value = Q(value)
    return '-inf' if not value else f'{log2(value.numerator)-log2(value.denominator):.8f}'


def authenticated_bch_dimensions(last_lp=104, refined=True):
    """Freshly check retained BCH inputs and uniform shortened-dimension caps."""
    integer(last_lp, 0, 256, 'last_lp')
    if type(refined) is not bool:
        raise ValueError('refined must be bool')
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from bch_joint_support import authenticated_caps
    from shortened_bound import dimension_caps
    spectrum = authenticated_caps()
    if next(i for i in range(1, 257) if spectrum[i]) != 38:
        raise ArithmeticError('unexpected authenticated minimum support premise')
    dimensions = dimension_caps(last_lp=last_lp)
    if refined:
        from shortening_polynomial import improve_dimensions
        from dual_shortening import improve_dimensions as dual_dimensions
        dimensions = dual_dimensions(improve_dimensions(dimensions))
    return dimensions


def authenticated_bch_cdf(last_lp=104, refined=True, width=8):
    """Freshly check retained BCH inputs, then return expected CDF and premises.

    This authenticates only the local outer-count interface. A caller that
    composes it with the inner must independently verify the whole first
    moment and all setup-distribution assumptions. Default width8 preserves
    the original GL32 interface; other widths use full GL(4*width,2) maps.
    """
    integer(width, 1, 256, 'width')
    if 256 % width:
        raise ValueError('width must divide256')
    dimensions = authenticated_bch_dimensions(last_lp, refined)
    ranks = canonical_rank_caps(dimensions, width)
    caps = total_caps(ranks)
    transformed = transport_cdf(caps, full_block(width))
    return transformed, dict(schema=('packed-canonical-full32-expected-cdf-1' if width == 8
                                    else 'packed-canonical-fullblock-expected-cdf-1'),
        local_distribution=f'independent uniformGL{4*width} or nonzeroGF2^{4*width} scalar per group/block',
        block_width=width, block_rows=4, groups_per_outer_word=256//width,
        last_lp=last_lp, refined=refined, shortened_dimensions=list(dimensions),
        canonical_rank_cdfs=[list(row) for row in ranks], canonical_cdf=list(caps),
        note='Expected cumulative nonzero-message counts, not shell differences; no complete distance certificate.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--last-lp', type=int, default=0, help='Exact checked LP proposals through this support; zero uses elementary bounds only')
    parser.add_argument('--refined', action='store_true', help='Also run retained primal/dual dimension refinements')
    args = parser.parse_args()
    if not 0 <= args.last_lp <= 256:
        parser.error('last-lp must lie in0..256')
    transformed, premises = authenticated_bch_cdf(args.last_lp, args.refined)
    dimensions, ranks, caps = (premises[key] for key in
        ('shortened_dimensions', 'canonical_rank_cdfs', 'canonical_cdf'))
    print('Canonical fixed8-column blocks; fresh independent fullGL32/GF2^32 pergroup/block.')
    print('H,shortened_dimension,log2_total_CDF,log2_rank4_CDF')
    for h in (4, 5, 6, 8, 10, 12, 14, 15, 16, 18, 20, 24, 28, 32):
        print(h, dimensions[8*h], display_log(caps[h]), display_log(ranks[3][h]), sep=',')
    print('output_support,log2_expected_CDF_cap')
    for u in (32, 38, 64, 80, 96, 104, 112, 115, 120, 128, 160, 192, 224, 256):
        print(u, display_log(transformed[u]), sep=',')
    print('Only a one-group expected CDF: no inner composition or full-distance claim.')


if __name__ == '__main__':
    main()

"""Exact incidence refinements of canonical block counts; no whole-code claim.

The pure interface requires uniform shortened-dimension bounds as premises.
The authenticated interface freshly checks those bounds and the production
five-octet exclusion. Neither interface changes the binary code or its basis.
"""
from math import comb

from canonical_counts import (authenticated_bch_dimensions, canonical_rank_caps,
                              rank_total, total_caps, transport_cdf)
from local_models import full_block
from outer_hill_octets import check_production


def incidence_rank_caps(dimensions, width=8, rows=4, original_rank_cdfs=None):
    """Bound rank-r nonzero row tuples supported on at most h whole blocks.

    dimensions[u] must bound the shortened-code dimension on EVERY u-subset.
    A tuple with s<=h active blocks lies in C(g-s,t-s) sets of t blocks.
    For h<=t, this multiplicity is at least C(g-h,t-h). Counting incidences
    with all t-block sets therefore gives

      F_r(h) <= floor(C(g,t)*rank_total(dimensions[width*t],rows,r)
                      / C(g-h,t-h)).

    The floor is valid because F_r(h) counts deterministic tuples. We minimize
    over every t>=h and the previous canonical cap. The result is no larger
    than that cap. No premise is inferred from the resulting cumulative caps;
    in particular they are not reused as per-subset dimensions or counts.
    """
    dimensions = tuple(dimensions)
    previous = canonical_rank_caps(dimensions, width, rows, original_rank_cdfs)
    groups = (len(dimensions)-1)//width
    result = []
    for rank, old in enumerate(previous, 1):
        incidences = [comb(groups, t)*rank_total(dimensions[width*t], rows, rank)
                      for t in range(groups+1)]
        row = [0]
        for h in range(1, groups+1):
            row.append(min(old[h], *(incidences[t]//comb(groups-h, t-h)
                                     for t in range(h, groups+1))))
        for h in range(groups-1, -1, -1):
            row[h] = min(row[h], row[h+1])
        result.append(tuple(row))
    return tuple(result)


def production_rank_caps(dimensions, proof):
    """Apply incidence caps and a freshly produced exact five-octet result.

    This pure composition does not authenticate ``proof``. Callers must pass
    the result of check_production(), not an unchecked saved receipt.
    """
    if (len(dimensions) != 257 or dimensions[-1] != 128
            or proof.get('schema') != 'canonical-octet-h5-exact-1'
            or proof.get('length') != 256 or proof.get('dimension') != 128
            or proof.get('block_width') != 8 or proof.get('minimum_distance') != 38):
        raise ValueError('matching production dimensions and fresh exact octet proof required')
    limits = proof.get('rank_cumulative_caps_at5')
    if (not isinstance(limits, (list, tuple)) or len(limits) != 4
            or any(type(value) is not int or value < 0 for value in limits)):
        raise ValueError('four nonnegative exact rank caps required')
    ranks = [list(row) for row in incidence_rank_caps(dimensions)]
    for row, limit in zip(ranks, limits):
        for h in range(6):
            row[h] = min(row[h], limit)
    return tuple(map(tuple, ranks))


def authenticated_bch_cdf(last_lp=104, refined=True):
    """Freshly authenticate H5+incidence counts and return (CDF, premises).

    The returned CDF bounds expected cumulative nonzero-message counts after
    independent uniform GL32 maps on consecutive eight-column blocks. It is
    not a shell measure or a certificate for any inner construction.
    """
    dimensions = authenticated_bch_dimensions(last_lp, refined)
    proof = check_production()
    ranks = production_rank_caps(dimensions, proof)
    caps = total_caps(ranks)
    transformed = transport_cdf(caps, full_block(8))
    premises = dict(
        schema='packed-canonical-full32-h5-incidence-expected-cdf-1',
        local_distribution='independent uniformGL32 or nonzeroGF2^32 scalar per group/block',
        block_width=8, block_rows=4, groups_per_outer_word=32,
        last_lp=last_lp, refined=refined, shortened_dimensions=list(dimensions),
        canonical_rank_cdfs=[list(row) for row in ranks], canonical_cdf=list(caps),
        count_refinements=dict(
            incidence='all supersets t>=h; exact integer multiplicity division',
            exact_h5=proof),
        note='Expected cumulative nonzero-message counts; no complete distance certificate.')
    return transformed, premises

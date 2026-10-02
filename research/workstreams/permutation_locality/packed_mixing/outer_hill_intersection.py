"""Optional EKR refinement of canonical counts; no whole-code certificate.

Only the family of maximum-dimension shortenings is proved intersecting.
The next-lower-dimension family is not assumed to be intersecting. This
module preserves the H5 exclusion and the original uniform dimension table.
"""
from math import comb

from canonical_counts import authenticated_bch_dimensions, rank_total, total_caps, transport_cdf
from local_models import full_block, integer
from outer_hill_incidence import production_rank_caps as incidence_production_rank_caps
from outer_hill_octets import check_production


EKR_REFERENCE = 'https://www.renyi.hu/~p_erdos/1961-07.pdf'


def ekr_incidence_numerators(groups, dimension, subset_size, high_dimension,
                             intersection_dimension, rows=4):
    """Bound rank-r tuple incidences over all subset_size-block sets.

    Premises: the ambient binary code has dimension ``dimension``; every
    selected block set has shortened dimension at most ``high_dimension``;
    every intersection of two such sets with disjoint complements has
    shortened dimension at most ``intersection_dimension``. The caller must
    authenticate these uniform premises for its particular code.

    Two high-dimensional shortenings intersect in dimension at least
    2*high_dimension-dimension. If this exceeds intersection_dimension,
    their complements cannot be disjoint. EKR bounds this intersecting
    family of k=groups-subset_size subsets by C(groups-1,k-1), provided
    1<=k<=groups/2. All remaining shortenings have dimension at most high-1.

    The output for each rank r is
      C(groups,subset_size)*R_r(high-1)
      + C(groups-1,k-1)*(R_r(high)-R_r(high-1)).
    EKR is Theorem1 in the original paper linked by EKR_REFERENCE. Its
    proof is in Section5; see also OUTER_REFINEMENTS.md for this application.
    """
    integer(groups, 2, 256, 'groups')
    integer(dimension, 1, 256, 'dimension')
    integer(subset_size, 1, groups-1, 'subset_size')
    integer(high_dimension, 1, dimension, 'high_dimension')
    integer(intersection_dimension, 0, dimension, 'intersection_dimension')
    integer(rows, 1, 16, 'rows')
    complement_size = groups-subset_size
    if 2*complement_size > groups:
        raise ValueError('EKR requires complement size at most half the universe')
    if 2*high_dimension-dimension <= intersection_dimension:
        raise ValueError('strict dimension contradiction required for intersecting complements')
    subsets = comb(groups, subset_size)
    family_cap = comb(groups-1, complement_size-1)
    return tuple(subsets*rank_total(high_dimension-1, rows, rank)
                 + family_cap*(rank_total(high_dimension, rows, rank)
                               - rank_total(high_dimension-1, rows, rank))
                 for rank in range(1, rows+1))


def production_rank_caps(dimensions, proof):
    """Intersect fresh H5/incidence caps with the optional production EKR cap.

    As in the existing pure composition, dimensions and proof are premises,
    not authenticated saved receipts. Use authenticated_bch_cdf for fresh
    verification. The global uniform dimension table is never modified.
    """
    dimensions = tuple(dimensions)
    previous = incidence_production_rank_caps(dimensions, proof)
    # The preceding call checks length257, ambient dimension128, exact H5
    # metadata, and the integer/monotone dimension interface.
    if dimensions[64] > 2 or dimensions[160] > 66:
        raise ValueError('EKR production premises require d(64)<=2 and d(160)<=66')
    numerators = ekr_incidence_numerators(32, 128, 20, 66, dimensions[64])
    result = []
    for old, numerator in zip(previous, numerators):
        row = list(old)
        for h in range(1, 21):
            # Every tuple with support s<=h occurs in at least this many
            # 20-subsets. Counts are integral before random GL32 transport.
            row[h] = min(row[h], numerator//comb(32-h, 20-h))
        for h in range(31, -1, -1):
            row[h] = min(row[h], row[h+1])
        result.append(tuple(row))
    return tuple(result)


def authenticated_bch_cdf(last_lp=104, refined=True):
    """Return freshly authenticated expected CDF and H5/incidence/EKR premises.

    Maps are independent uniform GL32 maps on canonical four-row/eight-column
    blocks. This is an expected cumulative nonzero-message count, not a shell
    measure. A whole-code claim still requires fresh independent inner proof.
    """
    dimensions = authenticated_bch_dimensions(last_lp, refined)
    proof = check_production()
    ranks = production_rank_caps(dimensions, proof)
    caps = total_caps(ranks)
    transformed = transport_cdf(caps, full_block(8))
    numerators = ekr_incidence_numerators(32, 128, 20, 66, dimensions[64])
    premises = dict(
        schema='packed-canonical-full32-h5-incidence-ekr-expected-cdf-1',
        local_distribution='independent uniformGL32 or nonzeroGF2^32 scalar per group/block',
        block_width=8, block_rows=4, groups_per_outer_word=32,
        last_lp=last_lp, refined=refined, shortened_dimensions=list(dimensions),
        canonical_rank_cdfs=[list(row) for row in ranks], canonical_cdf=list(caps),
        count_refinements=dict(
            incidence='all supersets t>=h; exact integer multiplicity division',
            exact_h5=proof,
            ekr=dict(schema='canonical-t20-intersection-incidence-1',
                theorem='Erdos-Ko-Rado uniform intersecting families, Theorem1',
                reference=EKR_REFERENCE, ambient_dimension=128,
                groups=32, subset_size=20, complement_size=12,
                subset_dimension_cap=dimensions[160], high_dimension=66,
                other_dimension_cap=65, intersection_blocks=8,
                intersection_dimension_cap=dimensions[64],
                forced_intersection_dimension=4,
                maximum_high_dimension_subsets=comb(31, 11),
                rank_incidence_numerators=list(numerators),
                assumes_lower_dimension_family_intersecting=False)),
        note='Expected cumulative nonzero-message counts; no complete distance certificate.')
    return transformed, premises

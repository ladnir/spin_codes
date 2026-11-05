"""Exact deterministic shortening cap from all regional cardinalities.

There are B=2k regions and M outer rows. Region j selects a_j of its M
coordinates. For any resulting binary incidence matrix, let T_i be row i's
occupancy. A direct sum of MDS [B,k] rows has shortening

    D = sum_i max(T_i-k, 0).

The maximum over all binary incidence matrices with these column counts is
the sum of the k smallest a_j. Fix the k regions with largest counts. Each
row's excess is at most its selected coordinates outside those regions.
Summing proves the bound. Conversely, put each column's selected entries
in its first a_j rows. After sorting columns by decreasing a_j, each row's
excess equals its number of selected entries after the first k columns.
These nested supports attain the bound.

This is an exact maximum for the incidence-matrix relaxation. The bound is
valid for every packet route, but its attaining matrix need not be attainable
by a particular packet route. It is not a distance certificate or a route
probability estimate. No field coefficients or inner ranks enter the result.

The old full/touched-region secant bounds the same matrices, so it cannot
improve this cap. Completing the k largest regions to M and applying that
secant gives exactly this cap. Thus optimal virtual completion gives the
same result without searching over completion sets.
"""


def _counts(region_counts, row_count):
    counts = tuple(region_counts)
    if not isinstance(row_count, int) or row_count <= 0:
        raise ValueError("row_count must be a positive integer")
    if not counts or len(counts) % 2:
        raise ValueError("the number of regions must be positive and even")
    if any(not isinstance(a, int) or not 0 <= a <= row_count for a in counts):
        raise ValueError("regional counts must be integers in [0,row_count]")
    return counts


def regional_shortening_cap(region_counts, row_count):
    """Return the exact incidence-matrix maximum, in O(B log B) time."""
    counts = _counts(region_counts, row_count)
    return sum(sorted(counts)[:len(counts) // 2])


def secant_shortening_cap(region_counts, row_count):
    """Return the older cap retaining only total/full/touched counts."""
    counts = _counts(region_counts, row_count)
    half = len(counts) // 2
    full = sum(a == row_count for a in counts)
    touched = sum(a > 0 for a in counts)
    selected = sum(counts)
    if touched <= half:
        return 0
    if full >= half:
        return selected - row_count * half
    return (selected - row_count * full) * (touched - half) // (touched - full)


def virtual_completion_cap(region_counts, row_count, completed_regions):
    """Complete specified regions, then apply the old secant cap.

    Completion enlarges the selected coordinate set and hence its shortening.
    This is a valid upper bound on the original support, not an exact value.
    """
    counts = list(_counts(region_counts, row_count))
    for j in completed_regions:
        if not isinstance(j, int) or not 0 <= j < len(counts):
            raise ValueError("completed region index is out of range")
        counts[j] = row_count
    return secant_shortening_cap(counts, row_count)


def largest_half_regions(region_counts, row_count):
    """Return a deterministic maximizing information set of regions."""
    counts = _counts(region_counts, row_count)
    return tuple(sorted(range(len(counts)), key=lambda j: (-counts[j], j))
                 [:len(counts) // 2])

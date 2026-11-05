"""Exact source-group strata for alternating adjacent pair matchings.

For an indexed even cycle of size n and a support G of weight r, let a,b
count full pairs in the two matchings and F=max(a,b). This helper counts
the F strata, without quotienting rotations or sampling new geometry.
Its explicit reduced-weight budget prevents an accidental dense sweep.
"""
from functools import lru_cache
from math import comb

from hd_shifted_general_count import exact_joint_count


def _checked(size, weight, max_reduced_weight):
    if type(size) is not int or size < 4 or size % 2:
        raise ValueError('An even indexed cycle of at least four positions is required')
    if type(weight) is not int or not 0 <= weight <= size:
        raise ValueError('An integer support weight in the cycle is required')
    if type(max_reduced_weight) is not int or max_reduced_weight < 0:
        raise ValueError('An explicit nonnegative integer reduced-weight budget is required')
    reduced = min(weight, size-weight)
    if reduced > max_reduced_weight:
        raise ValueError('Reduced support weight exceeds the explicit work budget')
    return reduced


def _joint_count(size, weight, a, b):
    """Sum the exact oriented classes; weight is already budget-checked."""
    pairs = size//2
    low = max(0, weight-pairs, a, b)
    high = min(pairs, weight, weight-a, weight-b)
    return sum(exact_joint_count(pairs, even, weight-even, a, b)
               for even in range(low, high+1))


@lru_cache(maxsize=32)
def _small_joint(size, weight):
    rows = []
    for a in range(weight//2+1):
        for b in range(a+1):
            count = _joint_count(size, weight, a, b)
            if count:
                rows.append(((a, b), count))
                if a != b:
                    rows.append(((b, a), count))
    if sum(count for _, count in rows) != comb(size, weight):
        raise ArithmeticError('Joint matching strata do not partition indexed supports')
    return tuple(sorted(rows))


def matching_joint_counts(size, weight, *, max_reduced_weight=32):
    """Return {(a,b): exact indexed-support count} for both fixed matchings.

    The first matching is (0,1),(2,3),...; the second is
    (1,2),(3,4),...,(size-1,0). Count symmetry does not permit reversing
    a transfer-matrix schedule. Empty classes are omitted.
    """
    reduced = _checked(size, weight, max_reduced_weight)
    shift = weight-size//2 if weight > size//2 else 0
    return {(a+shift, b+shift): count for (a, b), count in _small_joint(size, reduced)}


@lru_cache(maxsize=32)
def _small_cumulative(size, weight):
    # Rotation by one position swaps a and b for SUPPORT COUNTS. No
    # symmetry of initial-state transfer operators is assumed here.
    cumulative, total = [], 0
    for full in range(weight//2+1):
        border = _joint_count(size, weight, full, full)
        border += 2*sum(_joint_count(size, weight, full, other)
                        for other in range(full))
        total += border
        cumulative.append((full, total))
    if total != comb(size, weight):
        raise ArithmeticError('Matching-max strata do not partition indexed supports')
    return tuple(cumulative)


def matching_max_cumulative_counts(size, weight, *, max_reduced_weight=32):
    """Return {u: #{G: |G|=weight and max(a(G),b(G))<=u}} exactly.

    Keys run from the smallest feasible threshold through floor(weight/2).
    Counts are zero below the first key and binomial(size,weight) above the
    last. Complementing a support shifts both a and b by the same amount.
    """
    reduced = _checked(size, weight, max_reduced_weight)
    shift = weight-size//2 if weight > size//2 else 0
    return {full+shift: count for full, count in _small_cumulative(size, reduced)}


def matching_max_counts(size, weight, *, max_reduced_weight=32):
    """Return {F: exact number of indexed supports with max(a,b)=F}.

    Empty strata are omitted. This is a support count, not a probability
    law for resampling the active groups between regions or RS columns.
    """
    cumulative = matching_max_cumulative_counts(size, weight,
                                                max_reduced_weight=max_reduced_weight)
    result, previous = {}, 0
    for full, count in cumulative.items():
        difference = count-previous
        if difference < 0:
            raise ArithmeticError('Cumulative support counts decreased')
        if difference:
            result[full] = difference
        previous = count
    return result


def subset_full_pair_range(source_weight, source_full, subset_weight):
    """All possible full-pair counts of U subset G in one fixed matching.

    A source G has F double pairs and r-2F singleton pairs. A subset of
    size m has f full pairs only if m-2f <= (F-f)+(r-2F). Together with
    0<=f<=F and 2f<=m this gives the returned interval. Every integer in
    that interval can be attained by selecting suitable positions of G.
    """
    if (any(type(value) is not int for value in (source_weight, source_full, subset_weight))
            or source_weight < 0 or not 0 <= source_full <= source_weight//2
            or not 0 <= subset_weight <= source_weight):
        raise ValueError('Valid integer source and subset matching counts are required')
    return max(0, subset_weight-source_weight+source_full), min(source_full, subset_weight//2)


if __name__ == '__main__':
    import argparse
    import json
    from math import log2
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=8192)
    parser.add_argument('--weight', type=int, default=17)
    parser.add_argument('--max-reduced-weight', type=int, default=32)
    args = parser.parse_args()
    counts = matching_max_counts(args.size, args.weight,
                                 max_reduced_weight=args.max_reduced_weight)
    total = comb(args.size, args.weight)
    print(json.dumps(dict(size=args.size, weight=args.weight, indexed_supports=True,
        total_count=total, strata=[dict(full_max=full, count=count,
            log2_count=log2(count), fraction=count/total) for full, count in counts.items()],
        exact_integer_counts=True, directed_distance_certificate=False)), flush=True)

"""Exact setup-failure bound for the two-moment map of randomized RS rows.

This is a rank lemma, NOT a distance certificate or a numerical gate.
Fix a set A of a distinct row labels u_i, and let every row be RS[n,k].
Before any native gauge normalization, d_ij are independent uniform NZ
elements of GF(q).  The map from the a*k message coordinates to 2*n
syndrome coordinates is

    c -> ((sum_i d_ij*c_ij, sum_i u_i*d_ij*c_ij))_j.

A nonzero left-null candidate lambda determines row vectors
v_i(j)=lambda_0j+u_i*lambda_1j.  At most one v_i is identically zero.
For every other row, diag(d_i)*v_i must belong to RS[n,k]^perp.  Its
probability is at most (q-1)^(-k): for exact support w, the shortened
dual code has dimension w-k and at most (q-1)^(w-k) full-support words.
The event is impossible when 1<=w<=k.  Different rows use independent
diagonals.  Count projective lambda once to obtain the bound below.

On full rank, all q^(2*n) syndromes have exactly q^(a*k-2*n) preimages
in the space allowing zero constituent rows.  This is a FIXED-SETUP
linear-algebra statement.  It is not a statement of independent fiber
supports, or of independence between columns conditioned on a fiber.

Multiplying all d_ij in a column by the same nonzero scalar rescales
that column's two syndrome coordinates.  The rank event is invariant,
so the ideal-law bound transfers through that native gauge.  It does
not assert independence of the normalized diagonal coefficients.
"""
from fractions import Fraction
from math import comb


def _parameters(q, length, dimension):
    if (any(type(v) is not int for v in (q, length, dimension))
            or q < 2 or q & (q-1) or not 1 <= dimension < length <= q):
        raise ValueError('GF(2^b) and proper RS[n,k] parameters required')


def projective_classes(arity, *, q=256, length=16):
    """Return counts (no zero induced row, one zero induced row).

    Row labels are distinct.  Nonzero lambda cannot induce zero on two
    rows, so the arity exceptional sets are disjoint.
    """
    if (any(type(v) is not int for v in (arity, q, length))
            or q < 2 or q & (q-1) or not 1 <= arity <= q
            or not 1 <= length <= q):
        raise ValueError('Distinct row labels and a positive length required')
    exceptional = arity * (q**length-1)
    generic = q**(2*length)-1-exceptional
    generic_count, generic_remainder = divmod(generic, q-1)
    exceptional_count, exceptional_remainder = divmod(exceptional, q-1)
    if min(generic_count, exceptional_count) < 0 or generic_remainder or exceptional_remainder:
        raise ArithmeticError('Projective class count failed')
    return generic_count, exceptional_count


def fixed_arity_rank_failure(arity, *, q=256, length=16, dimension=8):
    """An exact rational upper bound for rank strictly below 2*n.

    Return the raw projective union bound; it can exceed one for small
    arities.  No independence between null candidates is assumed.
    """
    _parameters(q, length, dimension)
    generic, exceptional = projective_classes(arity, q=q, length=length)
    return (Fraction(generic, (q-1)**(dimension*arity))
            + Fraction(exceptional, (q-1)**(dimension*(arity-1))))


def all_subset_rank_failure(*, groups=8192, row_labels=16, minimum_arity=5,
                            q=256, length=16, dimension=8):
    """Union over all groups and all row subsets of the requested arities.

    This uses a union bound, not independence between overlapping row
    sets.  The only claimed event is simultaneous full syndrome rank.
    """
    _parameters(q, length, dimension)
    if (any(type(v) is not int for v in (groups, row_labels, minimum_arity))
            or groups < 1 or not 1 <= minimum_arity <= row_labels <= q):
        raise ValueError('A positive family of distinct row subsets required')
    rows = tuple((a, comb(row_labels, a), fixed_arity_rank_failure(
        a, q=q, length=length, dimension=dimension))
        for a in range(minimum_arity, row_labels+1))
    union = groups * sum((count * min(Fraction(1), bound)
                         for _, count, bound in rows), Fraction(0))
    return dict(q=q, length=length, dimension=dimension, groups=groups,
                row_labels=row_labels, minimum_arity=minimum_arity,
                per_arity=rows, union_upper=min(Fraction(1), union),
                raw_union=union, distance_certificate=False)


def full_rank_fiber_dimension(arity, *, length=16, dimension=8):
    """Dimension of EACH fiber, conditional on the stated full-rank event."""
    if (any(type(v) is not int for v in (arity, length, dimension))
            or not 1 <= dimension < length or arity*dimension < 2*length):
        raise ValueError('Full syndrome rank must be dimensionally possible')
    return arity*dimension-2*length


if __name__ == '__main__':
    import json
    record = all_subset_rank_failure()
    record['per_arity'] = [dict(arity=a, subsets=count, rank_failure_upper=str(bound))
                           for a, count, bound in record['per_arity']]
    record['union_upper'] = str(record['union_upper'])
    record['raw_union'] = str(record['raw_union'])
    print(json.dumps(record, separators=(',', ':')))

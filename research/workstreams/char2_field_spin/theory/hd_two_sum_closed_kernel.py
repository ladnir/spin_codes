"""Exact uniform-input kernel for inverse-disjoint affine two-sum labels.

Applies to the checked g8 labels8..15 and g16 labels06+GF16. These are
affine GF2 spaces disjoint from their inverses. It is a counting law for
uniform input on a FIXED coordinate space, including zero coordinates.
It is not the exact nonzero-input law or a whole-code certificate.
"""
from functools import lru_cache
from math import comb


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(None)
def zero_one_moment(size, q=256):
    """NZ assignments with one full-support linear equation equal to zero."""
    if type(size) is not int or size < 0:
        raise ValueError('Nonnegative support size required')
    return ((q-1)**size+(q-1)*(-1)**size)//q


@lru_cache(None)
def zero_two_moments(size, q=256):
    """NZ assignments to two Vandermonde checks with distinct labels."""
    if type(size) is not int or not 0 <= size <= 16:
        raise ValueError('Support size0..16 required')
    return sum((-1)**(size-j)*comb(size, j)*q**max(j-2, 0) for j in range(size+1))


def branch_counts(arity, inside_weight, group=8):
    if group not in (8, 16) or not 0 <= inside_weight <= arity <= group:
        raise ValueError('Valid g8/g16 arity and inside weight required')
    quiet = zero_two_moments(inside_weight)
    hole = zero_one_moment(inside_weight)-quiet
    full = 255**inside_weight-quiet-(group-arity)*hole
    if min(quiet, hole, full) < 0:
        raise ArithmeticError('A purported support branch has negative count')
    return quiet, hole, full


def uniform_mask_count(active, target, group=8):
    """One exact N_A(T); overlapping branches are ADDED for complement0/1."""
    if (group not in (8, 16) or any(type(x) is not int or not 0 <= x < 1 << group for x in (active, target))):
        raise ValueError('Valid g8/g16 coordinate masks required')
    outside = ((1 << group)-1) ^ active
    selected = target & outside
    quiet, hole, full = branch_counts(active.bit_count(), (active & target).bit_count(), group)
    return ((quiet if selected == 0 else 0)
            +(hole if (outside ^ selected).bit_count() == 1 else 0)
            +(full if selected == outside else 0))


def uniform_mask_counts(active, group=8):
    result = tuple(uniform_mask_count(active, target, group) for target in range(1 << group))
    if result[0] != 1 or sum(result) != 256**active.bit_count():
        raise ArithmeticError('Uniform kernel normalization failed')
    return result


def arity_power_sum(arity, weight, group=8, columns=32):
    """Exact sum_{|A|=arity} N_A(T)^columns for ANY fixed |T|=weight.

    The result depends only on weight, not the positions in T. Branches
    have disjoint outside support sizes when group-arity>=2. The two edge
    arities combine coincident branches BEFORE taking the power.
    """
    if (group not in (8, 16) or type(columns) is not int or columns < 1
            or not 0 <= arity <= group or not 0 <= weight <= group):
        raise ValueError('Valid arity, target weight and positive power required')
    if arity == group:
        return 255**(weight*columns)
    if arity == group-1:
        missing = (group-weight)*zero_one_moment(weight)**columns
        present = 0 if weight == 0 else weight*(255**(weight-1)-zero_one_moment(weight-1))**columns
        return missing+present
    result = choose(group-weight, arity-weight)*zero_two_moments(weight)**columns
    b = arity-group+weight
    if 0 <= b <= arity:
        result += choose(weight, b)*branch_counts(arity, b, group)[2]**columns
    b += 1
    if 0 <= b <= arity:
        result += (group-weight)*choose(weight, b)*branch_counts(arity, b, group)[1]**columns
    return result


def arity_holder_moments(arity, group=8, columns=32, dimension=16):
    """W_a(T)^columns, represented by(group+1) numerators and one denominator.

    Row sets remain fixed across columns under Holder. Their once-only
    multiplicity is inside the power sum; do not multiply by C(group,a).
    These weights are not a normalized probability distribution.
    """
    if not 0 < dimension < columns:
        raise ValueError('A positive proper constituent dimension is required')
    return (tuple(arity_power_sum(arity, w, group, columns) for w in range(group+1)),
            255**((columns-dimension)*arity))

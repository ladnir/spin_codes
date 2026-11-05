"""Exact labelled support counts for H0=I+11^T on an even-size block.

The default size is16 and field order256.  This is a new outer mixer,
not the previous H* ensemble.  Counts refer either to all vectors whose
support is CONTAINED in A, or to vectors whose support is EXACTLY A.
No field-amplitude enumeration, regional provider, or distance gate is
performed.  The output dictionaries retain coordinate labels.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb

CONSTRUCTION = 'G16-RS16x8-D2-preferred-random-ratio-pair-H0-D1-research'


def _parameters(mask, q, size):
    if (type(q) is not int or q < 2 or q & (q-1) or type(size) is not int
            or size < 2 or size % 2 or type(mask) is not int or not 0 <= mask < 1 << size):
        raise ValueError('A characteristic-two field, even size, and labelled support mask required')


def nonzero_sum_count(length, *, target_zero=True, q=256):
    """Count NZ field tuples with sum0, or sum equal to any fixed NZ value."""
    if type(length) is not int or length < 0 or type(q) is not int or q < 2 or q & (q-1):
        raise ValueError('A nonnegative length and power-of-two field order required')
    numerator = (q-1)**length+((q-1) if target_zero else -1)*(-1)**length
    answer, remainder = divmod(numerator, q)
    if remainder or answer < 0:
        raise AssertionError('NZ additive-character count is not a nonnegative integer')
    return answer


def nonzero_nonone_sum_count(length, *, target_zero=True, q=256):
    """Count tuples outside{0,1} summing to0 or1; these agree for length>0."""
    if type(length) is not int or length < 0 or type(q) is not int or q < 2 or q & (q-1):
        raise ValueError('A nonnegative length and power-of-two field order required')
    if length == 0:
        return int(target_zero)
    numerator = (q-2)**length+(q//2-1)*(-2)**length
    answer, remainder = divmod(numerator, q)
    if remainder or answer < 0:
        raise AssertionError('Punctured additive-character count is invalid')
    return answer


@lru_cache(maxsize=128)
def weight_coefficients(arity, *, mode='allowed', q=256):
    """Compact per-labelled-subset counts indexed by INSIDE output weight.

    These are not weight-enumerator coefficients: multiply by choose(a,b)
    only when explicitly summing all labelled inside subsets of size b.
    """
    if type(arity) is not int or arity < 0 or mode not in ('allowed', 'forced'):
        raise ValueError('A nonnegative arity and allowed/forced mode required')
    quiet = tuple(nonzero_sum_count(b, q=q) if mode == 'allowed' or b == arity else 0
                  for b in range(arity+1))
    count = nonzero_sum_count if mode == 'allowed' else nonzero_nonone_sum_count
    nonquiet = tuple((q-1)*count(b, target_zero=bool(arity % 2), q=q) for b in range(arity+1))
    return quiet, nonquiet


def _subsets(mask):
    subset = mask
    while True:
        yield subset
        if not subset:
            return
        subset = (subset-1) & mask


def support_branches(allowed_mask, *, mode='allowed', q=256, size=16):
    """Sparse exact quiet/nonquiet output-mask dictionaries for fixed A.

    allowed: x belongs to F^A; zero input is counted once.
    forced: each x_i for i in A is NZ; coordinates outside A are zero.
    For A=empty both modes contain the one zero input.

    The parity-sum branches must be ADDED when they have the same output
    support, notably when A is the whole coordinate set.
    """
    _parameters(allowed_mask, q, size)
    if mode not in ('allowed', 'forced'):
        raise ValueError("mode must be 'allowed' or 'forced'")
    a = allowed_mask.bit_count()
    outside = ((1 << size)-1) ^ allowed_mask
    quiet_coefficients, nonquiet_coefficients = weight_coefficients(a, mode=mode, q=q)
    quiet, nonquiet = {}, {}
    for inside in _subsets(allowed_mask):
        b = inside.bit_count()
        count = quiet_coefficients[b]
        if count:
            quiet[inside] = count
        count = nonquiet_coefficients[b]
        if count:
            nonquiet[outside | inside] = count
    mass = q**a if mode == 'allowed' else (q-1)**a
    if sum(quiet.values())+sum(nonquiet.values()) != mass:
        raise AssertionError('The exact input mass changed')
    return {'quiet': dict(sorted(quiet.items())), 'nonquiet': dict(sorted(nonquiet.items()))}


def support_counts(allowed_mask, *, mode='allowed', q=256, size=16):
    """Sparse outputmask->integercount; at most2^(|A|+1) closed-form terms."""
    branches = support_branches(allowed_mask, mode=mode, q=q, size=size)
    result = Counter(branches['quiet'])
    result.update(branches['nonquiet'])
    return dict(sorted(result.items()))


def forced_support_law(mask, *, q=256, size=16):
    """Law after independent uniform-NZ D1 on a fixed EXACT input support."""
    counts = support_counts(mask, mode='forced', q=q, size=size)
    denominator = (q-1)**mask.bit_count()
    return {output: Fraction(count, denominator) for output, count in counts.items()}


def direct(values):
    """Apply the parity mixer; XOR is field addition in every binary basis."""
    values = tuple(values)
    if len(values) < 2 or len(values) % 2 or any(type(v) is not int or v < 0 for v in values):
        raise ValueError('An even-length tuple of nonnegative field representations required')
    parity = 0
    for value in values:
        parity ^= value
    return tuple(value ^ parity for value in values)


def delta_full_pairs(mask, delta, *, size=16):
    """Number of xor-delta row pairs fully contained in fixed A."""
    if (type(size) is not int or size < 2 or size & (size-1) or type(mask) is not int
            or not 0 <= mask < 1 << size or type(delta) is not int or not 0 < delta < size):
        raise ValueError('A binary-labelled row set and nonzero xor difference required')
    return sum(bool(mask >> i & 1) and bool(mask >> (i ^ delta) & 1)
               for i in range(size) if i < (i ^ delta))


def pair_profile(mask, *, ordered=True):
    """Fixed-A overlap invariants for preferred one-layer differences8..15.

    For one difference, scalar equal-marker kernels depend only on
    (arity, full_pairs).  ordered=False is legitimate only when the same
    scalar pair functional/markers apply to every column pair.  This
    function does not authorize collapsing chronological matrix operators.
    """
    _parameters(mask, 256, 16)
    profile = tuple(delta_full_pairs(mask, delta) for delta in range(8, 16))
    return profile if ordered else tuple(sorted(Counter(profile).items()))


def profile_classes(arities, *, ordered=True, max_masks=65535):
    """Exact labelled multiplicities; no numerical kernel is evaluated."""
    arities = tuple(sorted(set(arities)))
    if (not arities or any(type(a) is not int or not 1 <= a <= 16 for a in arities)
            or type(max_masks) is not int or max_masks < 1
            or sum(comb(16, a) for a in arities) > max_masks):
        raise ValueError('Explicit arities and a sufficient labelled-mask budget required')
    classes = Counter()
    for a in arities:
        for selected in combinations(range(16), a):
            mask = sum(1 << i for i in selected)
            classes[a, pair_profile(mask, ordered=ordered)] += 1
    return classes


def scalar_column_kernel(mask, marker, *, mode='allowed', q=256, size=16):
    """Exact one-column scalar WE, supplied only as a tiny integration aid."""
    marker = Fraction(marker)
    if marker < 0:
        raise ValueError('A nonnegative scalar marker required')
    return sum(count*marker**output.bit_count()
               for output, count in support_counts(mask, mode=mode, q=q, size=size).items())


def primitive_counts(size=16):
    if type(size) is not int or size < 2 or size % 2:
        raise ValueError('An even coordinate count required')
    return dict(binary_xors=2*size-1, variable_products=0, fixed_products=0,
                excludes_diagonals=True, excludes_pair_layer=True, excludes_RS=True)

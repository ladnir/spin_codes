"""Exact affine-flat geometry for the fixed bit-syndrome outer.

This is local kernel preparation only. It evaluates no regional provider,
distance contraction, or native encoder. Coordinates are the16 binary
labels; field arithmetic is AES GF256. All returned support masks retain
their actual labels. No affine symmetry of H* is assumed.

rank_cross_masks(O,S) means rank(V_O B V_S^T), INCLUDING overlapping O,S.
It equals rank(H*[O,S]) only when those masks are disjoint. Each argument
can be replaced by its binary affine hull. The resulting pair cache is
lazy and symmetric; no complete rank table is prepared at import.

For |A|<=7, A^c spans all five affine rows. Thus outside evaluation is
injective on W_A=B im(V_A^T). A coefficient vector in W_A has a global
zero set that is empty or an affine flat in F2^4. Exact containment
inversion over the308 possibilities counts the actual outside image,
without enumerating2^|A^c| supports or q^rank coefficient vectors.
The input-to-outside fiber factor q^(|A|-rank(V_A)) is NOT included.
"""
from functools import lru_cache

from gf8_packet_maps import rank
from hd_bit_syndrome_outer import affine_rank, matrix

BITS, SIZE, Q = 4, 16, 256
FULL = (1 << SIZE)-1


def _checked_mask(mask, bits=4):
    if (type(bits) is not int or not 1 <= bits <= 4 or type(mask) is not int
            or not 0 <= mask < (1 << (1 << bits))):
        raise ValueError('An integer mask of the declared binary label set is required')
    return mask


def _labels(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask ^= bit


def flat_hull_mask(mask, bits=4):
    """Smallest binary affine flat containing mask; hull(empty)=empty."""
    _checked_mask(mask, bits)
    return _flat_hull_cached(mask, bits)


@lru_cache(maxsize=65536)
def _flat_hull_cached(mask, bits):
    if not mask:
        return 0
    anchor = (mask & -mask).bit_length()-1
    pivots = [0] * bits
    basis = []
    for label in _labels(mask):
        value = label ^ anchor
        while value:
            lead = value.bit_length()-1
            if pivots[lead]:
                value ^= pivots[lead]
            else:
                pivots[lead] = value
                basis.append(value)
                break
    span = [0]
    for direction in basis:
        span += [value ^ direction for value in span]
    return sum(1 << (anchor ^ value) for value in span)


def affine_flats(bits=4):
    """All affine flats plus empty, sorted by cardinality then mask.

    At bits4 there are16 points,120 lines,140 planes,30 hyperplanes,
    and the full16-point space:307 nonempty flats,308 including empty.
    Binary subspaces are generated once and then translated into cosets.
    """
    _checked_mask(0, bits)
    return _affine_flats_cached(bits)


@lru_cache(maxsize=4)
def _affine_flats_cached(bits):
    size = 1 << bits
    subspaces, queue = {1}, [1]  # The zero subspace contains label0.
    for subspace in queue:
        elements = tuple(_labels(subspace))
        for direction in range(size):
            if subspace >> direction & 1:
                continue
            enlarged = subspace | sum(1 << (value ^ direction) for value in elements)
            if enlarged not in subspaces:
                subspaces.add(enlarged)
                queue.append(enlarged)
    flats = {0}
    for subspace in subspaces:
        remaining = (1 << size)-1
        elements = tuple(_labels(subspace))
        while remaining:
            anchor = (remaining & -remaining).bit_length()-1
            coset = sum(1 << (anchor ^ value) for value in elements)
            flats.add(coset)
            remaining &= ~coset
    return tuple(sorted(flats, key=lambda value: (value.bit_count(), value)))


@lru_cache(maxsize=308)
def _flat_basis(flat):
    """Labels giving a binary/GF256 basis of the augmented affine rows."""
    pivots = [0] * 5
    basis = []
    for label in _labels(flat):
        value = SIZE | label
        while value:
            lead = value.bit_length()-1
            if pivots[lead]:
                value ^= pivots[lead]
            else:
                pivots[lead] = value
                basis.append(label)
                break
    return tuple(basis)


@lru_cache(maxsize=1)
def _bilinear_matrix():
    # H*-I=V B V^T, including its diagonal. Characteristic two: -=+=xor.
    h = matrix()
    return tuple(tuple(value ^ int(i == j) for j, value in enumerate(row))
                 for i, row in enumerate(h))


@lru_cache(maxsize=47586)
def _rank_flat_pair(left, right):
    # 308*309/2 entries suffice for the whole symmetric cache, if requested.
    bilinear = _bilinear_matrix()
    return rank(tuple(tuple(bilinear[i][j] for j in _flat_basis(right))
                      for i in _flat_basis(left)))


def rank_cross_masks(outside_mask, source_mask):
    """rank(V_O B V_S^T); H*[O,S] only for disjoint masks.

    Overlap is intentionally permitted for global cross-column basis tests.
    Empty arguments return zero. The cache key uses both affine hulls, in
    sorted order, because the fixed B is symmetric. Values are exact.
    """
    _checked_mask(outside_mask)
    _checked_mask(source_mask)
    if not outside_mask or not source_mask:
        return 0
    left, right = flat_hull_mask(outside_mask), flat_hull_mask(source_mask)
    return _rank_flat_pair(min(left, right), max(left, right))


@lru_cache(maxsize=308)
def _zero_counts_for_hull(hull):
    """Frozen (global zero-set,count) pairs for coefficient space W_hull."""
    dimension = affine_rank(hull)
    exact = []
    for flat in reversed(affine_flats()):
        vanish = Q ** (dimension-rank_cross_masks(flat, hull))
        count = vanish-sum(value for larger, value in exact if not (flat & ~larger))
        if count < 0:
            raise AssertionError('Exact zero-flat containment inversion became negative')
        if count:
            exact.append((flat, count))
    if sum(value for _, value in exact) != Q ** dimension or dict(exact).get(FULL) != 1:
        raise AssertionError('Global affine-evaluation counts have the wrong mass')
    return tuple(sorted(exact))


def correction_zero_set_counts(source_mask):
    """Exact global zero sets of Vw, w in B im(V_A^T), each w counted once.

    Counts depend on the affine hull of A, but on the ACTUAL fixed B.
    They are not transported under an assumed affine symmetry of H*.
    Source empty is permitted and gives the identically zero correction.
    """
    _checked_mask(source_mask)
    return dict(_zero_counts_for_hull(flat_hull_mask(source_mask)))


def outside_image_support_counts(source_mask):
    """Exact {global outside support:count} for nonempty |A|<=7.

    This counts im(H*[A^c,A]), not the full input counting measure. For
    |A|<=7, the outside evaluation is injective on W_A; distinct correction
    coefficients are therefore distinct outside words. Aggregate global
    zero flats only by their intersection with A^c. Outside zero is kept.
    """
    _checked_mask(source_mask)
    if not 1 <= source_mask.bit_count() <= 7:
        raise ValueError('This injective outside-image interface requires arity1..7')
    outside = FULL ^ source_mask
    dimension = affine_rank(source_mask)
    if affine_rank(outside) != 5 or rank_cross_masks(outside, source_mask) != dimension:
        raise AssertionError('Low-arity outside evaluation is not injective')
    answer = {}
    for zeros, count in correction_zero_set_counts(source_mask).items():
        support = outside & ~zeros
        answer[support] = answer.get(support, 0)+count
    if sum(answer.values()) != Q ** dimension or answer.get(0) != 1:
        raise AssertionError('Actual outside-image count or zero word changed')
    return dict(sorted(answer.items()))

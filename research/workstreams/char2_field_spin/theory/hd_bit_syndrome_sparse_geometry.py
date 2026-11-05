"""Exact labelled geometry for bit-syndrome outer arities2..7.

For fixed A and O=A^c, the nine or more outside labels span all five
affine syndrome coordinates.  Consequently the actual cross map
H*[O,A] has rank rank(V_A), and its zero-outside fiber is precisely the
quiet kernel ker(V_A^T).  No ambient outside-image substitution is made.

The existing flat helper counts each correction coefficient once by its
GLOBAL zero set.  Restrict those zero sets to O and aggregate identical
outside supports.  Per-weight maxima give positive coefficient caps C_A;
their constant is exactly1.  Deleting that constant gives the nonquiet
outside branch.  Information covers bound every translated inside fiber.

Only quiet ranks and covers are transported along binary affine orbits.
Every actual outside image is counted at its labelled A.  All row sets,
cover positions, and ordered gaps survive preparation.  There is no
regional provider, contraction, RS normalization, source count, or claim
of a distance certificate here.  Integer counts and rational covers are
exact; NumPy arrays contain only integer class and gap indices.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb

import numpy as np

from hd_bit_syndrome_arity9_gate import gap_geometry, quiet_cover_orbits
from hd_bit_syndrome_flat_geometry import (
    FULL, correction_zero_set_counts, rank_cross_masks)
from hd_bit_syndrome_kernel import information_cover
from hd_bit_syndrome_outer import CONSTRUCTION, affine_rank

SIZE, Q = 16, 256


def _source(source):
    source = tuple(source)
    if (not 2 <= len(source) <= 7 or source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < SIZE for i in source)):
        raise ValueError('A sorted labelled source set of arity2..7 is required')
    return source


def aggregate_outside_supports(outside, zero_sets, coefficient_dimension, cross_rank,
                               *, size=SIZE, q=Q):
    """Restrict exact global zero-set counts, dividing the constant fiber.

    The coefficient-to-global-word map is injective.  Its restriction
    to the outside coordinates has the supplied cross rank.  These are
    mathematical premises; production obtains them from the flat helper.
    The generic size/q parameters serve tiny exact tests only.
    """
    if (type(size) is not int or not 1 <= size <= SIZE or type(q) is not int
            or q < 2 or q & (q-1) or type(outside) is not int
            or not 0 <= outside < 1 << size
            or type(coefficient_dimension) is not int or type(cross_rank) is not int
            or not 0 <= cross_rank <= coefficient_dimension <= size):
        raise ValueError('Valid finite-field dimensions and labelled outside mask required')
    zero_sets = dict(zero_sets)
    full = (1 << size)-1
    if (not zero_sets or any(type(mask) is not int or not 0 <= mask <= full
                            or type(count) is not int or count <= 0
                            for mask, count in zero_sets.items())
            or sum(zero_sets.values()) != q**coefficient_dimension
            or zero_sets.get(full) != 1):
        raise ValueError('Exact injective global zero-set counts are required')
    answer = {}
    for zeros, count in zero_sets.items():
        support = outside & ~zeros
        answer[support] = answer.get(support, 0)+count
    divisor = q**(coefficient_dimension-cross_rank)
    if any(count % divisor for count in answer.values()):
        raise AssertionError('Outside restriction does not have constant-size fibers')
    answer = {support: count // divisor for support, count in answer.items()}
    if answer.get(0) != 1 or sum(answer.values()) != q**cross_rank:
        raise AssertionError('Actual outside image mass or zero word changed')
    return dict(sorted(answer.items()))


def caps_by_weight(support_counts, length):
    """Maximum count per FIXED support of each weight, not total weight mass."""
    if type(length) is not int or not 0 <= length <= SIZE:
        raise ValueError('A support length0..16 is required')
    caps = [0]*(length+1)
    for mask, count in support_counts.items():
        if (type(mask) is not int or not 0 <= mask <= FULL or mask.bit_count() > length
                or type(count) is not int or count < 0):
            raise ValueError('Nonnegative exact global support counts required')
        weight = mask.bit_count()
        caps[weight] = max(caps[weight], count)
    return tuple(caps)


def actual_image_data(source):
    """One actual labelled outside image, its caps, and exact rank metadata."""
    source = _source(source)
    mask = sum(1 << i for i in source)
    outside = FULL ^ mask
    rank = affine_rank(mask)
    cross_rank = rank_cross_masks(outside, mask)
    if affine_rank(outside) != 5 or cross_rank != rank:
        raise AssertionError('Sparse outside evaluation must be injective')
    counts = aggregate_outside_supports(outside, correction_zero_set_counts(mask), rank, cross_rank)
    caps = caps_by_weight(counts, SIZE-len(source))
    if caps[0] != 1:
        raise AssertionError('The full outside cap must contain one zero word')
    return dict(support_counts=counts, full_caps=caps, nonquiet_caps=(0,)+caps[1:],
                cross_rank=cross_rank, fiber_dimension=len(source)-rank,
                zero_outside_equals_quiet=True, restriction_fiber_size=1)


def quiet_data(source, *, q=Q):
    """Exact quiet support counts on one A, with at most128 binary subsets.

    This small binary-matroid inversion is used only on representatives
    during complete preparation.  It is not a new field-rank DP.
    """
    source = _source(source)
    if type(q) is not int or q < 2 or q & (q-1):
        raise ValueError('A characteristic-two field order is required')
    masks = [0]*(1 << len(source))
    for local in range(1, len(masks)):
        bit = local & -local
        masks[local] = masks[local ^ bit] | (1 << source[bit.bit_length()-1])
    values = [q**(mask.bit_count()-affine_rank(mask)) for mask in masks]
    for i in range(len(source)):
        step = 1 << i
        for base in range(0, len(values), 2*step):
            for offset in range(step):
                values[base+step+offset] -= values[base+offset]
    dimension = len(source)-affine_rank(masks[-1])
    if min(values) < 0 or values[0] != 1 or sum(values) != q**dimension:
        raise AssertionError('Quiet support inversion failed')
    counts = {mask: count for mask, count in zip(masks, values) if count}
    return dict(support_counts=counts, caps=caps_by_weight(counts, len(source)), dimension=dimension)


def prepare_geometry(arity, *, rowsets=None, check_time=lambda: None):
    """Per-A covers/classes and pooled ordered gaps, without numeric kernels.

    The complete path reuses quiet_cover_orbits for exact cover transport.
    An explicit rowsets tuple is a small-test interface and is marked
    partial.  Outside-image caps always use the actual labelled source.
    """
    if type(arity) is not int or not 2 <= arity <= 7:
        raise ValueError('Arity2..7 is required')
    check_time()
    expected = tuple(combinations(range(SIZE), arity))
    if rowsets is None:
        orbit = quiet_cover_orbits(SIZE, arity, check_time=check_time)
        rows, covers = orbit['rowsets'], orbit['covers']
        quiet_classes = orbit['classes'].copy()
        quiet_caps = tuple(quiet_data(row)['caps'] for row in orbit['representatives'])
        representatives, orbit_sizes = orbit['representatives'], orbit['orbit_sizes']
    else:
        rows = tuple(_source(row) for row in rowsets)
        if (not rows or len(set(rows)) != len(rows) or any(len(row) != arity for row in rows)):
            raise ValueError('Distinct nonempty rowsets matching the declared arity required')
        exact_covers = tuple(information_cover(row) for row in rows)
        covers = tuple(tuple(cover[i] for i in row) for row, cover in zip(rows, exact_covers))
        representatives, orbit_sizes = rows, (1,)*len(rows)
        quiet_caps = tuple(quiet_data(row)['caps'] for row in rows)
        quiet_classes = np.arange(len(rows), dtype=np.int16)
    cap_classes, unique_caps, classes, ranks, dimensions = {}, [], [], [], []
    for index, (source, cover) in enumerate(zip(rows, covers)):
        if index % 64 == 0:
            check_time()
        actual = actual_image_data(source)
        caps = actual['full_caps']
        cls = cap_classes.get(caps)
        if cls is None:
            cls = len(unique_caps)
            cap_classes[caps] = cls
            unique_caps.append(caps)
        rank, dimension = actual['cross_rank'], actual['fiber_dimension']
        if sum(cover, F(0)) != dimension:
            raise AssertionError('Quiet cover and actual fiber dimension differ')
        classes.append(cls)
        ranks.append(rank)
        dimensions.append(dimension)
    check_time()
    gap_values, gap_parents, gap_labels, gap_codes = gap_geometry(rows, covers, SIZE)
    if rows == expected and sum(orbit_sizes) != comb(SIZE, arity):
        raise AssertionError('The complete labelled partition changed')
    return dict(rowsets=rows, covers=covers, outside_caps=tuple(unique_caps),
                nonquiet_caps=tuple((0,)+caps[1:] for caps in unique_caps),
                classes=np.asarray(classes, dtype=np.int16),
                quiet_caps=quiet_caps, quiet_classes=np.asarray(quiet_classes, dtype=np.int16),
                representatives=representatives, orbit_sizes=orbit_sizes,
                outside_ranks=tuple(ranks), fiber_dimensions=tuple(dimensions),
                gap_values=gap_values, gap_parents=gap_parents, gap_labels=gap_labels, gap_codes=gap_codes,
                metadata=dict(construction=CONSTRUCTION, field_order=Q, arity=arity,
                    rowset_count=len(rows), all_labelled_rowsets=rows == expected,
                    marked_outside_positions=SIZE-arity, outside_zero_equals_quiet=True,
                    actual_outside_images=True, global_flat_counts=True,
                    field_rank_subset_dp=False, labels_retained=True,
                    orbit_use='quiet covers and quiet coefficient counts only',
                    cap_class_rule='identical actual coefficient-cap vectors only',
                    full_constant=1, nonquiet_constant=0, exact_geometry=True,
                    regional_provider_evaluated=False, whole_code_certificate=False))


def prepare_all(arities=range(2, 8), *, check_time=lambda: None):
    """Complete geometry for each requested arity, with a caller-owned budget."""
    arities = tuple(arities)
    if (not arities or len(set(arities)) != len(arities)
            or any(type(a) is not int or not 2 <= a <= 7 for a in arities)):
        raise ValueError('Distinct arities2..7 required')
    return {a: prepare_geometry(a, check_time=check_time) for a in arities}

"""Exact outside-weight4 fibers for a fixed arity5, affine-rank4 set.

A nonzero affine correction with four nonzero outside positions has at
least seven global zeros. Its global zero set is therefore a hyperplane,
and the correction is a nonzero scalar multiple of a binary affine form.
Enumerate the30 nonconstant forms, retain actual feasible syndromes, and
check the complete result against exact outside-image support counts.

The input-to-outside kernel has dimension1 and is quiet. At outside scale1,
input/output fibers are x0+beta*k and y0+beta*k, beta inGF256. Joint support
counts use their exact zero breakpoints, not a floating subtraction or an
amplitude sweep. Each normalized fiber has256 words; the255 outside scales
are then counted once. A is an allowed support, not a forced input mask.

This module is exact local algebra only. No regional operator, row-set
enumeration, RS factor, probability certificate, or native work occurs.
For the CURRENT B, the30-by30 hyperplane-pairing certificate below is
everywhere nonzero. Consequently the requested weight4 branch is empty
for every affine-rank4 A5; the fiber API checks this against actual counts.
"""
from collections import Counter

from gf8_packet_maps import INV, MUL, matvec
from hd_bit_syndrome_cross_fibers import _nullspace, _solve
from hd_bit_syndrome_flat_geometry import outside_image_support_counts
from hd_bit_syndrome_outer import (
    CONSTRUCTION, affine_matrix, affine_rank, coefficient_matrix, matrix,
)

Q, SIZE = 256, 16


def _source(source):
    source = tuple(source)
    if (len(source) != 5 or source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < SIZE for i in source)
            or affine_rank(sum(1 << i for i in source)) != 4):
        raise ValueError('A sorted labelled arity5 affine-rank4 set is required')
    return source


def _support(values):
    return sum(1 << i for i, value in enumerate(values) if value)


def hyperplane_pairings():
    """Exact n^T B^-1 c for all30 nonconstant binary affine forms n,c.

    A rank4 source is contained in the unique hyperplane n^T V_A^T=0.
    A scalar binary correction c is feasible only if this pairing is0.
    All entries for the fixed B are nonzero; thus no rank4 correction has
    a global zero hyperplane. For A5 its outside weight is at least7.
    This is a finite900-entry certificate, not an enumeration of A sets.
    """
    forms = tuple(tuple(mask >> j & 1 for j in range(5)) for mask in range(2, 32))
    inverse_images = []
    for c in forms:
        value, rank = _solve(coefficient_matrix(), c, 5)
        if value is None or rank != 5:
            raise AssertionError('B must be invertible')
        inverse_images.append(value)
    pairings = tuple(tuple(matvec((n,), c)[0] for c in inverse_images) for n in forms)
    if any(value == 0 for row in pairings for value in row):
        raise AssertionError('The current rank4 no-hyperplane certificate changed')
    return pairings


def line_support_counts(input_origin, output_origin, direction):
    """Exact joint GLOBAL support counts along a shared GF256 affine line.

    The returned counts refer to the256 line parameters only. The caller
    separately accounts for nonzero scaling of the fixed outside vector.
    A nonzero direction makes these parameters distinct input vectors.
    """
    x, y, k = tuple(input_origin), tuple(output_origin), tuple(direction)
    if (any(len(values) != SIZE for values in (x, y, k)) or not any(k)
            or any(type(value) is not int or not 0 <= value < Q
                   for values in (x, y, k) for value in values)):
        raise ValueError('Two sixteen-byte origins and a nonzero shared direction required')
    exceptional = {MUL[value][INV[step]]
                   for origin in (x, y) for value, step in zip(origin, k) if step}
    parameters = [(value, 1) for value in sorted(exceptional)]
    if len(exceptional) < Q:
        generic = next(value for value in range(Q) if value not in exceptional)
        parameters.append((generic, Q-len(exceptional)))
    result = Counter()
    for beta, count in parameters:
        dx = tuple(MUL[beta][value] for value in k)
        key = (_support(tuple(a ^ d for a, d in zip(x, dx))),
               _support(tuple(a ^ d for a, d in zip(y, dx))))
        result[key] += count
    if sum(result.values()) != Q:
        raise AssertionError('Affine-line parameter partition changed')
    return dict(sorted(result.items()))


def rank4_rays(source):
    """Return all actual weight4 outside rays and their full inside fibers.

    Empty rays are a valid result. Exhaustiveness is checked at every call.
    Pattern counts include BOTH the256 fiber values and255 outside scales.
    """
    source = _source(source)
    source_mask = sum(1 << i for i in source)
    outside_mask = ((1 << SIZE)-1) ^ source_mask
    v, b, h = affine_matrix(), coefficient_matrix(), matrix()
    vat = tuple(tuple(v[i][j] for i in source) for j in range(5))
    kernel, pivots = _nullspace(vat, 5)
    if len(kernel) != 1 or len(pivots) != 4:
        raise AssertionError('The rank4 input-to-outside kernel must be one-dimensional')
    direction = tuple(kernel[0][source.index(i)] if i in source else 0 for i in range(SIZE))
    if matvec(h, direction) != direction:
        raise AssertionError('The common fiber direction must be quiet')
    rays = []
    for c_mask in range(2, 32):
        c = tuple(c_mask >> j & 1 for j in range(5))
        correction = matvec(v, c)
        outside_support = _support(correction) & outside_mask
        if outside_support.bit_count() != 4:
            continue
        syndrome, coefficient_rank = _solve(b, c, 5)
        if syndrome is None or coefficient_rank != 5:
            raise AssertionError('B must be invertible')
        local_input, source_rank = _solve(vat, syndrome, 5)
        if source_rank != 4:
            raise AssertionError('The labelled source rank changed')
        if local_input is None:
            continue
        x = tuple(local_input[source.index(i)] if i in source else 0 for i in range(SIZE))
        y = matvec(h, x)
        if y != tuple(a ^ d for a, d in zip(x, correction)) or _support(y) & outside_mask != outside_support:
            raise AssertionError('The actual H* map disagrees with the normalized fiber')
        normalized = line_support_counts(x, y, direction)
        patterns = tuple(dict(input_mask=im, input_weight=im.bit_count(),
                              output_mask=om, output_weight=om.bit_count(),
                              inside_output_mask=om & source_mask,
                              normalized_count=count, count=(Q-1)*count)
                         for (im, om), count in normalized.items())
        if any(row['input_mask'] & ~source_mask or row['output_mask'] & outside_mask != outside_support
               for row in patterns):
            raise AssertionError('A line pattern escaped its fixed input/outside support')
        rays.append(dict(affine_coefficient_mask=c_mask, affine_coefficients=c,
                         correction_values=correction, input_origin=x, output_origin=y,
                         direction=direction, outside_mask=outside_support, outside_weight=4,
                         outside_vector_count=Q-1, normalized_fiber_count=Q,
                         count=Q*(Q-1), patterns=patterns))
    expected = {mask: count for mask, count in outside_image_support_counts(source_mask).items()
                if mask.bit_count() == 4}
    found = {ray['outside_mask']: ray['outside_vector_count'] for ray in rays}
    if len(found) != len(rays) or found != expected:
        raise AssertionError('Binary-form rays do not exhaust the actual outsideweight4 counts')
    return dict(source=source, source_mask=source_mask, outside_positions_mask=outside_mask,
                rays=tuple(rays), outside_weight4_counts=expected, ray_count=len(rays),
                nonzero_vector_count=sum(ray['count'] for ray in rays),
                metadata=dict(construction=CONSTRUCTION, field_order=Q, arity=5,
                    affine_rank=4, cross_rank=4, input_to_outside_fiber_dimension=1,
                    labelled_source_fixed=True, normalized_outside_value=1,
                    scalar_multiplicity=Q-1, normalized_fiber_cardinality=Q,
                    input_support_is_allowed_not_forced=True, exhaustiveness_checked=True,
                    complete_output_masks=True, exact_local_counts=True,
                    regional_provider=False, complete_outer_kernel=False,
                    whole_code_certificate=False))


def exact_terms(source):
    """Positive {(input_weight, full_output_mask):count} on all returned fibers."""
    result = Counter()
    for ray in rank4_rays(source)['rays']:
        for row in ray['patterns']:
            result[row['input_weight'], row['output_mask']] += row['count']
    return dict(sorted(result.items()))

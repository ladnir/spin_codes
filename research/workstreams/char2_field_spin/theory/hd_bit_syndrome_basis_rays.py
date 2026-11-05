"""Exact minimum-outside-weight rays for an affine-basis input set.

For one fixed labelled A of size5 and affine rank5, the eleven outside
coordinates evaluate every affine syndrome injectively. Outside weight4
occurs exactly on five projective rays: the binary affine form that is0
at one selected A coordinate and1 at the other four. Each ray contains
255 nonzero scalar multiples. The implementation checks exhaustiveness
against the existing exact outside-image support counts, rather than
assuming the five-ray characterization.

Every outside value determines a unique input supported within A.
Small AES-field solves recover this input and its COMPLETE output mask.
There is no affine transport of H*, amplitude sampling, regional provider,
RS factor, row-set sum, or distance certificate here. Input weight need
not equal5; A is an allowed input support, not a forced exact support.
"""
from gf8_packet_maps import matvec
from hd_bit_syndrome_cross_fibers import _solve
from hd_bit_syndrome_flat_geometry import outside_image_support_counts
from hd_bit_syndrome_outer import (
    CONSTRUCTION, affine_matrix, affine_rank, coefficient_matrix, matrix,
)

Q, SIZE = 256, 16


def _source(source):
    source = tuple(source)
    if (len(source) != 5 or source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < SIZE for i in source)
            or affine_rank(sum(1 << i for i in source)) != 5):
        raise ValueError('A sorted labelled arity5 affine basis is required')
    return source


def _support(values):
    return sum(1 << i for i, value in enumerate(values) if value)


def _unique_solve(equations, target):
    values, rank = _solve(equations, target, 5)
    if values is None or rank != 5:
        raise AssertionError('The affine-basis system must have a unique solution')
    return values


def basis_rays(source):
    """Return five exact labelled projective rays, each with count255.

    Each representative has outside amplitudes in{0,1}, with four ones.
    The returned input/output vectors use all sixteen original positions.
    The isolated_source field identifies the unique A coordinate at which
    the binary affine correction vanishes. All other A corrections are1.
    """
    source = _source(source)
    source_mask = sum(1 << i for i in source)
    outside_mask = ((1 << SIZE)-1) ^ source_mask
    v, b, h = affine_matrix(), coefficient_matrix(), matrix()
    va = tuple(v[i] for i in source)
    vat = tuple(zip(*va))
    rays = []
    for isolated in source:
        c = _unique_solve(va, tuple(int(i != isolated) for i in source))
        if any(value not in (0, 1) for value in c):
            raise AssertionError('The prescribed binary affine values must give binary coefficients')
        syndrome = _unique_solve(b, c)
        local_input = _unique_solve(vat, syndrome)
        x = [0]*SIZE
        for i, value in zip(source, local_input):
            x[i] = value
        x = tuple(x)
        correction = matvec(v, c)
        y = matvec(h, x)
        if y != tuple(a ^ d for a, d in zip(x, correction)):
            raise AssertionError('The actual H* map disagrees with the syndrome solve')
        input_mask, output_mask = _support(x), _support(y)
        outside_support = output_mask & outside_mask
        if (input_mask & ~source_mask or outside_support.bit_count() != 4
                or any(y[i] not in (0, 1) for i in range(SIZE) if outside_mask >> i & 1)):
            raise AssertionError('The normalized ray has incorrect input or outside support')
        rays.append(dict(isolated_source=isolated, affine_coefficients=c,
                         correction_values=correction, input_values=x, output_values=y,
                         input_mask=input_mask, input_weight=input_mask.bit_count(),
                         outside_mask=outside_support, outside_weight=4,
                         inside_output_mask=output_mask & source_mask,
                         output_mask=output_mask, output_weight=output_mask.bit_count(),
                         count=Q-1))
    expected = {mask: count for mask, count in outside_image_support_counts(source_mask).items()
                if mask.bit_count() == 4}
    found = {ray['outside_mask']: ray['count'] for ray in rays}
    if len(found) != 5 or expected != found:
        raise AssertionError('Five rays must exhaust the exact outside-weight4 image count')
    return dict(source=source, source_mask=source_mask, outside_positions_mask=outside_mask,
                rays=tuple(rays), outside_weight4_counts=expected,
                ray_count=len(rays), nonzero_vector_count=sum(found.values()),
                metadata=dict(construction=CONSTRUCTION, field_order=Q, arity=5,
                    affine_rank=5, cross_rank=5, input_to_outside_fiber_dimension=0,
                    labelled_source_fixed=True, outside_value_normalization=1,
                    projective_ray_count_included=Q-1, exhaustiveness_checked=True,
                    input_support_is_allowed_not_forced=True,
                    complete_output_masks=True, exact_local_counts=True,
                    regional_provider=False, complete_outer_kernel=False,
                    whole_code_certificate=False))


def exact_terms(source):
    """Positive {(input_weight, full_output_mask): integer_count} on these rays."""
    result = {}
    for ray in basis_rays(source)['rays']:
        key = (ray['input_weight'], ray['output_mask'])
        result[key] = result.get(key, 0)+ray['count']
    return result

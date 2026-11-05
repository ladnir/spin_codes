"""Exact cross-map fibers for a fixed labelled input set of H*.

This module provides local linear algebra, not complete outer kernels or
a distance bound.  It supports arities1..16 algebraically.  Exhaustive
inside-support extraction has an explicit subset budget (default512).
No regional provider, probability normalization, RS floor, or row-set
multiplicity is included.

For A fixed and O its complement, C=H*[O,A] is the ACTUAL outside map.
Its image need not equal im(V_O).  The zero-outside input kernel ker C
is distinguished from the quiet kernel ker(V_A^T).  The inside output
fiber is obtained through H*, and independently constrained through H*^-1.

Fix outside output z.  For an inside allowed support S, H*^-1 implies
    H*^-1[O,S] y_S = H*^-1[O,O] z.
The sign is plus in characteristic two.  The exact count is
256^(|S|-rank H*^-1[O,S]) when consistent, and zero otherwise.
Subset inversion then counts exact inside supports, including all coset
shifts.  No affine-equivariance assumption about H* is used.
"""
from functools import lru_cache
from itertools import combinations
from math import comb

from gf8_packet_maps import INV, MUL, matvec
from hd_bit_syndrome_outer import CONSTRUCTION, affine_matrix, matrix

Q = 256


def _source(source):
    source = tuple(source)
    if (not 1 <= len(source) <= 16 or source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('A must be a sorted labelled subset of size1..16')
    return source


def _rref(rows, variables):
    """Exact AES-field elimination; appended columns are not pivot columns."""
    a = [list(row) for row in rows]
    pivots = []
    used = 0
    for column in range(variables):
        pivot = next((i for i in range(used, len(a)) if a[i][column]), None)
        if pivot is None:
            continue
        a[used], a[pivot] = a[pivot], a[used]
        inverse = INV[a[used][column]]
        a[used] = [MUL[inverse][x] for x in a[used]]
        for i, row in enumerate(a):
            if i != used and row[column]:
                scale = row[column]
                a[i] = [x ^ MUL[scale][y] for x, y in zip(row, a[used])]
        pivots.append(column)
        used += 1
    return tuple(tuple(row) for row in a), tuple(pivots)


def _nullspace(rows, variables):
    reduced, pivots = _rref(rows, variables)
    basis = []
    for free in range(variables):
        if free not in pivots:
            vector = [0] * variables
            vector[free] = 1
            for row, pivot in enumerate(pivots):
                vector[pivot] = reduced[row][free]
            basis.append(tuple(vector))
    return tuple(basis), pivots


def _solve(rows, rhs, variables):
    reduced, pivots = _rref([tuple(row) + (value,) for row, value in zip(rows, rhs)], variables)
    if any(not any(row[:variables]) and row[variables] for row in reduced):
        return None, len(pivots)
    result = [0] * variables
    for row, pivot in enumerate(pivots):
        result[pivot] = reduced[row][variables]
    return tuple(result), len(pivots)


@lru_cache(maxsize=1)
def inverse_matrix():
    h = matrix()
    augmented = [row + tuple(int(i == j) for j in range(16)) for i, row in enumerate(h)]
    reduced, pivots = _rref(augmented, 16)
    if pivots != tuple(range(16)):
        raise AssertionError('The fixed bit-syndrome matrix is singular')
    result = tuple(row[16:] for row in reduced)
    for j in range(16):
        if matvec(h, tuple(result[i][j] for i in range(16))) != tuple(int(i == j) for i in range(16)):
            raise AssertionError('Inverse construction failed')
    return result


def cross_fiber(source):
    """Return exact labelled cross image, kernels, and inside fiber bases.

    Basis vectors are rows in the returned tuples.  Outside vectors use
    the recorded O ordering; input and inside vectors use A ordering.
    The actual outside-image basis consists of independent C columns.
    """
    source = _source(source)
    outside = tuple(i for i in range(16) if i not in source)
    h, hi = matrix(), inverse_matrix()
    cross = tuple(tuple(h[i][j] for j in source) for i in outside)
    inside = tuple(tuple(h[i][j] for j in source) for i in source)
    input_kernel, pivots = _nullspace(cross, len(source))
    outside_basis = tuple(tuple(row[j] for row in cross) for j in pivots)
    inside_basis = tuple(matvec(inside, vector) for vector in input_kernel)
    inverse_cross = tuple(tuple(hi[i][j] for j in source) for i in outside)
    inverse_kernel, inverse_pivots = _nullspace(inverse_cross, len(source))
    dimension = len(source) - len(pivots)
    if (len(inverse_kernel) != dimension
            or len(_rref(inside_basis, len(source))[1]) != dimension
            or any(any(matvec(inverse_cross, vector)) for vector in inside_basis)):
        raise AssertionError('Forward and inverse inside fibers disagree')
    v = affine_matrix()
    quiet_constraints = tuple(tuple(v[i][j] for i in source) for j in range(5))
    quiet_basis, quiet_pivots = _nullspace(quiet_constraints, len(source))
    if (any(any(matvec(cross, vector)) for vector in quiet_basis)
            or any(matvec(inside, vector) != vector for vector in quiet_basis)):
        raise AssertionError('The quiet space is not fixed inside the zero-outside fiber')
    return dict(source=source, outside=outside, arity=len(source),
                cross_matrix=cross, cross_rank=len(pivots),
                outside_image_basis=outside_basis,
                outside_image_dimension=len(pivots),
                input_kernel_basis=input_kernel,
                input_kernel_dimension=dimension,
                inside_zero_fiber_basis=inside_basis,
                inside_zero_fiber_dimension=dimension,
                inverse_cross_matrix=inverse_cross,
                inverse_cross_rank=len(inverse_pivots),
                inverse_inside_kernel_basis=inverse_kernel,
                quiet_constraints=quiet_constraints, quiet_basis=quiet_basis,
                quiet_dimension=len(quiet_basis), quiet_constraint_rank=len(quiet_pivots),
                zero_outside_excess_dimension=dimension-len(quiet_basis),
                zero_outside_equals_quiet=dimension == len(quiet_basis),
                metadata=dict(construction=CONSTRUCTION, field_order=Q,
                              labels_retained=True, exact_local_linear_algebra=True,
                              complete_outer_kernel=False, regional_provider=False,
                              whole_code_certificate=False))


def _outside_values(fiber, values):
    result = (0,) * len(fiber['outside']) if values is None else tuple(values)
    if (len(result) != len(fiber['outside'])
            or any(type(x) is not int or not 0 <= x < Q for x in result)):
        raise ValueError('One AES-field byte per recorded outside coordinate required')
    return result


def outside_fiber(fiber, values=None):
    """Return a particular x_A and y_A, or None for impossible outside z."""
    values = _outside_values(fiber, values)
    x, _ = _solve(fiber['cross_matrix'], values, fiber['arity'])
    if x is None:
        return None
    h = matrix()
    y = matvec(tuple(tuple(h[i][j] for j in fiber['source']) for i in fiber['source']), x)
    return dict(outside_values=values, input_particular=x,
                inside_particular=y, inside_difference_basis=fiber['inside_zero_fiber_basis'],
                cardinality=Q ** fiber['inside_zero_fiber_dimension'])


def information_sets(fiber, *, max_candidates=12870):
    """All inside information sets for the difference space, as GLOBAL labels."""
    dimension = fiber['inside_zero_fiber_dimension']
    source = fiber['source']
    if type(max_candidates) is not int or max_candidates < 1 or comb(len(source), dimension) > max_candidates:
        raise ValueError('Information-set candidate budget exceeded')
    generator = fiber['inside_zero_fiber_basis']
    result = []
    for choice in combinations(range(len(source)), dimension):
        selected = tuple(tuple(row[i] for i in choice) for row in generator)
        if len(_rref(selected, dimension)[1]) == dimension:
            result.append(tuple(source[i] for i in choice))
    if not result:
        raise AssertionError('A finite-dimensional space needs an information set')
    return tuple(result)


def supported_count(fiber, support, values=None):
    """Exact number with inside support CONTAINED IN support, at fixed z."""
    support = tuple(support)
    if (support != tuple(sorted(set(support)))
            or any(type(i) is not int or i not in fiber['source'] for i in support)):
        raise ValueError('Inside allowed support must be a sorted subset of A')
    values = _outside_values(fiber, values)
    hi = inverse_matrix()
    equations = tuple(tuple(hi[i][j] for j in support) for i in fiber['outside'])
    rhs = matvec(tuple(tuple(hi[i][j] for j in fiber['outside']) for i in fiber['outside']), values)
    solution, r = _solve(equations, rhs, len(support))
    return 0 if solution is None else Q ** (len(support) - r)


def inside_support_counts(fiber, values=None, *, max_subsets=512):
    """Exact {GLOBAL inside-support mask:count}; explicit exhaustive budget."""
    source = fiber['source']
    if (type(max_subsets) is not int or not 1 <= max_subsets <= 65536
            or 1 << len(source) > max_subsets):
        raise ValueError('Inside-support enumeration exceeds explicit subset budget')
    values = _outside_values(fiber, values)
    if outside_fiber(fiber, values) is None:
        return {}
    table = []
    masks = []
    for local in range(1 << len(source)):
        support = tuple(source[i] for i in range(len(source)) if local >> i & 1)
        table.append(supported_count(fiber, support, values))
        masks.append(sum(1 << i for i in support))
    for bit in range(len(source)):
        step = 1 << bit
        for base in range(0, len(table), 2 * step):
            for offset in range(step):
                table[base + step + offset] -= table[base + offset]
    if min(table) < 0 or sum(table) != Q ** fiber['inside_zero_fiber_dimension']:
        raise AssertionError('Exact inside-support inversion failed')
    return {mask: count for mask, count in zip(masks, table) if count}

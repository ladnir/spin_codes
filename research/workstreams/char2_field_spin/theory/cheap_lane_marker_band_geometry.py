"""Exact vertices for an ordinary-marker band cell.

The coordinates are (k, a, b), not actual size r when holes are used.
The caller chooses that coordinate convention and supplies closed rational
intervals.  This module enumerates the REAL polytope

    klo <= k <= khi, alo <= a <= ahi, blo <= b <= bhi,
    2*a <= k, 2*b <= k, a+b <= k-1.

The last inequality removes zero-run profiles.  Nonnegative intervals and
1 <= k <= 4096 exclude the constant-support endpoints.  No source factor,
profile count, denominator, or probability is computed here.

There are nine planes and 84 triples.  Every vertex has three independent
active normals, including when the polytope has lower affine dimension.
Otherwise a nonzero vector orthogonal to all active normals permits two
opposite sufficiently small feasible displacements, contradicting extremity.
Conversely, a feasible point fixed by three independent active planes is a
vertex.  Thus exact triple enumeration gives all vertices.  Boundedness
ensures that a nonempty cell has a vertex; the empty tuple means an empty
cell.  Concave denominator logarithms may be minimized on these vertices,
but that analytic argument and directed evaluation belong to the caller.

No helper imports, floating-point operations, or work on module import.
The exported exact_tests function is run only by an authorized caller.
"""

from fractions import Fraction
from itertools import combinations, product


NORMALS = (
    (-1, 0, 0), (1, 0, 0),
    (0, -1, 0), (0, 1, 0),
    (0, 0, -1), (0, 0, 1),
    (-1, 2, 0), (-1, 0, 2), (-1, 1, 1),
)


def _rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("cell endpoints must be int or Fraction, not float")
    return Fraction(value)


def _bounds(klo, khi, alo, ahi, blo, bhi):
    endpoints = tuple(_rational(x) for x in (klo, khi, alo, ahi, blo, bhi))
    kl, kh, al, ah, bl, bh = endpoints
    if not (1 <= kl <= kh <= 4096):
        raise ValueError("require 1 <= klo <= khi <= 4096")
    if not (0 <= al <= ah <= 4096 and 0 <= bl <= bh <= 4096):
        raise ValueError("marker intervals must be ordered and within [0,4096]")
    return endpoints


def _cross(u, v):
    return (u[1]*v[2]-u[2]*v[1],
            u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0])


def vertices(klo, khi, alo, ahi, blo, bhi, guard):
    """Return a sorted tuple of distinct Fraction triples (k,a,b).

    Closed intervals may be singletons.  A well-formed but infeasible cell
    returns ().  Reversed/out-of-domain intervals and non-rational endpoint
    types raise ValueError/TypeError.  guard() is called for each triple;
    the caller owns time and memory limits.  Fractional vertices are never
    rounded to nearby integer profiles.
    """
    kl, kh, al, ah, bl, bh = _bounds(klo, khi, alo, ahi, blo, bhi)
    rhs = (-kl, kh, -al, ah, -bl, bh, Fraction(0), Fraction(0), Fraction(-1))
    found = set()
    for i, j, ell in combinations(range(9), 3):
        guard()
        row0, row1, row2 = NORMALS[i], NORMALS[j], NORMALS[ell]
        col0 = _cross(row1, row2)
        determinant = row0[0]*col0[0] + row0[1]*col0[1] + row0[2]*col0[2]
        if not determinant:
            continue
        col1, col2 = _cross(row2, row0), _cross(row0, row1)
        d0, d1, d2 = rhs[i], rhs[j], rhs[ell]
        k = (d0*col0[0] + d1*col1[0] + d2*col2[0]) / determinant
        a = (d0*col0[1] + d1*col1[1] + d2*col2[1]) / determinant
        b = (d0*col0[2] + d1*col1[2] + d2*col2[2]) / determinant
        # These are exactly all nine inequalities, without a tolerance.
        if (kl <= k <= kh and al <= a <= ah and bl <= b <= bh
                and 2*a <= k and 2*b <= k and a+b <= k-1):
            found.add((k, a, b))
    return tuple(sorted(found))


def _grid_points(bounds, guard):
    """Independent half-integer feasible-grid enumeration for tiny tests."""
    kl, kh, al, ah, bl, bh = bounds
    result = []
    for twice_k in range(2*kl, 2*kh+1):
        guard()
        k = Fraction(twice_k, 2)
        for twice_a in range(2*al, 2*ah+1):
            a = Fraction(twice_a, 2)
            for twice_b in range(2*bl, 2*bh+1):
                b = Fraction(twice_b, 2)
                if a <= k/2 and b <= k/2 and a+b+1 <= k:
                    result.append((k, a, b))
    return tuple(result)


def _dot(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]


def exact_tests(guard):
    """Run small exact geometry controls; return counts, not proof evidence.

    The finite tests supplement the active-normal proof above.  They do
    not establish a source-size bound or replace the caller's replay.
    """
    checks = 0

    def check(condition):
        nonlocal checks
        assert condition
        checks += 1

    def run(bounds):
        return vertices(*bounds, guard)

    q = Fraction
    check(run((1, 1, 0, 2, 0, 2)) == ((q(1), q(0), q(0)),))
    check(run((2, 2, 1, 1, 1, 1)) == ())
    check(run((1, 4, 0, 0, 0, 0)) == ((q(1), q(0), q(0)), (q(4), q(0), q(0))))
    check(run((4, 4, 1, 1, 1, 1)) == ((q(4), q(1), q(1)),))
    check(run((2, 5, 1, 1, 0, 0)) == ((q(2), q(1), q(0)), (q(5), q(1), q(0))))
    # The fixed odd-size section has genuine half-integer vertices.
    expected = ((q(3), q(0), q(0)), (q(3), q(0), q(3, 2)),
                (q(3), q(1, 2), q(3, 2)), (q(3), q(3, 2), q(0)),
                (q(3), q(3, 2), q(1, 2)))
    check(run((3, 3, 0, 2, 0, 2)) == expected)
    rational_bounds = (q(3, 2), q(5, 2), 0, 1, 0, 1)
    rational_vertices = run(rational_bounds)
    check((q(3, 2), q(1, 2), q(0)) in rational_vertices)
    check((q(5, 2), q(1), q(1, 2)) in rational_vertices)
    check((q(2), q(1), q(0)) in rational_vertices)
    for point in rational_vertices:
        check(all(isinstance(x, Fraction) for x in point))

    # For these integer boxes, half-grid enumeration includes every vertex.
    # Compare extrema of every {-1,0,1} linear objective against the FULL
    # feasible half-grid, independently of triple intersections.
    directions = tuple(x for x in product((-1, 0, 1), repeat=3) if x != (0, 0, 0))
    k_intervals = ((1, 1), (2, 2), (3, 3), (1, 3), (2, 4), (3, 6))
    intervals = ((0, 0), (0, 3), (1, 1), (2, 3))
    cells = grid_points = 0
    for kbox, abox, bbox in product(k_intervals, intervals, intervals):
        guard()
        box = kbox + abox + bbox
        actual, grid = run(box), _grid_points(box, guard)
        cells += 1
        grid_points += len(grid)
        check(bool(actual) == bool(grid))
        check(set(actual).issubset(set(grid)))
        swapped = run(kbox + bbox + abox)
        check(tuple(sorted((k, b, a) for k, a, b in actual)) == swapped)
        if not grid:
            continue
        for direction in directions:
            guard()
            check(max(_dot(direction, x) for x in actual)
                  == max(_dot(direction, x) for x in grid))
        # Every reported vertex has an independent active triple.
        kl, kh, al, ah, bl, bh = box
        rhs = (-kl, kh, -al, ah, -bl, bh, 0, 0, -1)
        for point in actual:
            active = [NORMALS[i] for i in range(9) if _dot(NORMALS[i], point) == rhs[i]]
            check(any(_dot(a, _cross(b, c)) != 0 for a, b, c in combinations(active, 3)))

    for bad in ((0, 1, 0, 1, 0, 1), (1, 4097, 0, 1, 0, 1),
                (2, 1, 0, 1, 0, 1), (1, 2, -1, 1, 0, 1),
                (1, 2, 2, 1, 0, 1), (1, 2, 0, 1, 2, 1)):
        try:
            run(bad)
        except ValueError:
            checks += 1
        else:
            raise AssertionError("invalid interval was accepted")
    for bad_value in (1.0, True, "1"):
        try:
            run((bad_value, 2, 0, 1, 0, 1))
        except TypeError:
            checks += 1
        else:
            raise AssertionError("non-exact endpoint type was accepted")
    return {"checks": checks, "grid_cells": cells, "feasible_grid_points": grid_points,
            "linear_objectives": len(directions), "planes": 9, "triples_per_cell": 84,
            "numerical_replay": False, "source_size_certificate": False}

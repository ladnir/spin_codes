"""Exact vertices for size-normalized ordinary-marker cells.

For k in a closed positive band, define x=2048*a/k and y=2048*b/k.
The exported vertices function returns the real polytope with x and y in
specified closed rational intervals, 2*a<=k, 2*b<=k, and a+b<=k-1.
It returns coordinates (k,a,b), which the existing denominator evaluator
uses.  It neither transports marker values nor changes the source factor.

Reference integer cells [lo,hi] must first become CONTINUOUS cells: use
cuts hi+1/2 between neighboring cells and outer endpoints0 and1024.
These intervals cover every normalized coordinate, including those from
odd integer sizes.  Closed boundaries can overlap; a caller may assign a
boundary deterministically or maximize over all cells.  No union factor
for cells is required when taking a maximum.  The old unexpanded integer
intervals would leave gaps after normalization and must not be used.

All nine constraints are affine in (k,a,b).  At a vertex their active
normals span three dimensions, including for lower-dimensional cells:
otherwise a sufficiently small displacement in either direction along
their common nullspace remains feasible.  Conversely, a feasible point
fixed by three independent active planes is a vertex.  Exact enumeration
of all84 triples therefore returns all vertices.  The set is bounded, so
an empty result means an empty cell.  Concave-denominator evaluation is
the caller's responsibility; this module produces no numerical bound.

There is no main, proof-provider import, or numerical work on import.
"""

from fractions import Fraction
from itertools import combinations, product


REFERENCE_WEIGHT = 2048
MAX_NORMALIZED_FULL = 1024


def _exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("endpoints must be int or Fraction, not float")
    return Fraction(value)


def _bounds(klo, khi, xlo, xhi, ylo, yhi):
    values = tuple(_exact(v) for v in (klo, khi, xlo, xhi, ylo, yhi))
    kl, kh, xl, xh, yl, yh = values
    if not 1 <= kl <= kh <= 4096:
        raise ValueError("require 1 <= klo <= khi <= 4096")
    if not (0 <= xl <= xh <= 1024 and 0 <= yl <= yh <= 1024):
        raise ValueError("normalized intervals must be ordered within [0,1024]")
    return values


def _cross(u, v):
    return (u[1]*v[2]-u[2]*v[1],
            u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0])


def _dot(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]


def vertices(klo, khi, xlo, xhi, ylo, yhi, guard):
    """Return a sorted tuple of distinct exact Fraction triples (k,a,b).

    Singleton bands and interval endpoints are allowed.  A well-formed
    infeasible cell returns ().  Invalid types or interval domains raise
    TypeError or ValueError.  The caller supplies guard() for every triple
    and owns runtime/memory limits.  No vertex is rounded to an integer.
    """
    kl, kh, xl, xh, yl, yh = _bounds(klo, khi, xlo, xhi, ylo, yhi)
    normals = ((-1, 0, 0), (1, 0, 0),
               (xl, -2048, 0), (-xh, 2048, 0),
               (yl, 0, -2048), (-yh, 0, 2048),
               (-1, 2, 0), (-1, 0, 2), (-1, 1, 1))
    zero = Fraction(0)
    rhs = (-kl, kh, zero, zero, zero, zero, zero, zero, Fraction(-1))
    found = set()
    for i, j, ell in combinations(range(9), 3):
        guard()
        row0, row1, row2 = normals[i], normals[j], normals[ell]
        col0 = _cross(row1, row2)
        determinant = _dot(row0, col0)
        if not determinant:
            continue
        col1, col2 = _cross(row2, row0), _cross(row0, row1)
        d0, d1, d2 = rhs[i], rhs[j], rhs[ell]
        k = (d0*col0[0]+d1*col1[0]+d2*col2[0])/determinant
        a = (d0*col0[1]+d1*col1[1]+d2*col2[1])/determinant
        b = (d0*col0[2]+d1*col1[2]+d2*col2[2])/determinant
        if (kl <= k <= kh and xl*k <= 2048*a <= xh*k
                and yl*k <= 2048*b <= yh*k
                and 2*a <= k and 2*b <= k and a+b <= k-1):
            found.add((k, a, b))
    return tuple(sorted(found))


def _section(k, xl, xh, yl, yh):
    """Independent rectangle clipped by a+b<=k-1, for exact tests."""
    al, ah, bl, bh = (k*x/2048 for x in (xl, xh, yl, yh))
    limit = k-1
    candidates = [(a, b) for a in (al, ah) for b in (bl, bh)]
    candidates += [(a, limit-a) for a in (al, ah)]
    candidates += [(limit-b, b) for b in (bl, bh)]
    return {(k, a, b) for a, b in candidates
            if al <= a <= ah and bl <= b <= bh and a+b <= limit}


def _reference(bounds, guard):
    """Independent frustum-edge/cutting-plane enumeration for tests.

    Before the run cut, the polytope is the frustum over a rectangle of
    ray directions.  Its edges lie in the end sections or follow one of
    the four corner rays.  Clipping adds only intersections with these
    edges.  This test reference does not solve plane triples.
    """
    kl, kh, xl, xh, yl, yh = map(Fraction, bounds)
    result = _section(kl, xl, xh, yl, yh) | _section(kh, xl, xh, yl, yh)
    for x, y in product((xl, xh), (yl, yh)):
        guard()
        denominator = 2048-x-y
        if denominator <= 0:
            continue
        k = Fraction(2048)/denominator
        if kl <= k <= kh:
            result.add((k, k*x/2048, k*y/2048))
    return tuple(sorted(result))


def exact_tests(guard):
    """Return tiny exact-test counts; a caller must authorize execution."""
    checks = 0

    def check(value):
        nonlocal checks
        assert value
        checks += 1

    def run(bounds):
        return vertices(*bounds, guard)

    q = Fraction
    check(run((1, 1, 0, 1024, 0, 1024)) == ((q(1), q(0), q(0)),))
    check(run((1, 4096, 1024, 1024, 1024, 1024)) == ())
    check(run((1, 5, 0, 0, 0, 0)) == ((q(1), q(0), q(0)), (q(5), q(0), q(0))))
    # The run plane cuts a corner ray at a nonintegral size, with no rounding.
    check(run((1, 4, 256, 256, 256, 256))
          == ((q(4, 3), q(1, 6), q(1, 6)), (q(4), q(1, 2), q(1, 2))))
    check(run((3, 3, 0, 1024, 0, 1024))
          == ((q(3), q(0), q(0)), (q(3), q(0), q(3, 2)),
              (q(3), q(1, 2), q(3, 2)), (q(3), q(3, 2), q(0)),
              (q(3), q(3, 2), q(1, 2))))

    kboxes = ((1, 1), (3, 3), (1, 4), (q(3, 2), q(7, 2)), (2, 7))
    intervals = ((0, 0), (0, 1024), (256, 768), (q(2045, 2), 1024))
    directions = tuple(x for x in product((-1, 0, 1), repeat=3) if any(x))
    cells = sampled_points = 0
    for kbox, xbox, ybox in product(kboxes, intervals, intervals):
        guard()
        box = kbox+xbox+ybox
        actual, reference = run(box), _reference(box, guard)
        cells += 1
        check(actual == reference)
        swapped = run(kbox+ybox+xbox)
        check(tuple(sorted((k, b, a) for k, a, b in actual)) == swapped)
        for point in actual:
            check(all(isinstance(v, Fraction) for v in point))
        if not actual:
            continue
        # Every exact grid point satisfies all linear objectives' vertex
        # bounds.  The independent edge reference checks all extrema too.
        kl, kh, xl, xh, yl, yh = map(q, box)
        samples = []
        for step, sx, sy in product(range(5), range(5), range(5)):
            k = kl+(kh-kl)*step/4
            x, y = xl+(xh-xl)*sx/4, yl+(yh-yl)*sy/4
            point = (k, k*x/2048, k*y/2048)
            if point[1]+point[2] <= k-1:
                samples.append(point)
        sampled_points += len(samples)
        for direction in directions:
            guard()
            maximum = max(_dot(direction, p) for p in actual)
            check(maximum == max(_dot(direction, p) for p in reference))
            check(all(_dot(direction, p) <= maximum for p in samples))

    # Continuous half-integer cuts cover normalized targets of every parity.
    # These are a small reference partition, not a replacement marker bank.
    cuts = (q(0), q(1, 2), q(1023, 2), q(2047, 2), q(1024))
    covered_targets = 0
    odd_targets = 0
    for k in (*range(1, 18), 2047, 2048, 2049, 4096):
        guard()
        for a in range(k//2+1):
            x = q(2048*a, k)
            check(any(lo <= x <= hi for lo, hi in zip(cuts, cuts[1:])))
            covered_targets += 1
            odd_targets += k % 2
    # Original integer labels are preserved by hi+1/2 cuts at k=2048.
    integer_cells = ((0, 0), (1, 511), (512, 1023), (1024, 1024))
    for index, (lo, hi) in enumerate(integer_cells):
        guard()
        for a in range(lo, hi+1):
            check(cuts[index] <= a <= cuts[index+1])
    # The original offender no longer belongs to the almost-full marker.
    check(q(2048*1023, 2304) < q(2045, 2))
    check(q(2304)-2*q(2304)*q(2045, 2)/2048 == q(27, 8))

    for bad in ((0, 1, 0, 1, 0, 1), (1, 4097, 0, 1, 0, 1),
                (2, 1, 0, 1, 0, 1), (1, 2, -1, 1, 0, 1),
                (1, 2, 0, 1025, 0, 1), (1, 2, 2, 1, 0, 1)):
        try:
            run(bad)
        except ValueError:
            checks += 1
        else:
            raise AssertionError("invalid normalized cell accepted")
    for value in (1.0, True, "1"):
        try:
            run((value, 2, 0, 1, 0, 1))
        except TypeError:
            checks += 1
        else:
            raise AssertionError("non-exact endpoint accepted")
    return {"checks": checks, "cells": cells, "sampled_feasible_points": sampled_points,
            "covered_integer_marginals": covered_targets, "odd_weight_marginals": odd_targets,
            "linear_objectives": len(directions), "planes": 9, "triples_per_cell": 84,
            "reference_weight": 2048, "fractional_vertices_preserved": True,
            "source_size_certificate": False}

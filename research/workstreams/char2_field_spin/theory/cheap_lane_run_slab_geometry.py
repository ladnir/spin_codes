"""Exact normalized marker-cell geometry with an inclusive run slab.

The coordinate size is k, with normalized matching counts x=2048*a/k
and y=2048*b/k.  The caller supplies closed intervals for k,x,y and the
run count c=k-a-b.  The original positive-run constraint c>=1 remains.
All returned vertices are exact (k,a,b) Fractions, not rounded supports.

The two run planes augment the frozen nine-plane cell to eleven planes.
Every vertex has three independent active normals, even in a singleton
size band or run slice: otherwise a small displacement in either common
nullspace direction remains feasible.  Conversely three independent
active planes fix a feasible vertex.  Enumerating all165 triples and
checking every inequality therefore suffices.  The cell is bounded; an
empty vertex tuple means the well-formed cell is empty.

Integer run slabs [1,w],[w+1,2*w],... partition positive integer runs.
Their closed REAL polytopes can have gaps between neighboring slabs;
this does not omit an integer target.  Counts must be intersected with
each slab separately.  This module supplies geometry only, with no
source factor, route bound, support count, or numerical certificate.

No provider imports, main, or work on import.  All execution, including
exact_tests, belongs to the caller and uses its resource guard.
"""

from fractions import Fraction as F
from itertools import combinations, product


REFERENCE_WEIGHT = 2048
PLANES = 11
TRIPLES = 165


def _exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise TypeError("endpoints must be int or Fraction, not float")
    return F(value)


def _bounds(klo, khi, leftlo, lefthi, rightlo, righthi, clo, chi):
    values = tuple(_exact(value) for value in
                   (klo, khi, leftlo, lefthi, rightlo, righthi, clo, chi))
    kl, kh, xl, xh, yl, yh, cl, ch = values
    if not 1 <= kl <= kh <= 4096:
        raise ValueError("require1<=klo<=khi<=4096")
    if not (0 <= xl <= xh <= 1024 and 0 <= yl <= yh <= 1024):
        raise ValueError("normalized intervals must be ordered within[0,1024]")
    if not 1 <= cl <= ch <= 4096:
        raise ValueError("require1<=clo<=chi<=4096")
    return values


def _cross(u, v):
    return (u[1]*v[2]-u[2]*v[1],
            u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0])


def _dot(u, v):
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]


def _inside(point, bounds):
    k, a, b = point
    kl, kh, xl, xh, yl, yh, cl, ch = bounds
    return (kl <= k <= kh and xl*k <= 2048*a <= xh*k
            and yl*k <= 2048*b <= yh*k and 2*a <= k and 2*b <= k
            and a+b <= k-1 and cl <= k-a-b <= ch)


def vertices(klo, khi, leftlo, lefthi, rightlo, righthi, clo, chi, guard):
    """Return all distinct exact vertices, sorted by(k,a,b).

    Singleton or fractional intervals are allowed.  Infeasible cells
    return(); invalid types/domains raise.  guard() runs for each of the
    165 triples.  The caller owns all time and memory bounds.
    """
    bounds = _bounds(klo, khi, leftlo, lefthi, rightlo, righthi, clo, chi)
    kl, kh, xl, xh, yl, yh, cl, ch = bounds
    normals = ((-1, 0, 0), (1, 0, 0),
               (xl, -2048, 0), (-xh, 2048, 0),
               (yl, 0, -2048), (-yh, 0, 2048),
               (-1, 2, 0), (-1, 0, 2), (-1, 1, 1),
               (-1, 1, 1), (1, -1, -1))
    zero = F(0)
    rhs = (-kl, kh, zero, zero, zero, zero, zero, zero, F(-1), -cl, ch)
    found = set()
    for i, j, ell in combinations(range(PLANES), 3):
        guard()
        row0, row1, row2 = normals[i], normals[j], normals[ell]
        col0 = _cross(row1, row2)
        determinant = _dot(row0, col0)
        if not determinant:
            continue
        col1, col2 = _cross(row2, row0), _cross(row0, row1)
        d0, d1, d2 = rhs[i], rhs[j], rhs[ell]
        point = ((d0*col0[0]+d1*col1[0]+d2*col2[0])/determinant,
                 (d0*col0[1]+d1*col1[1]+d2*col2[1])/determinant,
                 (d0*col0[2]+d1*col1[2]+d2*col2[2])/determinant)
        if _inside(point, bounds):
            found.add(point)
    return tuple(sorted(found))


def _section(k, xl, xh, yl, yh, cl, ch):
    """Independent fixed-k rectangle intersected with two sum halfplanes."""
    al, ah, bl, bh = (k*x/2048 for x in (xl, xh, yl, yh))
    lower, upper = k-ch, k-cl
    candidates = [(a, b) for a in (al, ah) for b in (bl, bh)]
    for limit in (lower, upper):
        candidates.extend((a, limit-a) for a in (al, ah))
        candidates.extend((limit-b, b) for b in (bl, bh))
    return {(k, a, b) for a, b in candidates
            if al <= a <= ah and bl <= b <= bh and lower <= a+b <= upper}


def _reference(bounds, guard):
    """Independent frustum clipping reference; no plane-triple solutions.

    The normalized rectangle forms a frustum between two k sections.
    Its remaining edges are four corner rays.  Parallel run-plane cuts
    add intersections on those edges or the section edges; they cannot
    intersect each other in an additional edge.  Coincident run planes
    give the same section vertices and still need no special rounding.
    """
    kl, kh, xl, xh, yl, yh, cl, ch = map(F, bounds)
    found = (_section(kl, xl, xh, yl, yh, cl, ch)
             | _section(kh, xl, xh, yl, yh, cl, ch))
    for x, y, c in product((xl, xh), (yl, yh), (cl, ch)):
        guard()
        denominator = 2048-x-y
        if denominator <= 0:
            continue
        k = 2048*c/denominator
        point = (k, k*x/2048, k*y/2048)
        if _inside(point, (kl, kh, xl, xh, yl, yh, cl, ch)):
            found.add(point)
    return tuple(sorted(found))


def exact_tests(base_geometry, guard):
    """Small exact controls only; the parent authorizes their execution.

    base_geometry is the frozen normalized geometry module, supplied by
    the caller.  No provider is loaded here.  Full-range tests compare
    linear extrema, without requiring identical redundant vertex lists.
    """
    checks = 0

    def check(condition):
        nonlocal checks
        assert condition
        checks += 1

    def run(bounds):
        return vertices(*bounds, guard)

    check(run((1, 1, 0, 1024, 0, 1024, 1, 1)) == ((F(1), F(0), F(0)),))
    check(run((1, 5, 0, 0, 0, 0, 2, 3))
          == ((F(2), F(0), F(0)), (F(3), F(0), F(0))))
    check(run((3, 3, 0, 1024, 0, 1024, 1, 1))
          == ((F(3), F(1, 2), F(3, 2)), (F(3), F(3, 2), F(1, 2))))
    check(run((3, 3, 0, 1024, 0, 1024, 3, 3)) == ((F(3), F(0), F(0)),))
    check(run((1, 4, 256, 256, 256, 256, 1, F(5, 4)))
          == ((F(4, 3), F(1, 6), F(1, 6)),
              (F(5, 3), F(5, 24), F(5, 24))))
    check(run((3, 3, 0, 1024, 0, 1024, 4, 4)) == ())
    check(run((1, 7, 1024, 1024, 1024, 1024, 1, 7)) == ())

    directions = tuple(d for d in product((-1, 0, 1), repeat=3) if any(d))
    kboxes = ((1, 1), (3, 3), (1, 4), (F(3, 2), F(7, 2)), (2, 7))
    intervals = ((0, 0), (0, 1024), (256, 768))
    slabs = ((1, 1), (1, 3), (2, 4), (F(3, 2), F(5, 2)))
    cells = sampled = full_range = 0
    for kbox, xbox, ybox in product(kboxes, intervals, intervals):
        guard()
        base_bounds = kbox+xbox+ybox
        baseline = base_geometry.vertices(*base_bounds, guard)
        broad = run(base_bounds+(1, 4096))
        check(bool(baseline) == bool(broad)); full_range += 1
        for direction in directions:
            guard()
            if baseline:
                check(max(_dot(direction, p) for p in baseline)
                      == max(_dot(direction, p) for p in broad))
        for slab in slabs:
            bounds = base_bounds+slab
            actual, reference = run(bounds), _reference(bounds, guard)
            cells += 1; check(actual == reference)
            swapped = run(kbox+ybox+xbox+slab)
            check(tuple(sorted((k, b, a) for k, a, b in actual)) == swapped)
            for point in actual:
                check(_inside(point, tuple(map(F, bounds))))
                check(all(isinstance(value, F) for value in point))
            kl, kh, xl, xh, yl, yh, cl, ch = map(F, bounds)
            samples = []
            for step, sx, sy in product(range(4), repeat=3):
                guard()
                k = kl+(kh-kl)*step/3
                x, y = xl+(xh-xl)*sx/3, yl+(yh-yl)*sy/3
                point = (k, k*x/2048, k*y/2048)
                if _inside(point, (kl, kh, xl, xh, yl, yh, cl, ch)):
                    samples.append(point)
            sampled += len(samples)
            check(not samples or bool(actual))
            if actual:
                for direction in directions:
                    guard()
                    maximum = max(_dot(direction, p) for p in actual)
                    check(maximum == max(_dot(direction, p) for p in reference))
                    check(all(_dot(direction, p) <= maximum for p in samples))

    # Exhaust every small integer(k,a,b).  Slabs partition integer runs,
    # not the full real interval between consecutive integer endpoints.
    integer_targets = odd_targets = partitions = 0
    for k in range(1, 18):
        for width in (1, 2, 3, 5):
            guard()
            partition = tuple((lo, min(k, lo+width-1)) for lo in range(1, k+1, width))
            check([c for lo, hi in partition for c in range(lo, hi+1)] == list(range(1, k+1)))
            partitions += 1
            for a, b in product(range(k//2+1), repeat=2):
                guard(); c = k-a-b
                if c < 1:
                    continue
                matching = [(lo, hi) for lo, hi in partition if lo <= c <= hi]
                check(len(matching) == 1)
                lo, hi = matching[0]
                x, y = F(2048*a, k), F(2048*b, k)
                bounds = tuple(map(F, (k, k, x, x, y, y, lo, hi)))
                check(_inside((F(k), F(a), F(b)), bounds))
                integer_targets += 1; odd_targets += k % 2
    # Two integer-separated slabs intentionally need not cover k=5/2
    # on the a=b=0 ray; integer sizes1..4 remain covered exactly once.
    check(not _inside((F(5, 2), F(0), F(0)), tuple(map(F, (1, 4, 0, 0, 0, 0, 1, 2)))))
    check(not _inside((F(5, 2), F(0), F(0)), tuple(map(F, (1, 4, 0, 0, 0, 0, 3, 4)))))

    invalid = ((0, 1, 0, 1, 0, 1, 1, 1), (1, 4097, 0, 1, 0, 1, 1, 2),
               (2, 1, 0, 1, 0, 1, 1, 2), (1, 2, -1, 1, 0, 1, 1, 2),
               (1, 2, 0, 1025, 0, 1, 1, 2), (1, 2, 2, 1, 0, 1, 1, 2),
               (1, 2, 0, 1, 0, 1, 0, 1), (1, 2, 0, 1, 0, 1, 2, 1),
               (1, 2, 0, 1, 0, 1, 1, 4097))
    for bounds in invalid:
        try: run(bounds)
        except ValueError: checks += 1
        else: raise AssertionError("invalid run-slab bounds accepted")
    for value in (1.0, True, "1"):
        for index in (0, 2, 6):
            bounds = [1, 2, 0, 1024, 0, 1024, 1, 2]; bounds[index] = value
            try: run(bounds)
            except TypeError: checks += 1
            else: raise AssertionError("non-exact run-slab endpoint accepted")
    return {"checks": checks, "cells": cells, "full_run_range_cells": full_range,
            "sampled_feasible_points": sampled, "linear_objectives": len(directions),
            "integer_run_partitions": partitions, "covered_integer_targets": integer_targets,
            "odd_coordinate_targets": odd_targets, "planes": PLANES,
            "triples_per_cell": TRIPLES, "fractional_vertices_preserved": True,
            "singleton_size_and_run_controls": True, "source_size_certificate": False}

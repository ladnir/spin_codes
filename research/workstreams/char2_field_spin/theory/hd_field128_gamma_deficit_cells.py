"""Exact parity/deficit interior cells and projected support-count bounds.

Fix eps in {0,1}. Coordinates (u,d,e) mean the profile
(r,a,b)=(2u+eps,u-d,u-e). Even interiors require d,e>=1;
odd interiors require d,e>=0. Singleton-free even faces remain separate
objects of the frozen Cell type; this module neither represents nor counts
them as deficit interiors. No old numeric evaluator is modified or adapted.

Counting makes one bundle query per feasible u and no profile-pair loop.
The inherited bounded PrefixCache may build O(pairs)-length tables for a
cold weight; streaming many cold weights retains that per-weight cost.
"""

import ast
from dataclasses import dataclass
from fractions import Fraction as F
from functools import cached_property
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from hd_field128_gamma_cell_geometry import _exact_point, _intersection
from hd_field128_gamma_cell_counts import PrefixCache


@dataclass(frozen=True)
class DeficitCell:
    """Fixed-parity native domain intersected with an inclusive integer box."""

    eps: int
    ranges: tuple
    pairs: int = 4096

    kind = 'deficit_interior'
    dimension = 3
    coordinate_names = ('u', 'd', 'e')

    def __post_init__(self):
        if (type(self.eps) is not int or self.eps not in (0, 1)
                or type(self.pairs) is not int or self.pairs < 1):
            raise ValueError('Expected parity0/1 and a positive integer pair count')
        try:
            ranges = tuple(tuple(row) for row in self.ranges)
        except TypeError as error:
            raise ValueError('Expected three inclusive integer ranges') from error
        if (len(ranges) != 3 or any(len(row) != 2
                or any(type(value) is not int for value in row) or row[0] > row[1] for row in ranges)):
            raise ValueError('Expected three ordered inclusive integer ranges')
        object.__setattr__(self, 'ranges', ranges)

    @cached_property
    def constraints(self):
        """Exact inequalities (coefficient_tuple,rhs), meaning dot<=rhs."""
        lower, p, eps = 1-self.eps, self.pairs, self.eps
        rows = [((-1, 0, 0), -lower), ((1, 0, 0), p-1),
                ((0, -1, 0), -lower), ((0, 0, -1), -lower),
                ((-1, 1, 0), 0), ((-1, 0, 1), 0),
                ((1, 1, 0), p-eps), ((1, 0, 1), p-eps)]
        for axis, (low, high) in enumerate(self.ranges):
            normal = tuple(int(axis == index) for index in range(3))
            rows.append((normal, high))
            rows.append((tuple(-value for value in normal), -low))
        tight = {}
        for normal, rhs in rows:
            tight[normal] = min(rhs, tight.get(normal, rhs))
        return tuple(sorted(tight.items()))

    def contains_real(self, point):
        point = _exact_point(point, 3)
        return all(sum(a*x for a, x in zip(normal, point)) <= rhs for normal, rhs in self.constraints)

    def contains(self, point):
        point = _exact_point(point, 3)
        return all(F(value).denominator == 1 for value in point) and self.contains_real(point)

    def profile(self, point):
        """Affine profile map; this operation alone does not assert membership."""
        u, d, e = _exact_point(point, 3)
        return 2*u+self.eps, u-d, u-e

    @cached_property
    def _vertices(self):
        vertices = set()
        for rows in combinations(self.constraints, 3):
            intersection = _intersection(rows, 3)
            if intersection is None:
                continue
            numerators, denominator = intersection
            if all(sum(a*x for a, x in zip(normal, numerators)) <= rhs*denominator
                   for normal, rhs in self.constraints):
                vertices.add(tuple(F(value, denominator) for value in numerators))
        return tuple(sorted(vertices))

    def vertices(self):
        """Complete rational vertices, including lower-dimensional cells.

        Every vertex of this bounded polytope has three independent active
        normals, even when the polytope is a segment or point. Empty output
        therefore means real emptiness, not merely absence of integer vertices.
        """
        return self._vertices

    def real_empty(self):
        return not self.vertices()

    def integer_point(self):
        """O(1) integer-feasibility witness using the least allowed deficits."""
        (ul, uh), (dl, dh), (el, eh) = self.ranges
        lower = 1-self.eps
        dl, el = max(lower, dl), max(lower, el)
        ul = max(lower, ul, dl, el)
        uh = min(self.pairs-1, uh, self.pairs-self.eps-dl, self.pairs-self.eps-el)
        return None if dl > dh or el > eh or ul > uh else (ul, dl, el)

    def integer_empty(self):
        return self.integer_point() is None

    def split(self, axis, integer_cut):
        """Disjoint lattice children <=cut and >=cut+1, not a real partition."""
        if type(axis) is not int or not 0 <= axis < 3 or type(integer_cut) is not int:
            raise ValueError('Expected a valid axis and an integer cut')
        low, high = self.ranges[axis]
        if not low <= integer_cut < high:
            raise ValueError('Cut must give two nonempty integer box ranges')
        left, right = list(self.ranges), list(self.ranges)
        left[axis], right[axis] = (low, integer_cut), (integer_cut+1, high)
        return DeficitCell(self.eps, tuple(left), self.pairs), DeficitCell(self.eps, tuple(right), self.pairs)

    def histogram_bounds(self, parity, *, integer=True):
        """Extrema for (single,full); integer envelopes contain all lattice images."""
        if type(parity) is not int or parity not in (0, 1) or type(integer) is not bool:
            raise ValueError('Expected parity0/1 and a boolean envelope flag')
        vertices = self.vertices()
        if not vertices:
            return None
        singles = tuple(2*point[parity+1]+self.eps for point in vertices)
        full = tuple(point[0]-point[parity+1] for point in vertices)
        bounds = ((min(singles), max(singles)), (min(full), max(full)))
        if not integer:
            return bounds
        return tuple((-((-low.numerator)//low.denominator), high.numerator//high.denominator)
                     for low, high in bounds)

    def descriptor(self):
        return dict(kind=self.kind, eps=self.eps, ranges=[list(row) for row in self.ranges], pairs=self.pairs)

    @classmethod
    def from_descriptor(cls, descriptor):
        if (not isinstance(descriptor, dict) or set(descriptor) != {'kind', 'eps', 'ranges', 'pairs'}
                or descriptor['kind'] != cls.kind):
            raise ValueError('Expected a complete deficit-interior descriptor')
        return cls(descriptor['eps'], descriptor['ranges'], descriptor['pairs'])


_CACHE = PrefixCache()


def _count(cell, budget, cache):
    if not isinstance(cell, DeficitCell):
        raise ValueError('Expected an explicit DeficitCell, not a legacy Cell')
    cache = _CACHE if cache is None else cache
    if not isinstance(cache, PrefixCache):
        raise ValueError('Expected the frozen bounded PrefixCache')
    budget.check()
    (ul, uh), (dl, dh), (el, eh) = cell.ranges
    p, eps, lower = cell.pairs, cell.eps, 1-cell.eps
    dl, el = max(lower, dl), max(lower, el)
    ul = max(lower, ul, dl, el)
    uh = min(p-1, uh, p-eps-dl, p-eps-el)
    even, odd, runs, slices = 0, 0, 0, 0
    if dl <= dh and el <= eh:
        for u in range(ul, uh+1):
            budget.check()
            maximum = min(u, p-u-eps)
            dhi, ehi = min(dh, maximum), min(eh, maximum)
            if dl > dhi or el > ehi:
                continue
            weight = 2*u+eps
            marginal, cycle = cache.get('interior', p, weight, budget)
            even += marginal.interval(u-dhi, u-dl)
            odd += marginal.interval(u-ehi, u-el)
            runs += cycle.interval(dl+el+eps, dhi+ehi+eps)
            slices += 1
    budget.check()
    return min(even, odd, runs), (even, odd, runs), slices


def support_count(cell, budget, *, cache=None):
    """Exact integer upper bound from projection sums, not an exact joint count."""
    return _count(cell, budget, cache)[0]


def count_diagnostics(cell, budget, *, cache=None):
    count, projections, slices = _count(cell, budget, cache)
    return dict(support_count_upper=count, exact_support_count=False,
                method='parity-deficit-minimum-of-exact-projection-sums',
                projection_sums=dict(zip(('even_marginal', 'odd_marginal', 'cycle_runs'), projections)),
                weight_parity=cell.eps, nonempty_integer_slices=slices,
                profile_pairs_enumerated=0, omitted_coordinate_factor=False,
                singleton_free_faces_included=False)


def source_pins():
    directory, pending, pins = Path(__file__).resolve().parent, [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in pins:
            continue
        raw = (directory/name).read_bytes()
        pins[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in pins:
                    pending.append(target)
    return dict(sorted(pins.items()))

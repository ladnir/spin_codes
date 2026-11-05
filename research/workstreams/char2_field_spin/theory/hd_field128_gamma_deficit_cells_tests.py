"""Tiny exact geometry/count oracles, never production numerical gates."""

from dataclasses import FrozenInstanceError
from fractions import Fraction as F
from itertools import combinations, product
import json
from types import SimpleNamespace
import unittest

import hd_field128_gamma_deficit_cells as d
from hd_field128_gamma_cell_geometry import Cell
import hd_field128_gamma_cell_counts as old_counts


def budget():
    return SimpleNamespace(check=lambda: None)


def lattice(cell):
    return {point for point in product(*(range(lo, hi+1) for lo, hi in cell.ranges)) if cell.contains(point)}


def supports(pairs):
    rows = []
    for mask in range(1, (1 << (2*pairs))-1):
        active = [(mask >> i) & 1 for i in range(2*pairs)]
        r = sum(active)
        a = sum(active[2*i]*active[2*i+1] for i in range(pairs))
        b = sum(active[(2*i+1) % (2*pairs)]*active[(2*i+2) % (2*pairs)] for i in range(pairs))
        rows.append((r, a, b))
    return rows


def gaussian_vertices(cell):
    """Independent Fraction elimination, not the production Cramer helper."""
    found = set()
    for rows in combinations(cell.constraints, 3):
        matrix = [[F(x) for x in normal]+[F(rhs)] for normal, rhs in rows]
        for col in range(3):
            pivot = next((row for row in range(col, 3) if matrix[row][col]), None)
            if pivot is None:
                break
            matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
            scale = matrix[col][col]
            matrix[col] = [x/scale for x in matrix[col]]
            for row in range(3):
                if row != col:
                    scale = matrix[row][col]
                    matrix[row] = [x-scale*y for x, y in zip(matrix[row], matrix[col])]
        else:
            point = tuple(row[-1] for row in matrix)
            if cell.contains_real(point):
                found.add(point)
    return tuple(sorted(found))


class DeficitCellsTests(unittest.TestCase):
    def test_actual_support_parity_partition_with_separate_legacy_faces(self):
        for pairs in range(1, 6):
            roots = [d.DeficitCell(eps, ((0, pairs), (0, pairs), (0, pairs)), pairs) for eps in (0, 1)]
            faces = [Cell(kind, ((0, pairs), (0, pairs)), pairs) for kind in ('even_face', 'odd_face')]
            for r, a, b in supports(pairs):
                u, eps = divmod(r, 2)
                point = u, u-a, u-b
                native = roots[eps].contains(point)
                membership = [native, eps == 0 and a == u and faces[0].contains((u, r-a-b)),
                              eps == 0 and b == u and faces[1].contains((u, r-a-b))]
                self.assertEqual(sum(membership), 1)
                if native:
                    self.assertEqual(roots[eps].profile(point), (r, a, b))
                    self.assertEqual(r-a-b, point[1]+point[2]+eps)
            for eps, cell in enumerate(roots):
                for point in product(range(pairs+1), repeat=3):
                    u, x, y = point
                    expected = (1-eps <= u <= pairs-1 and 1-eps <= x <= min(u, pairs-u-eps)
                                and 1-eps <= y <= min(u, pairs-u-eps))
                    self.assertEqual(cell.contains(point), expected)

    def test_all_small_boxes_integer_feasibility_and_splits(self):
        pairs = 2
        intervals = tuple((lo, hi) for lo in range(-1, pairs+1) for hi in range(lo, pairs+1))
        for eps in (0, 1):
            for index, ranges in enumerate(product(intervals, repeat=3)):
                cell = d.DeficitCell(eps, ranges, pairs)
                points = lattice(cell)
                self.assertEqual(cell.integer_empty(), not points)
                self.assertIn(cell.integer_point(), points or {None})
                if not index % 47:
                    self.assertEqual(cell.real_empty(), not points)
                    self.assertEqual(cell.vertices(), gaussian_vertices(cell))
                for axis, (lo, hi) in enumerate(ranges):
                    if lo < hi:
                        left, right = cell.split(axis, (lo+hi)//2)
                        self.assertFalse(lattice(left) & lattice(right))
                        self.assertEqual(lattice(left) | lattice(right), points)

    def test_fractional_vertices_lower_dimensions_and_histograms(self):
        cells = [d.DeficitCell(0, ((1, 2), (1, 2), (1, 2)), 3),
                 d.DeficitCell(1, ((0, 3), (0, 3), (0, 3)), 4),
                 d.DeficitCell(1, ((1, 1), (0, 1), (0, 1)), 4),
                 d.DeficitCell(1, ((1, 1), (0, 1), (1, 1)), 4),
                 d.DeficitCell(1, ((0, 0), (0, 0), (0, 0)), 4),
                 d.DeficitCell(0, ((1, 1), (2, 2), (1, 1)), 4)]
        for cell in cells:
            self.assertEqual(cell.vertices(), gaussian_vertices(cell))
            for parity in (0, 1):
                bounds = cell.histogram_bounds(parity)
                if bounds is None:
                    self.assertTrue(cell.real_empty())
                    continue
                for point in lattice(cell):
                    r, a, b = cell.profile(point)
                    full = (a, b)[parity]
                    self.assertTrue(all(lo <= value <= hi for (lo, hi), value in zip(bounds, (r-2*full, full))))
        self.assertIn((F(3, 2), F(3, 2), F(3, 2)), cells[0].vertices())
        self.assertIn((F(3, 2), F(3, 2), F(3, 2)), cells[1].vertices())
        # Convex real samples, including noninteger weights, remain inside the original real domain.
        for first in cells[0].vertices():
            for second in cells[0].vertices():
                point = tuple((x+2*y)/3 for x, y in zip(first, second))
                r, a, b = cells[0].profile(point)
                self.assertTrue(cells[0].contains_real(point))
                self.assertTrue(1 <= r <= 5 and a >= r-3 and b >= r-3 and 2*a <= r-2 and 2*b <= r-2)

    def test_even_deficits_remove_old_phantom_half_profiles(self):
        old = Cell('interior', ((4, 4), (0, 2), (0, 2)), 4)
        fresh = d.DeficitCell(0, ((2, 2), (1, 2), (1, 2)), 4)
        self.assertEqual({fresh.profile(p) for p in lattice(fresh)}, lattice(old))
        self.assertEqual(old.histogram_bounds(0, integer=False)[0][0], 1)
        self.assertEqual(fresh.histogram_bounds(0, integer=False)[0][0], 2)
        self.assertNotIn((F(4), F(3, 2), F(3, 2)), tuple(fresh.profile(v) for v in fresh.vertices()))

    def test_projection_sums_equal_actual_support_partition_oracles(self):
        cache = d.PrefixCache(4)
        for pairs in range(2, 6):
            actual = supports(pairs)
            for eps in (0, 1):
                cells = [d.DeficitCell(eps, ((0, pairs), (0, pairs), (0, pairs)), pairs),
                         d.DeficitCell(eps, ((1, pairs-1), (0, 1), (0, 2)), pairs),
                         d.DeficitCell(eps, ((0, pairs), (1, 2), (0, 1)), pairs)]
                for cell in cells:
                    profiles = {cell.profile(point) for point in lattice(cell)}
                    even = {(r, a) for r, a, b in profiles}
                    odd = {(r, b) for r, a, b in profiles}
                    runs = {(r, r-a-b) for r, a, b in profiles}
                    oracle = [sum((r, a) in even for r, a, b in actual),
                              sum((r, b) in odd for r, a, b in actual),
                              sum((r, r-a-b) in runs for r, a, b in actual)]
                    report = d.count_diagnostics(cell, budget(), cache=cache)
                    self.assertEqual(list(report['projection_sums'].values()), oracle)
                    self.assertEqual(report['support_count_upper'], min(oracle))
                    self.assertGreaterEqual(report['support_count_upper'], sum(row in profiles for row in actual))
                    self.assertEqual(report['nonempty_integer_slices'], len({r for r, a, b in profiles}))
                    self.assertFalse(report['exact_support_count'])
                    self.assertFalse(report['singleton_free_faces_included'])
                    self.assertEqual(report['profile_pairs_enumerated'], 0)
                    self.assertLessEqual(len(cache), 4)

    def test_single_profile_replays_frozen_interior_count(self):
        for pairs in range(2, 6):
            for eps in (0, 1):
                root = d.DeficitCell(eps, ((0, pairs), (0, pairs), (0, pairs)), pairs)
                for point in lattice(root):
                    cell = d.DeficitCell(eps, tuple((x, x) for x in point), pairs)
                    profile = cell.profile(point)
                    old = Cell('interior', tuple((x, x) for x in profile), pairs)
                    self.assertEqual(d.support_count(cell, budget()), old_counts.support_count(old, budget()))

    def test_descriptors_guards_budget_and_type_separation(self):
        cell = d.DeficitCell(1, [[0, 3], [0, 2], [0, 2]], 4)
        self.assertEqual(d.DeficitCell.from_descriptor(json.loads(json.dumps(cell.descriptor()))), cell)
        self.assertNotIsInstance(cell, Cell)
        with self.assertRaises(FrozenInstanceError):
            cell.eps = 0
        for args in ((True, cell.ranges, 4), (2, cell.ranges, 4), (1, ((0., 3), (0, 2), (0, 2)), 4)):
            with self.assertRaises(ValueError):
                d.DeficitCell(*args)
        with self.assertRaises(ValueError):
            cell.contains((0., 0, 0))
        with self.assertRaises(ValueError):
            d.DeficitCell.from_descriptor(Cell('interior', ((1, 1), (0, 0), (0, 0)), 4).descriptor())
        with self.assertRaises(ValueError):
            d.support_count(Cell('interior', ((1, 1), (0, 0), (0, 0)), 4), budget())
        expired = SimpleNamespace(check=lambda: (_ for _ in ()).throw(TimeoutError('tiny budget')))
        with self.assertRaises(TimeoutError):
            d.support_count(cell, expired)
        pins = d.source_pins()
        self.assertEqual(set(pins), {'hd_field128_gamma_deficit_cells.py',
                                    'hd_field128_gamma_cell_counts.py', 'hd_field128_gamma_cell_geometry.py'})


if __name__ == '__main__':
    unittest.main()

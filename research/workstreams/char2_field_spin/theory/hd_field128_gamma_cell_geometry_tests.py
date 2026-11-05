"""Tiny exhaustive exact-geometry tests; no coefficient banks or gates."""

from dataclasses import FrozenInstanceError
from fractions import Fraction as F
from itertools import combinations, product
import json
import unittest

from hd_field128_gamma_cell_geometry import Cell, KINDS


def domain(kind, point, pairs):
    if kind == 'interior':
        r, a, b = point
        return (1 <= r <= 2*pairs-1 and 0 <= a <= pairs and 0 <= b <= pairs
                and 2*a <= r-1 and 2*b <= r-1 and a >= r-pairs
                and b >= r-pairs and r-a-b >= 1)
    u, k = point
    return 1 <= u <= pairs-1 and 1 <= k <= min(u, pairs-u)


def lattice(cell):
    return {point for point in product(*(range(lo, hi+1) for lo, hi in cell.ranges))
            if domain(cell.kind, point, cell.pairs)}


def independent_vertices(cell):
    """Slow exact Gaussian oracle, independent of production Cramer formulas."""
    result = set()
    n = cell.dimension
    for equations in combinations(cell.constraints, n):
        matrix = [[F(x) for x in coefficients]+[F(rhs)] for coefficients, rhs in equations]
        for col in range(n):
            pivot = next((row for row in range(col, n) if matrix[row][col]), None)
            if pivot is None:
                break
            matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
            scale = matrix[col][col]
            matrix[col] = [x/scale for x in matrix[col]]
            for row in range(n):
                if row != col:
                    scale = matrix[row][col]
                    matrix[row] = [x-scale*y for x, y in zip(matrix[row], matrix[col])]
        else:
            point = tuple(row[-1] for row in matrix)
            if all(sum(a*x for a, x in zip(coefficients, point)) <= rhs
                   for coefficients, rhs in cell.constraints):
                result.add(point)
    return tuple(sorted(result))


class CellGeometryTests(unittest.TestCase):
    def test_native_domains_and_integer_witnesses_exhaustively(self):
        for pairs in range(1, 5):
            for kind in KINDS:
                ranges = ((0, 2*pairs), (0, pairs), (0, pairs)) if kind == 'interior' else ((0, pairs), (0, pairs))
                cell = Cell(kind, ranges, pairs)
                expected = lattice(cell)
                for point in product(*(range(lo, hi+1) for lo, hi in ranges)):
                    self.assertEqual(cell.contains(point), point in expected)
                self.assertEqual(cell.integer_empty(), not expected)
                self.assertEqual(cell.real_empty(), not expected)
                self.assertIn(cell.integer_point(), expected or {None})
                # Every possible singleton box independently checks the O(1) witness path.
                for point in product(*(range(lo, hi+1) for lo, hi in ranges)):
                    single = Cell(kind, tuple((x, x) for x in point), pairs)
                    self.assertEqual(single.integer_point(), point if point in expected else None)

    def test_vertices_against_independent_exact_oracle(self):
        cells = [
            Cell('interior', ((1, 7), (0, 4), (0, 4)), 4),
            Cell('interior', ((2, 6), (0, 3), (1, 3)), 4),
            Cell('interior', ((2, 2), (0, 1), (0, 1)), 4),
            Cell('interior', ((3, 3), (0, 1), (1, 1)), 4),
            Cell('interior', ((3, 3), (1, 1), (1, 1)), 4),
            Cell('interior', ((2, 2), (1, 1), (0, 1)), 4),
            Cell('even_face', ((0, 5), (0, 5)), 5),
            Cell('odd_face', ((2, 2), (0, 4)), 5),
            Cell('even_face', ((1, 1), (1, 1)), 5),
            Cell('odd_face', ((1, 2), (3, 4)), 5),
        ]
        for cell in cells:
            self.assertEqual(cell.vertices(), independent_vertices(cell))
            self.assertTrue(all(cell.contains_real(v) for v in cell.vertices()))
        half_square = cells[2].vertices()
        self.assertEqual(half_square, tuple((F(2), a, b) for a, b in product((F(0), F(1, 2)), repeat=2)))
        self.assertFalse(cells[2].contains((2, F(1, 2), 0)))
        self.assertEqual(cells[6].vertices(), ((F(1), F(1)), (F(5, 2), F(5, 2)), (F(4), F(1))))

    def test_every_small_box_integer_feasibility(self):
        pairs = 2
        intervals = lambda high: tuple((lo, hi) for lo in range(-1, high+1) for hi in range(lo, high+1))
        for kind in KINDS:
            choices = ((intervals(2*pairs), intervals(pairs), intervals(pairs))
                       if kind == 'interior' else (intervals(pairs), intervals(pairs)))
            for index, ranges in enumerate(product(*choices)):
                cell = Cell(kind, ranges, pairs)
                points = lattice(cell)
                self.assertEqual(cell.integer_empty(), not points)
                self.assertIn(cell.integer_point(), points or {None})
                if not index % 41:
                    self.assertEqual(cell.real_empty(), not points)
                    self.assertEqual(cell.vertices(), independent_vertices(cell))

    def test_splits_partition_integer_points_not_real_points(self):
        for pairs in range(2, 5):
            for kind in KINDS:
                ranges = ((1, 2*pairs-1), (0, pairs), (0, pairs)) if kind == 'interior' else ((1, pairs-1), (1, pairs))
                cell = Cell(kind, ranges, pairs)
                expected = lattice(cell)
                for axis, (low, high) in enumerate(ranges):
                    for cut in range(low, high):
                        left, right = cell.split(axis, cut)
                        a, b = lattice(left), lattice(right)
                        self.assertFalse(a & b)
                        self.assertEqual(a | b, expected)
                        for child, points in ((left, a), (right, b)):
                            self.assertEqual(child.integer_empty(), not points)
                            self.assertIn(child.integer_point(), points or {None})
        cell = Cell('interior', ((2, 3), (0, 0), (0, 0)), 4)
        left, right = cell.split(0, 2)
        midpoint = (F(5, 2), 0, 0)
        self.assertTrue(cell.contains_real(midpoint))
        self.assertFalse(left.contains_real(midpoint) or right.contains_real(midpoint))

    def test_all_tiny_supports_partition_and_face_orientation(self):
        for pairs in range(2, 6):
            interior = Cell('interior', ((1, 2*pairs-1), (0, pairs), (0, pairs)), pairs)
            faces = [Cell(kind, ((1, pairs-1), (1, pairs)), pairs) for kind in KINDS[1:]]
            for mask in range(1, (1 << (2*pairs))-1):
                active = [(mask >> i) & 1 for i in range(2*pairs)]
                r = sum(active)
                a = sum(active[2*i]*active[2*i+1] for i in range(pairs))
                b = sum(active[(2*i+1) % (2*pairs)]*active[(2*i+2) % (2*pairs)] for i in range(pairs))
                k = r-a-b
                membership = [interior.contains((r, a, b)),
                              r % 2 == 0 and a == r//2 and faces[0].contains((r//2, k)),
                              r % 2 == 0 and b == r//2 and faces[1].contains((r//2, k))]
                self.assertEqual(sum(membership), 1)
                for face, member in zip(faces, membership[1:]):
                    if member:
                        self.assertEqual(face.profile((r//2, k)), (r, a, b))

    def test_histogram_envelopes_are_exact_and_contain_all_lattice_images(self):
        for kind in KINDS:
            ranges = ((2, 5), (0, 2), (0, 3)) if kind == 'interior' else ((1, 4), (1, 3))
            cell = Cell(kind, ranges, 5)
            points = lattice(cell)
            for parity in (0, 1):
                exact = cell.histogram_bounds(parity, integer=False)
                envelope = cell.histogram_bounds(parity)
                images = [(cell.profile(v)[0]-2*cell.profile(v)[parity+1], cell.profile(v)[parity+1])
                          for v in cell.vertices()]
                self.assertEqual(exact, tuple((min(row[i] for row in images), max(row[i] for row in images)) for i in (0, 1)))
                for point in points:
                    r, a, b = cell.profile(point)
                    full = (a, b)[parity]
                    self.assertTrue(all(lo <= value <= hi for (lo, hi), value in zip(envelope, (r-2*full, full))))
        half = Cell('interior', ((2, 2), (0, 1), (0, 1)), 4)
        self.assertEqual(half.histogram_bounds(0, integer=False), ((F(1), F(2)), (F(0), F(1, 2))))
        self.assertEqual(half.histogram_bounds(0), ((1, 2), (0, 0)))
        self.assertIsNone(Cell('even_face', ((1, 1), (2, 2)), 4).histogram_bounds(0))

    def test_descriptor_freezing_and_input_guards(self):
        cell = Cell('interior', [[1, 4], [0, 2], [0, 2]], 4)
        self.assertEqual(Cell(**json.loads(json.dumps(cell.descriptor()))), cell)
        self.assertEqual(cell.ranges, ((1, 4), (0, 2), (0, 2)))
        with self.assertRaises(FrozenInstanceError):
            cell.pairs = 5
        for kind, ranges, pairs in [('unknown', (), 4), ('interior', ((1, 2),), 4),
                                    ('even_face', ((1, 0), (1, 2)), 4),
                                    ('even_face', ((1., 2), (1, 2)), 4),
                                    ('even_face', ((1, 2), (1, 2)), True)]:
            with self.assertRaises(ValueError):
                Cell(kind, ranges, pairs)
        for point in ((1., 0, 0), (True, 0, 0), ('1', 0, 0), (1, 0)):
            with self.assertRaises(ValueError):
                cell.contains(point)
        for args in ((0, 4), (0, 0), (3, 1), (True, 1), (0, F(2))):
            with self.assertRaises(ValueError):
                cell.split(*args)
        with self.assertRaises(ValueError):
            cell.histogram_bounds(True)


if __name__ == '__main__':
    unittest.main()

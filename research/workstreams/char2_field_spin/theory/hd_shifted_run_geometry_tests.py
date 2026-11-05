from fractions import Fraction as F
from itertools import product
from math import factorial
import unittest

from hd_shifted_run_geometry import (
    run_box_profile_count, run_box_receipt, run_box_support_count_upper,
    run_box_vertices, validate_run_box,
)


def toy_supports(pairs):
    values = []
    for bits in product((0, 1), repeat=2*pairs):
        r = sum(bits)
        a = sum(bits[2*i]*bits[2*i+1] for i in range(pairs))
        b = sum(bits[2*i+1]*bits[(2*i+2) % (2*pairs)] for i in range(pairs))
        values.append((r, r-a-b, a-b))
    return values


class RunGeometryTests(unittest.TestCase):
    def test_requested_box_retains_eight_exact_vertices_and_parity(self):
        box = ((992, 1056), (220, 228), (-4, 4))
        self.assertEqual(len(run_box_vertices(box)), 8)
        self.assertEqual(run_box_profile_count(box), 2633)
        half_box = ((1023, 1023), (224, 224), (0, 0))
        self.assertEqual(run_box_vertices(half_box), ((F(1023), F(799,2), F(799,2)),))
        self.assertEqual(run_box_profile_count(half_box), 0)

    def test_exhaustive_small_all_valid_boxes_count_exact_profiles(self):
        for pairs in range(1, 5):
            physical = set(toy_supports(pairs))
            intervals = lambda lo, hi: [(a, b) for a in range(lo, hi+1) for b in range(a, hi+1)]
            domains = (intervals(1, 2*pairs-1), intervals(1, pairs), intervals(-pairs, pairs))
            checked = 0
            for box in product(*domains):
                try:
                    validate_run_box(box, pairs)
                except ValueError:
                    continue
                expected = sum(all(lo <= value <= hi for value, (lo, hi) in zip(p, box))
                               for p in physical)
                self.assertEqual(run_box_profile_count(box, pairs), expected, (pairs, box))
                for r, a, b in run_box_vertices(box, pairs):
                    self.assertGreaterEqual(min(pairs-r+a, r-2*a, a, pairs-r+b, r-2*b, b), 0)
                checked += 1
            self.assertGreater(checked, 0)

    def test_rejects_boundary_crossings_and_constant_run_domain(self):
        for box in (((1, 3), (1, 2), (-2, 2)), ((1, 4), (0, 1), (0, 0)),
                    ((-1, 1), (1, 1), (0, 0))):
            with self.assertRaises(ValueError):
                run_box_vertices(box, 2)

    def test_whole_box_transfer_counts_supports_without_profile_multiplier(self):
        pairs = 4
        box = ((3, 5), (1, 2), (-1, 1))
        values = toy_supports(pairs)
        count = sum(all(lo <= value <= hi for value, (lo, hi) in zip(p, box)) for p in values)
        for x, y in ((F(1), F(1)), (F(1,2), F(2)), (F(2), F(1,3))):
            self.assertGreaterEqual(run_box_support_count_upper(box, x, y, pairs), count)
        self.assertEqual(run_box_support_count_upper(box, F(1), F(1), pairs), 2**(2*pairs))
        empty = ((3, 3), (2, 2), (0, 0))
        self.assertEqual(run_box_support_count_upper(empty, F(1), F(1), pairs), 0)

    def test_coupled_run_tilt_bounds_weighted_support_sum(self):
        from hd_shifted_run_count import transfer_value
        pairs = 4
        box = ((3, 5), (1, 2), (-1, 1))
        vertices = list(product(*box))
        mass = lambda r, c, d: F(factorial(r)*factorial(c)*factorial(d+pairs))
        values = [p for p in toy_supports(pairs)
                  if all(lo <= value <= hi for value, (lo, hi) in zip(p, box))]
        actual = sum((mass(*p) for p in values), F(0))
        for x, y in ((F(1,2), F(2)), (F(2), F(1,3)), (F(1), F(1))):
            generating = transfer_value(pairs, x, y)
            coupled = generating*max(mass(r,c,d)/(x**r*y**c) for r,c,d in vertices)
            separate = generating/min(x**r*y**c for r,c,d in vertices)*max(mass(*p) for p in vertices)
            self.assertGreaterEqual(coupled, actual)
            self.assertLessEqual(coupled, separate)

    def test_receipt_uses_alias_ranges_and_no_extra_profile_factor(self):
        from flint import arb
        targets = []

        class Atlas:
            def numerator(self, ranges):
                targets.append(ranges)
                return (arb(1), arb(0), arb(0), arb(1))

        class Witness:
            precision = 128
            grid = 128
            uniform_aliases = False
            regions = (Atlas(), Atlas())

            def scalar(self, point, count_matching=None):
                if count_matching is not None:
                    raise AssertionError('Run support sum must use the count-free scalar')
                return arb(0)

        result = run_box_receipt(Witness(), ((992,1056),(220,228),(-4,4)), F(1), F(1))
        self.assertEqual(targets, [((216,232),(380,420))]*2)
        self.assertEqual(result['integer_profile_triples'], 2633)
        ceiling = F(result['log2_upper_rational'])
        self.assertGreaterEqual(ceiling, 8192)
        self.assertLessEqual(ceiling, F(8192)+F(1,10**9))


if __name__ == '__main__':
    unittest.main()

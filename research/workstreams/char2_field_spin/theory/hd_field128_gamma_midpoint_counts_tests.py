"""Independent tiny support enumeration and all boundary/count identities."""

from fractions import Fraction as F
from math import log2
import unittest

import numpy as np

import hd_field128_gamma_midpoint_counts as mode
from hd_shifted_pairs_count import exact_profiles
from hd_field128_gamma_near_face_certificate_tests import fixture, multiply, power


class MidpointTests(unittest.TestCase):
    def test_exhaustive_supports_cap_and_factor(self):
        for p in range(1, 8):
            oracle = exact_profiles(p)
            for r in range(2*p+1):
                c = mode.MidpointCounts(p, r)
                for a in c.full:
                    for b in c.full:
                        actual = oracle.get((r, a, b), 0)
                        cap = c.midpoint_upper(a, b)
                        aa = max(a, b)-c.shift
                        self.assertLessEqual(actual, cap)
                        self.assertLessEqual(cap, (c.r-2*aa+1)*actual)
                        self.assertLessEqual(actual, c.count(a, b))
                        self.assertLessEqual(c.count(a, b), c.base.count(a, b))
                        self.assertEqual(c.count(a, b), c.count(b, a))

    def test_midpoint_and_monotonicity_against_fixed_flows(self):
        for p in range(1, 25):
            for r in range(p+1):
                c = mode.MidpointCounts(p, r)
                for a in c.full:
                    for b in range(a+1):
                        if r == p and a == b == 0:
                            self.assertEqual(c.midpoint_upper(a, b), 2)
                            continue
                        s = r-2*a
                        flows = [mode.exact.fixed_flow_count(p, r, a, b, u) for u in range(s//2+1)]
                        self.assertTrue(all(x <= y for x, y in zip(flows, flows[1:])))
                        self.assertEqual(c.midpoint_upper(a, b), (s+1)*flows[-1])

    def test_vector_logs_and_complement(self):
        for p in range(2, 17):
            for r in range(2*p+1):
                c = mode.MidpointCounts(p, r)
                inverse = mode.MidpointCounts(p, 2*p-r)
                full = np.array(list(c.full))
                values = c.log_counts(full[:, None], full[None, :])
                for i, a in enumerate(c.full):
                    for j, b in enumerate(c.full):
                        count = c.count(a, b)
                        expected = log2(count) if count else -np.inf
                        self.assertAlmostEqual(float(values[i, j]), expected, places=10)
                        self.assertEqual(count, inverse.count(p-r+a, p-r+b))
                self.assertEqual(float(c.log_counts(int(full[0]), int(full[0]))), float(values[0, 0]))

    def test_float_transfer_and_guards(self):
        for p, r in ((4, 3), (5, 5), (6, 8)):
            rows = fixture(p, r)
            c = mode.MidpointCounts(p, r)
            expected = F(0)
            for a, left in rows.items():
                for b, right in rows.items():
                    matrix = power(multiply(left, right), 4)
                    expected += c.count(a, b)*(matrix[0]+matrix[1])
            regions = np.array([[log2(v) for v in row] for row in rows.values()]).reshape(-1, 2, 2)
            result = mode.support_sum(p, r, regions, repetitions=8, counts=c)
            self.assertAlmostEqual(result['log2_support_sum_screen'], log2(expected), places=10)
            self.assertEqual(result['profile_pairs_evaluated'], len(rows)**2)
            self.assertFalse(result['directed'])
        with self.assertRaises(ValueError):
            mode.MidpointCounts(0, 0)
        with self.assertRaises(ValueError):
            c.count(True, 0)
        self.assertEqual(c.count(-1, 0), 0)
        self.assertIn('../HD_FIELD128_MIDPOINT_COUNTS.md', mode.source_pins())


if __name__ == '__main__':
    unittest.main()

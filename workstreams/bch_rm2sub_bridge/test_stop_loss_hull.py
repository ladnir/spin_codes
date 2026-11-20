"""An exact test that convexifying upper call bounds preserves domination."""
from fractions import Fraction as F
import unittest
from certify_johnson_convex_law import lower_convex_hull


class StopLossHullTests(unittest.TestCase):
    def test_nonconvex_bounds(self):
        laws = [[F(1, 2), 0, 0, 0, F(1, 2)], [0, F(1, 2), 0, F(1, 2), 0]]
        true_calls = [[sum(max(a-t, 0)*p for a, p in enumerate(law)) for t in range(5)] for law in laws]
        bounds = [max(call[t] for call in true_calls) for t in range(5)]
        bounds[2] += F(3, 4)
        _, hull = lower_convex_hull(bounds)
        self.assertTrue(all(call[t] <= hull[t] <= bounds[t] for call in true_calls for t in range(5)))
        self.assertEqual(hull[0], 2)
        self.assertEqual(hull[-1], 0)
        slopes = [hull[t+1]-hull[t] for t in range(4)]
        masses = [1+slopes[0]]+[slopes[t]-slopes[t-1] for t in range(1, 4)]+[-slopes[-1]]
        self.assertTrue(all(p >= 0 for p in masses))
        self.assertEqual(sum(masses), 1)
        self.assertEqual(sum(a*p for a, p in enumerate(masses)), 2)


if __name__ == '__main__':
    unittest.main()

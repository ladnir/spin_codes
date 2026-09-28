"""Exact small-domain tests and outward checks for shell/CDF folding."""
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from itertools import combinations_with_replacement, product
from math import comb, log
from types import SimpleNamespace
import unittest

import numpy as np
from flint import arb, arb_mat, ctx

from shell_cover import IntervalFolds, prefix_rank_bound, cover
from occupancy_cdf_cover import fold_function, fold_arb, geometry_test, retained_test


def rational(point):
    value = point.fmpq()
    return Q(int(value.p), int(value.q))


def cdf_bound(cdf, weights):
    increments = [cdf[0], *(b - a for a, b in zip(cdf, cdf[1:]))]
    return sum(c * max(weights[i:]) for i, c in enumerate(increments))


class ShellCoverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 192

    def test_prefix_bound_matches_tiny_exact_optimum(self):
        # Prefix constraints and individual caps have integral vertices.
        # These cases also check the nested-set identity independently by
        # enumerating every feasible integer vector, rather than using it.
        weights_list = [(Q(1), Q(2), Q(3)), (Q(3), Q(2), Q(1)),
                        (Q(3, 7), Q(1, 2), Q(2, 9)), (Q(4), Q(1), Q(4))]
        for cdf in combinations_with_replacement(range(5), 3):
            for shells in product(range(3), repeat=3):
                feasible = [x for x in product(*(range(s + 1) for s in shells))
                            if all(sum(x[:i + 1]) <= cdf[i] for i in range(3))]
                for weights in weights_list:
                    actual = max(sum(xi * wi for xi, wi in zip(x, weights)) for x in feasible)
                    bound = prefix_rank_bound(cdf, shells, weights)
                    self.assertEqual(actual, bound)
                    self.assertLessEqual(bound, cdf_bound(cdf, weights))
                    self.assertLessEqual(bound, sum(s * w for s, w in zip(shells, weights)))

    def test_joint_caps_can_beat_both_individual_folds(self):
        cdf, shells, weights = [1, 3, 4], [10, 1, 10], [3, 10, 2]
        bound = prefix_rank_bound(cdf, shells, weights)
        self.assertEqual(bound, 17)
        self.assertLess(bound, cdf_bound(cdf, weights))
        self.assertLess(bound, sum(s * w for s, w in zip(shells, weights)))

    def test_interval_counts_and_outward_arithmetic(self):
        n = 4
        for actual in product(range(3), repeat=3):
            actual = [0, *actual, 0]
            cdf = [Q(sum(actual[:i + 1])) + Q(i + 1, 3) for i in range(n + 1)]
            shells = [Q(x) + Q((i + 1) % 3, 5) for i, x in enumerate(actual)]
            for prefix in (False, True):
                folds = IntervalFolds(cdf, shells, n=n, prefix_rank=prefix)
                for lo in range(n + 1):
                    for hi in range(lo, n + 1):
                        for p in (Q(1, 10), Q(1, 2), Q(9, 10)):
                            weights = [1 / (comb(n, u) * p**u * (1 - p)**(n - u)) for u in range(lo, hi + 1)]
                            exact = sum(Q(actual[u]) * weights[u - lo] for u in range(lo, hi + 1))
                            bound = folds.outward(lo, hi, arb(p.numerator) / p.denominator)
                            self.assertGreaterEqual(rational(bound), exact)
                            proposed = folds.log(lo, hi, float(p))
                            if bound > 0:
                                self.assertAlmostEqual(proposed, log(rational(bound)), places=11)
                            else:
                                self.assertEqual(proposed, -np.inf)

    def test_no_shells_preserves_existing_fold(self):
        counts = [u * u + 1 for u in range(257)]
        folds = IntervalFolds(counts)
        for lo, hi in ((38, 256), (66, 92), (128, 128), (200, 256)):
            old = fold_function(counts, lo, hi)
            for p in (.0004, .1, .5, .9, .999999):
                self.assertAlmostEqual(folds.log(lo, hi, p), old(p), places=10)
                self.assertEqual(folds.outward(lo, hi, arb(p)), fold_arb(counts, lo, hi, arb(p)))

    def test_rational_caps_and_large_counts(self):
        counts = [Q(2**2000 * (u + 1), 7) for u in range(5)]
        shells = [Q(2**1998, 11)] * 5
        fold = IntervalFolds(counts, shells, n=4, prefix_rank=True)
        proposal = fold.log(0, 4, .5)
        bound = fold.outward(0, 4, arb(1) / 2)
        self.assertTrue(np.isfinite(proposal))
        exact_bound = rational(bound)
        self.assertAlmostEqual(proposal, log(exact_bound.numerator) - log(exact_bound.denominator), places=10)

    def test_invalid_caps_rejected(self):
        with self.assertRaises(AssertionError):
            IntervalFolds([0, 2, 1], n=2)
        with self.assertRaises(AssertionError):
            IntervalFolds([0, 1, 2], [-1, 0, 1], n=2)
        with self.assertRaises(AssertionError):
            IntervalFolds([0, 1], n=2)

    def test_complete_cover_replay_and_old_geometry(self):
        args = SimpleNamespace(groups=1, target_bits=40, max_splits=0,
                               retain_parents=True, joint_witness=False, joint_top=2,
                               probe_supports=[], probe_vector=[], screen_only=False,
                               precision=192)
        value = arb(2)**-16
        operators = {('.001', '1'): ([arb_mat([[value]]), arb_mat([[value]])],
                                    [np.array([[2.**-16]]), np.array([[2.**-16]])])}
        counts = {'1': [u + 1 for u in range(257)]}
        shells = {'1': [Q(1)] * 257}
        with redirect_stdout(StringIO()):
            geometry_test()
            retained_test()
            result = cover(args, operators, counts, [1], shells, prefix_rank=True)
        self.assertTrue(0 < result < arb(2)**-40)


if __name__ == '__main__':
    unittest.main()

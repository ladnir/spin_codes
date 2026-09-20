"""Combinatorial checks of the independent zero-state lower-bound argument."""
import itertools
import math
import unittest

import numpy as np
from scipy.stats import binom

import bch_zero_state_lower_v1 as lower
import refine_bch_dense_v1 as dense


class LowerTests(unittest.TestCase):
    def test_even_occupation_parity_bound_by_enumeration(self):
        for block in (3, 4, 5):
            for weight in range(1, block):
                words = [sum(1 << j for j in support) for support in itertools.combinations(range(block), weight)]
                for q in (2, 4):
                    counts = {0: 1}
                    for _ in range(q):
                        next_counts = {}
                        for value, mass in counts.items():
                            for word in words:
                                next_counts[value ^ word] = next_counts.get(value ^ word, 0)+mass
                        counts = next_counts
                    self.assertGreaterEqual(counts.get(0, 0)*2**(block-1), len(words)**q)

    def test_concentration_charge_against_exact_binomial_tails(self):
        for q, p, block in ((1000, .2, 64), (5000, .17, 128), (10000, .5, 128)):
            lo, hi, good = lower.concentrated_range(q, p, block)
            outside = np.logaddexp(binom.logcdf(lo-1, q, p), binom.logsf(hi, q, p))+math.log(block)
            allowed = -(block+4)*math.log(2)
            self.assertLessEqual(outside, allowed+1e-9)
            self.assertAlmostEqual(good, (1-block)*math.log(2)+math.log1p(-1/32))

    def test_kernel_coefficients_against_exact_integer_convolution(self):
        kernel = [1, 0, 2, 0, 1]
        actual = lower.kernel_region_logs(kernel, 3, 12)
        coefficients = np.polynomial.polynomial.polypow([1, 2, 1], 3)
        expected = np.array([math.log(int(n)/math.comb(12, 2*j)) for j, n in enumerate(coefficients)])
        np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-14)

    def test_convex_lower_hull_for_all_small_profiles(self):
        xs = np.arange(0, 9, 2); ys = np.array([0., -2., -1., -4., -3.])
        for indices in itertools.product(range(len(xs)), repeat=4):
            mean = sum(xs[j] for j in indices)/4
            value, _ = lower.convex_lower_value(xs, ys, mean)
            self.assertLessEqual(4*value, sum(ys[j] for j in indices)+1e-12)

    def test_exact_bounded_type_volume(self):
        for lo in ((0, 0, 0), (1, 0, 2)):
            for hi in ((3, 4, 4), (2, 2, 3)):
                for total in range(1, 9):
                    exact = sum(sum(x) == total for x in itertools.product(*(range(a, b+1) for a, b in zip(lo, hi))))
                    self.assertEqual(exact, dense.lattice_count(list(lo), list(hi), total))


if __name__ == '__main__':
    unittest.main()

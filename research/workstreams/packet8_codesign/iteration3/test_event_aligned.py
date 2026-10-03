"""Tiny exhaustive noncommuting checks for event-pair fractional placement."""
from itertools import product
from math import comb
import unittest

import numpy as np

import event_aligned as event


def direct_blocks(counts, local, alpha):
    """Independent sequence parser: cut after every second nonempty step."""
    size = local.shape[1]
    answer, block = np.eye(size), np.eye(size)
    occupied = 0
    for j in counts:
        block = block@local[j]
        occupied += j > 0
        if occupied == 2:
            answer = answer@(block**alpha)
            block, occupied = np.eye(size), 0
    return answer@(block**alpha)


class EventAlignedTests(unittest.TestCase):
    def setUp(self):
        self.local = np.array([[[1., 0.], [0., .27]],
            [[.12, .21], [.04, .18]], [[.03, .09], [.11, .07]]])

    def test_all_short_sequences_and_terminal_tails(self):
        W, E = 2, 5
        for alpha in (.4, 1.):
            actual, metadata = event.event_regional(self.local, W*E,
                epochs=E, alpha=alpha)
            expected = np.zeros_like(actual)
            for counts in product(range(W+1), repeat=E):
                q = sum(counts)
                coefficient = 1
                for j in counts:
                    coefficient *= comb(W, j)
                expected[q] += coefficient*direct_blocks(counts, self.local, alpha)
            expected /= np.array([comb(W*E, q) for q in range(W*E+1)])[:, None, None]
            np.testing.assert_allclose(np.exp(actual), expected, rtol=2e-13, atol=1e-15)
            self.assertLess(metadata['alpha1_regression']['max_regional_log_error'], 1e-12)

    def test_scalar_unit_mass_and_q0(self):
        regional, _ = event.event_regional(np.ones((3, 1, 1)), 12, epochs=6, alpha=.4)
        np.testing.assert_allclose(regional, 0., atol=3e-14)
        regional, _ = event.event_regional(self.local, 0, epochs=6, alpha=.4)
        np.testing.assert_allclose(np.exp(regional[0]), np.linalg.matrix_power(self.local[0], 6)**.4,
            rtol=1e-13, atol=1e-15)

    def test_positive_path_fractional_majorant(self):
        for counts in product(range(3), repeat=4):
            moment = np.eye(2)
            for j in counts:
                moment = moment@self.local[j]
            bound = direct_blocks(counts, self.local, .4)
            self.assertLessEqual(moment[0].sum()**.4, bound[0].sum()+2e-14)

    def test_log_blocks_preserve_values_that_linear_products_lose(self):
        local = np.array([[[1., 0.], [0., 1e-100]],
            [[.2, .1], [.1, .2]]])
        _, tail = event.event_blocks(local, 8, [.25])
        self.assertTrue(np.isfinite(tail[0, 8, 0, 1, 1]))
        self.assertAlmostEqual(tail[0, 8, 0, 1, 1], np.log(1e-100)*2, places=12)
        self.assertGreater(np.exp(tail[0, 8, 0, 1, 1]), 0.)
        self.assertEqual(np.linalg.matrix_power(local[0], 8)[1, 1], 0.)

    def test_guarded_product_agrees_with_full_log_and_falls_back(self):
        rng = np.random.default_rng(7401)
        for spread in (2., 500., 2000.):
            left = -rng.random((3, 2, 2))*spread
            right = -rng.random((3, 2, 2))*spread
            left[0, 0, 0] = -np.inf
            statistics = dict(full_log_products=0, guarded_linear_products=0,
                smallest_guarded_scalar_log_lower=0.)
            actual = event.log_product(left, right, statistics=statistics)
            expected = event.log_product(left, right, force_log=True)
            np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=1e-12)
        self.assertGreater(statistics['full_log_products'], 0)
        fast, _ = event.event_regional(self.local, 6, epochs=4, alpha=.4)
        slow, _ = event.event_regional(self.local, 6, epochs=4, alpha=.4, force_log_products=True)
        np.testing.assert_allclose(fast, slow, rtol=2e-14, atol=2e-13)

    def test_invalid_parameters(self):
        for alpha in (0., -1., 1.01, float('nan')):
            with self.assertRaises(ValueError):
                event.event_regional(self.local, 2, epochs=2, alpha=alpha)
        for q in (-1, 5):
            with self.assertRaises(ValueError):
                event.event_regional(self.local, q, epochs=2)


if __name__ == '__main__':
    unittest.main()

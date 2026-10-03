"""Exact comparison-moment identity for the potential-birth basis."""
from itertools import product
import unittest

import numpy as np

import potential_birth_gate as birth


class PotentialBirthTests(unittest.TestCase):
    def setUp(self):
        weighted = np.array([[1., .1, .02], [0., .02, .01],
            [0., .1, .03], [0., .03, .04]])
        emission = np.array([[1., .5, .25], [.3, .2, .15],
            [.7, .3, .12], [.4, .1, .05]])
        local = birth.prior.operators(weighted, emission)
        self.old = birth.marked.potential_operators(local, 0., packet_bits=2)
        self.new, self.change, self.check = birth.rebase_potential(self.old)

    def test_stochastic_intertwining(self):
        np.testing.assert_allclose(self.change.sum(axis=1), 1., atol=2e-15)
        np.testing.assert_array_equal(self.change[0], [1., 0., 0., 0.])
        np.testing.assert_allclose(self.new@self.change, self.change@self.old,
            rtol=2e-14, atol=1e-15)
        self.assertLess(self.check['intertwining_max_relative_error'], 2e-14)

    def test_exactly_one_birth_coordinate_per_occupancy(self):
        np.testing.assert_array_equal(self.new[0, 0, 2:], 0.)
        for j in (1, 2):
            self.assertGreater(self.new[j, 0, j+1], 0.)
            row = self.new[j, 0, 2:].copy()
            row[j-1] = 0.
            np.testing.assert_array_equal(row, 0.)

    def test_all_short_chronological_moments_match(self):
        for counts in product(range(3), repeat=5):
            old = np.array([1., 0., 0., 0.])
            new = old.copy()
            for j in counts:
                old = old@self.old[j]
                new = new@self.new[j]
            np.testing.assert_allclose(new@self.change, old, rtol=3e-14, atol=1e-15)
            self.assertAlmostEqual(new.sum(), old.sum(), places=14)

    def test_alpha1_conditional_placement_matches(self):
        old, _ = birth.prior.placement(self.old, 8, epochs=4, windows=2)
        new, _ = birth.prior.placement(self.new, 8, epochs=4, windows=2)
        for q in range(9):
            old_value = birth.prior.logarithmic.log_power_matrix(old[q], 3)
            new_value = birth.prior.logarithmic.log_power_matrix(new[q], 3)
            self.assertAlmostEqual(old_value, new_value, places=12)


if __name__ == '__main__':
    unittest.main()

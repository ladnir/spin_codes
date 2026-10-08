"""Composition checks; these do not certify a global distance bound."""

from math import isinf
import unittest

from hd_randomized_outer_refined_profile import occupancy_counts, refined_profiles
from hd_randomized_outer_state_profile import local_profiles


class RefinedProfileTests(unittest.TestCase):
    def test_occupancy_normalization(self):
        for j, counts in enumerate(occupancy_counts()):
            self.assertAlmostEqual(sum(2**weight for u, weight in counts), 1, places=12)
            self.assertTrue(all(0 <= u <= min(j, 8) for u, weight in counts))

    def test_disabled_refinements_are_identity(self):
        self.assertEqual(refined_profiles(-1/64, joint=False, flat=False, returns=False),
                         local_profiles(-1/64))

    def test_pointwise_minima_and_structural_zeros(self):
        for log_z in (-1/64, -8):
            old, new = local_profiles(log_z), refined_profiles(log_z)
            for j in range(65):
                for a in range(3):
                    for b in range(3):
                        self.assertLessEqual(new[j][a][b], old[j][a][b])
            self.assertEqual(new[0][0][0], 0)
            for a, b in ((0, 1), (0, 2), (1, 0), (2, 0)):
                self.assertTrue(isinf(new[0][a][b]) and new[0][a][b] < 0)
            self.assertLess(new[32][1][1], old[32][1][1])


if __name__ == '__main__':
    unittest.main()

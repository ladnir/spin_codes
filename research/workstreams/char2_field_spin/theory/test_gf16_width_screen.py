"""Checks for the optional shorter-step screen; no certificate assertions."""

import unittest
import numpy as np

from gf16_sizes import local_model, screen_bound
from gf16_width_screen import local, screen, singleton_weights


class WidthScreenTests(unittest.TestCase):
    def test_sixteen_matches_preserved_model(self):
        for state in (2, 4, 6, 8):
            for z in (.4, .9999):
                np.testing.assert_allclose(local(z, 16, state), local_model(z, state), rtol=1e-14, atol=0)

    def test_same_geometry_matches_exact_and_envelope_sources(self):
        for active in (1, 2, 5, 9):
            actual = screen(16, active, .999, 16, 4)
            expected = screen_bound(16, active, .999, state=4, source_rho=.871)
            self.assertAlmostEqual(actual, expected, places=9)

    def test_eight_point_additive_profiles(self):
        for state in range(1, 5):
            weights = singleton_weights(8, state)  # also checks CA=0 exactly
            self.assertEqual(len(weights), 8)
            self.assertTrue(all(1 <= w <= 8 for w in weights))
            np.testing.assert_allclose(np.sum(local(1., 8, state), axis=2), 1., atol=1e-14)
        with self.assertRaises(ValueError):
            singleton_weights(8, 5)


if __name__ == '__main__':
    unittest.main()

import unittest

import numpy as np
from packet_regional_log import array_logs, log_placement
from packet_regional_supported import placement_logs


class SupportedPlacementTests(unittest.TestCase):
    def test_noncommuting_matrix_family_matches_log(self):
        operators = np.array([[[.9, .04], [0, .6]], [[.01, .3], [.08, .5]], [[.06, 0], [.1, .3]]])
        expected = log_placement(array_logs(operators), 6, epochs=4, windows=2)
        actual = placement_logs(operators, 6, epochs=4, windows=2)
        np.testing.assert_array_equal(np.isneginf(actual), np.isneginf(expected))
        np.testing.assert_allclose(actual, expected, atol=1e-11)
        rowwise = placement_logs(operators, 6, epochs=4, windows=2, rowwise=True)
        np.testing.assert_array_equal(np.isneginf(rowwise), np.isneginf(expected))
        np.testing.assert_allclose(rowwise, expected, atol=1e-11)

    def test_structural_zeros_stay_zero(self):
        local = np.array([[[1., 0], [0, .5]], [[.2, .1], [0, .4]]])
        actual = placement_logs(local, 2, epochs=3, windows=1)
        self.assertTrue(np.all(np.isneginf(actual[:, 1, 0])))
        self.assertTrue(np.isneginf(actual[0, 0, 1]))
        self.assertTrue(np.all(np.isfinite(actual[1:, 0, 1])))

    def test_lost_positive_entry_is_rejected(self):
        local = np.array([[[1., 0], [0, 1e-200]], [[1., 0], [0, 1e-200]]])
        with self.assertRaises(FloatingPointError):
            placement_logs(local, 2, epochs=2, windows=1)

    def test_negligible_summand_does_not_remove_support(self):
        # Some products underflow, but every completed entry has another
        # positive contribution. There is no lost edge to propagate.
        local = np.array([[[1e-200, 1], [1, 1e-200]], [[1e-200, 1], [1, 1e-200]]])
        got = placement_logs(local, 2, epochs=3, windows=1)
        expected = log_placement(array_logs(local), 2, epochs=3, windows=1)
        np.testing.assert_allclose(got, expected, atol=1e-10)

    def test_row_scaling_preserves_tiny_nonzero_initial_state(self):
        local = np.array([[[1., 0], [0, .2]], [[1., 0], [0, .2]]])
        with self.assertRaises(FloatingPointError):
            placement_logs(local, 1, epochs=1000, windows=1)
        got = placement_logs(local, 1, epochs=1000, windows=1, rowwise=True)
        np.testing.assert_allclose(got[:, 0, 0], 0., atol=1e-10)
        np.testing.assert_allclose(got[:, 1, 1], 1000*np.log(.2), atol=1e-10)
        self.assertTrue(np.all(np.isneginf(got[:, 0, 1])))


if __name__ == '__main__':
    unittest.main()

import math
import unittest

import numpy as np

from gf8_packet_sparse_screen import log_coefficients, bound
from gf8_packet_width_screen import local_matrix


class PacketSparseScreenTests(unittest.TestCase):
    def test_one_packet_is_local_polynomial(self):
        y, tau = .9, 5.
        actual = np.exp(log_coefficients(8, y, tau, 8, 8, 4))
        expected = np.asarray([np.asarray(local_matrix(w, y, tau, 8, 4))*math.comb(8, w)
                               for w in range(9)])
        np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=1e-15)

    def test_truncation_preserves_needed_coefficients(self):
        full = log_coefficients(16, .95, 7., 16, 8, 4)
        truncated = log_coefficients(3, .95, 7., 16, 8, 4)
        np.testing.assert_allclose(truncated, full[:4], rtol=1e-14, atol=1e-14)

    def test_width16_control_matches_frozen_joint(self):
        from joint_outer_screen import log_coefficients as old_coefficients
        from state8_joint_return import minimum_return_weight
        minima = tuple(minimum_return_weight(w, 16) for w in range(1, 17))
        actual = log_coefficients(7, .97, 11., 64, 16, 8, minima)
        expected = old_coefficients(7, .97, 11., 64)
        np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=1e-13)

    def test_stronger_return_cannot_weaken_bound(self):
        point = (6.4, math.log(16))
        generic = bound(8, point)
        stronger = bound(8, point, minima=(5, 5, 4, 2, 1, 1, 1, 1))
        self.assertLessEqual(stronger, generic)

    def test_invalid_geometry(self):
        with self.assertRaises(ValueError):
            log_coefficients(2, .9, 2, 15, 8, 4)


if __name__ == '__main__':
    unittest.main()

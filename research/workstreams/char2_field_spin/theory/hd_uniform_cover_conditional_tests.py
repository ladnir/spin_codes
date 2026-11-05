import unittest
import numpy as np

from hd_bit_syndrome_gate import weighted_cover_matrices
from hd_post_refresh_quiet_rows import log_matrix_power
from hd_uniform_cover_conditional import column_parts, norm
from fractions import Fraction


class ConditionalCoverTests(unittest.TestCase):
    def setUp(self):
        self.measure = dict(log2_probability=-2., log2_atom00=-1.,
            log2_atom01=-.5, log2_atom10=-.75,
            log2_scale01=-2., log2_scale10=-3., log2_scale11=-1.)

    def test_positive_excess_matches_direct_power(self):
        for a in (1, 4, 7, 11, 16):
            for r in (1, 3):
                constant, excess = column_parts(self.measure, a, r)
                bank = weighted_cover_matrices(self.measure, r, (Fraction(a, 16),))
                direct = log_matrix_power(bank['cover_matrices'][0], 16)
                np.testing.assert_allclose(np.logaddexp2(constant, excess), direct,
                                           rtol=1e-13, atol=1e-12)

    def test_marker_endpoint_is_unmarked(self):
        parts = column_parts(self.measure, 4, 3)
        self.assertTrue(np.isfinite(norm(parts, 4, 3, 0.)))

    def test_reject_bad_dimensions(self):
        for a, r in ((0, 1), (17, 1), (1, 0)):
            with self.assertRaises(ValueError):
                column_parts(self.measure, a, r)


if __name__ == '__main__':
    unittest.main()

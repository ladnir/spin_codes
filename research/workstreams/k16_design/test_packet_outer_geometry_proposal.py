import unittest
from fractions import Fraction as Q
from math import comb, log
import numpy as np
from flint import arb, arb_mat
from packet_outer_geometry_proposal import estimate, upper_arrays
from rs_uniform_envelope import UniformInputEnvelope


class OuterGeometryTests(unittest.TestCase):
    def test_scalar_exact(self):
        # For T_j=a*b^j, placement does not depend on the occupied slots.
        a, b = 0.9, 0.8
        envelope = UniformInputEnvelope(2, 1, 4, 1)
        local = np.array([[[a * b**j]] for j in range(3)])
        receipt = estimate(local, K=32, envelope=envelope, occupancies=[1, 2, 4],
                           tilt=Q(1, 20), windows=2)
        for q in (1, 2, 4):
            moment = 2 * (4 * log(a) + q * log((1 + 15*b)/16))
            self.assertAlmostEqual(receipt['witnesses'][str(q)]['log_moment'], moment, places=11)

    def test_matrix_state_not_reset(self):
        matrix = np.array([[0.8, 0.1], [0.02, 0.6]])
        local = np.stack([matrix, matrix, matrix])
        envelope = UniformInputEnvelope(2, 1, 4, 1)
        result = estimate(local, K=16, envelope=envelope, occupancies=[1],
                          tilt=Q(1, 20), windows=2)
        expected = log(np.linalg.matrix_power(matrix, 4)[0].sum())
        reset = 2 * log(np.linalg.matrix_power(matrix, 2)[0].sum())
        self.assertAlmostEqual(result['witnesses']['1']['log_moment'], expected, places=12)
        self.assertGreater(abs(expected-reset), 0.001)

    def test_doubled_beta_identity(self):
        small, large = UniformInputEnvelope(16, 8, 4, 4), UniformInputEnvelope(16, 8, 4, 8)
        self.assertEqual(large.beta / small.beta**2, Q(65535, 65537)**8)
        self.assertEqual((large.message_bits, large.output_bits, large.regions), (256, 512, 128))

    def test_upper_conversion_and_invalid_geometry(self):
        a = upper_arrays([arb_mat([[1, 0], [arb(1)/3, 1]])])
        self.assertGreater(a[0, 1, 0], 0)
        with self.assertRaises(ValueError):
            estimate(np.ones((33, 1, 1)), K=4097,
                     envelope=UniformInputEnvelope(16, 8, 4, 8),
                     occupancies=[1], tilt=Q(1, 10))


if __name__ == '__main__':
    unittest.main()

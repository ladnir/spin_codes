import unittest
import numpy as np

import hd_field128_gamma_transfer_diagnostic as diagnostic
from hd_randomized_outer_state_profile import multiply, power


class TransferTests(unittest.TestCase):
    def test_batched_power_matches_scalar_log_oracle(self):
        rng = np.random.default_rng(92846)
        matrices = -rng.uniform(0, 100, (3, 4, 2, 2))
        for exponent in (0, 1, 2, 7, 128):
            actual = diagnostic.power(matrices, exponent)
            for i in range(3):
                for j in range(4):
                    expected = np.asarray(power(matrices[i, j], exponent))
                    np.testing.assert_allclose(actual[i, j], expected, atol=1e-10, rtol=0)

    def test_common_contraction_is_an_upper_bound(self):
        rng = np.random.default_rng(89173)
        for weight in (1, 2, 4, 8):
            full = diagnostic.tail.full_counts(4, weight)
            regions = -rng.uniform(0, 20, (len(full), 2, 2))
            row = diagnostic.support_sum(4, weight, regions, repetitions=8)
            for ratio in (-20., -3., 0., 7., 20.):
                upper = diagnostic.tail.support_sum_logs(
                    4, weight, regions, ratio, method='runs', repetitions=8)
                self.assertLessEqual(row['log2_support_sum_screen'], upper+1e-10)

    def test_independent_small_support_sum_and_chunk_invariance(self):
        pairs, weight = 4, 3
        regions = np.array([[[-1., -3.], [-5., -2.]], [[-3., -2.], [-4., -1.]]])
        # Cyclic length8, weight3: c=1 =>8; c=2 =>32; c=3 =>16 supports.
        counts = {1: 8, 2: 32, 3: 16}
        terms = []
        for even in range(2):
            for odd in range(2):
                product = power(multiply(regions[even], regions[odd]), 4)
                terms.append(np.log2(counts[weight-even-odd])+np.logaddexp2(*product[0]))
        expected = np.logaddexp2.reduce(terms)
        for chunk in (1, 2, 32):
            actual = diagnostic.support_sum(pairs, weight, regions, repetitions=8, chunk=chunk)
            self.assertAlmostEqual(actual['log2_support_sum_screen'], expected, places=10)


if __name__ == '__main__':
    unittest.main()

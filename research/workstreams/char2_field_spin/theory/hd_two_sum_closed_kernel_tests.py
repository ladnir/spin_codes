"""Cross-check the closed kernel against independent finite-field ranks."""
from itertools import combinations
from math import comb
import random
import unittest

from gf8_packet_maps import rank, MUL, xor_sum
import hd_two_sum_affine_outer as old
import hd_two_sum_gf16_sparse as large
from hd_two_sum_closed_kernel import (uniform_mask_count, uniform_mask_counts,
    arity_power_sum, arity_holder_moments, zero_one_moment, zero_two_moments)


class ClosedKernelTests(unittest.TestCase):
    def test_entire_g8_kernel_against_rank_and_mobius(self):
        for active in range(256):
            self.assertEqual(uniform_mask_counts(active), old.uniform_mask_counts(active))

    def test_all_g8_arity_moments_against_fixed_row_sets(self):
        counts = [old.uniform_mask_counts(active) for active in range(256)]
        for power in (1, 2, 32):
            for a in range(9):
                for target in range(256):
                    expected = sum(counts[active][target]**power for active in range(256) if active.bit_count() == a)
                    self.assertEqual(arity_power_sum(a, target.bit_count(), columns=power), expected)

    def test_g16_sampled_containment_ranks(self):
        rng = random.Random(166516)
        h = large.matrix()
        # The containment count sums the proposed exact mask kernel; the
        # independent oracle uses an arbitrary rectangular GF256 submatrix.
        for a in (0, 1, 2, 3, 4, 8, 12, 15, 16):
            active = sum(1 << i for i in rng.sample(range(16), a))
            columns = [i for i in range(16) if active >> i & 1]
            targets = [0, 65535]+[sum(1 << i for i in rng.sample(range(16), size)) for size in (3, 7, 11, 14)]
            for target in targets:
                expected = 256**(a-rank([[h[i][j] for j in columns] for i in range(16) if not target >> i & 1]))
                value, subset = 0, target
                while True:
                    value += uniform_mask_count(active, subset, 16)
                    if subset == 0:
                        break
                    subset = (subset-1) & target
                self.assertEqual(value, expected)

    def test_g16_arity_counts_by_two_independent_enumerations(self):
        # For representative targets enumerate all A, not just sampled A.
        for target in (0, 1, 0x5555, 0xfffe, 0xffff):
            sums = [[0]*17 for _ in range(2)]
            for active in range(65536):
                value = uniform_mask_count(active, target, 16)
                sums[0][active.bit_count()] += value
                sums[1][active.bit_count()] += value*value
            for a in range(17):
                for power in (1, 2):
                    self.assertEqual(sums[power-1][a], arity_power_sum(a, target.bit_count(), 16, power))
        for a in range(17):
            self.assertEqual(sum(comb(16, w)*arity_power_sum(a, w, 16, 1) for w in range(17)), comb(16, a)*256**a)

    def test_every_principal_determinant_and_bordered_coefficient(self):
        # Closed2x2 determinants and all outside/inside coefficients; a
        # zero coefficient would invalidate the one-equation NZ count.
        for labels in (old.LABELS, large.LABELS):
            g = len(labels)
            for active in range(1 << g):
                a = active.bit_count()
                u = xor_sum(labels[i] for i in range(g) if active >> i & 1)
                det = MUL[u][u] if a % 2 else 1
                self.assertNotEqual(det, 0)
                for j in range(g):
                    if active >> j & 1:
                        continue
                    alpha = 1 ^ MUL[u][u] ^ MUL[labels[j]][u]
                    beta = u ^ (labels[j] if not a % 2 else 0)
                    for i in range(g):
                        if active >> i & 1:
                            self.assertNotEqual(alpha ^ MUL[beta][labels[i]], 0)

    def test_holder_normalizer_and_small_counts(self):
        self.assertEqual([zero_one_moment(j) for j in range(3)], [1, 0, 255])
        self.assertEqual([zero_two_moments(j) for j in range(4)], [1, 0, 0, 255])
        for a in range(3, 17):
            values, denominator = arity_holder_moments(a, 16, 16, 8)
            self.assertEqual(denominator, 255**(8*a))
            self.assertEqual(values[0], comb(16, a))


if __name__ == '__main__':
    unittest.main()

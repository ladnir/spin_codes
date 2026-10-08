"""Small exact tests of the regional rank-envelope proof object."""

from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest

from hd_two_state_regional import (exact_bound, identity, matmul, polynomial_power,
    post_profile, post_transfers, precoder_profiles, regional_post_kernels,
    sharp_factor, source_rows, zero_matrix)


def gf4_mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


SCALE = [[gf4_mul(a, z & 3) | (gf4_mul(a, z >> 2) << 2)
          for z in range(16)] for a in range(4)]


class RegionalRankTests(unittest.TestCase):
    def test_sharp_live_pivot_lemma_all_gf4_forms(self):
        for tau in (F(1), F(2), F(3), F(4), F(7)):
            for live in range(1, 16):
                for other in range(16):
                    expected = F(0)
                    for a in range(1, 4):
                        equation = other ^ SCALE[a][live]
                        span = {SCALE[c][equation] for c in range(4)}
                        if live in span:
                            continue
                        expected += (F(1) if not equation else 1 / tau) / 3
                    self.assertLessEqual(expected, sharp_factor(tau, 4))

    def test_region_kernel_is_average_of_only_regional_masks(self):
        t0, t1 = post_transfers(F(2, 3), F(1, 5))
        for rows in range(1, 5):
            kernels = regional_post_kernels(rows, F(2, 3), F(1, 5))
            for h in range(rows + 1):
                average = zero_matrix()
                for subset in combinations(range(rows), h):
                    selected = set(subset)
                    matrix = identity()
                    for i in range(rows):
                        matrix = matmul(matrix, t1 if i in selected else t0)
                    for a, b in product(range(4), repeat=2):
                        average[a][b] += matrix[a][b] / comb(rows, h)
                self.assertEqual(kernels[h], average)

    def test_equal_total_occupancy_does_not_determine_post_bound(self):
        kernels = regional_post_kernels(2, F(1, 2), F(1, 3))
        separated = post_profile((2, 0), kernels)
        balanced = post_profile((1, 1), kernels)
        reversed_order = post_profile((0, 2), kernels)
        self.assertNotEqual(separated, balanced)
        self.assertNotEqual(separated, reversed_order)
        global_kernel = regional_post_kernels(4, F(1, 2), F(1, 3))
        global_value = post_profile((2,), global_kernel)
        # Uniform global masks mix regional occupancies hypergeometrically.
        expected = sum(F(comb(2, h) * comb(2, 2 - h), comb(4, 2))
                       * post_profile((h, 2 - h), kernels) for h in range(3))
        self.assertEqual(global_value, expected)

    def test_source_supports_use_shortening_not_message_counts(self):
        rows = dict(source_rows(4, 2))
        self.assertEqual(len(rows), 6)
        self.assertEqual(rows[0], 0)
        for mask, dimension in rows.items():
            if mask:
                self.assertEqual(dimension, mask.bit_count() - 2)
                self.assertIn(mask.bit_count(), (3, 4))

    def test_pre_global_carry_and_zero_source_exclusion(self):
        tau, g = F(3), F(1, 3)
        profiles = precoder_profiles(2, 1, 1, tau, g)
        self.assertEqual(profiles, {(1, 0): tau * g, (1, 1): tau * (1 + g)})
        self.assertNotIn((0, 0), profiles)

    def test_sharper_bound_is_monotone_in_return_factor(self):
        tau, q, y = F(2), 4, F(1, 2)
        old = F(1, q - 1) + (1 - F(1, q - 1)) / tau
        sharp = sharp_factor(tau, q)
        old_profiles = precoder_profiles(2, 1, 2, tau, old)
        new_profiles = precoder_profiles(2, 1, 2, tau, sharp)
        old_kernels = regional_post_kernels(2, y, old)
        new_kernels = regional_post_kernels(2, y, sharp)
        old_moment = sum(mass * post_profile(h, old_kernels) for h, mass in old_profiles.items())
        new_moment = sum(mass * post_profile(h, new_kernels) for h, mass in new_profiles.items())
        self.assertLess(new_moment, old_moment)

    def test_exact_small_complete_contraction(self):
        result = exact_bound()
        self.assertIsInstance(result["rank_moment"], F)
        self.assertEqual(result["g"], F(1, 3))
        self.assertGreater(result["profile_count"], 1)
        self.assertEqual(result["route"], "independent uniform regional")
        self.assertGreater(result["probability_upper_unclipped"], 1)
        self.assertFalse(result["whole_code_certified"])


if __name__ == "__main__":
    unittest.main()

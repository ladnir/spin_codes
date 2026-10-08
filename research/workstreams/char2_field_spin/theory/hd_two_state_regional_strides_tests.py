"""Exact and directed checks of the fixed affine-placement receipt."""

from fractions import Fraction as F
from itertools import product
from math import comb
import unittest

from hd_two_state_regional_strides import (ROWS, STRIDES, cancellation_tail_receipt,
    certify_spectrum, safe_column_bound)
from hd_two_state_regional_tests import gf4_mul


class StrideReceiptTests(unittest.TestCase):
    def test_fixed_affine_addresses_are_bijections(self):
        self.assertEqual(len(set(STRIDES)), 64)
        for stride in STRIDES:
            inverse = pow(stride, -1, ROWS)
            self.assertEqual(stride * inverse % ROWS, 1)
            self.assertEqual(len({stride * j % ROWS for j in range(ROWS)}), ROWS)

    def test_directed_all_mode_receipt(self):
        receipt = certify_spectrum()
        self.assertTrue(receipt["rigorous_spectral_inequality"])
        self.assertEqual(receipt["modes_checked"], ROWS // 2)
        self.assertEqual(receipt["lambda_upper"], "323/1000")

    def test_small_all_subsets_internal_edge_bound(self):
        # The symmetric generators +/-1,+/-3 on Z/8 give K4,4.
        # Its largest nonconstant positive eigenvalue is exactly zero.
        for bits in product((0, 1), repeat=8):
            alpha = F(sum(bits), 8)
            edges = F(sum(bits[x] * bits[(x + a) % 8]
                          for x in range(8) for a in (1, 3)), 16)
            triples = F(sum(bits[x] * bits[(x + a) % 8] * bits[(x + 2 * a) % 8]
                            for x in range(8) for a in (1, 3)), 16)
            self.assertLessEqual(triples, edges)
            self.assertLessEqual(edges, alpha * alpha)

    def test_adjacent_equality_events_are_independent(self):
        alphabet, length = 3, 5
        words = list(product(range(alphabet), repeat=length))
        for selected in product((0, 1), repeat=length - 1):
            probability = F(sum(all(not selected[i] or word[i] == word[i + 1]
                                    for i in range(length - 1)) for word in words),
                            len(words))
            self.assertEqual(probability, F(1, alphabet ** sum(selected)))

    def test_middle_inverse_entry_cancellation_identity(self):
        inverse = {a: next(b for b in range(1, 4) if gf4_mul(a, b) == 1)
                   for a in range(1, 4)}
        for first, second_left, second_right in product(range(1, 4), repeat=3):
            middle = gf4_mul(inverse[first], inverse[second_left] ^ inverse[second_right])
            self.assertEqual(middle == 0, second_left == second_right)

    def test_quarter_cut_and_exact_cancellation_tail(self):
        two = safe_column_bound(cancelled_middle=1023)
        one = safe_column_bound(pre_order=1)
        self.assertEqual(two["safe_column_upper"], 259231)
        self.assertEqual(one["safe_column_upper"], 258144)
        self.assertEqual(one["check_capacity"], 262144)
        self.assertTrue(one["fewer_safe_columns_than_checks"])
        self.assertTrue(two["fewer_safe_columns_than_checks"])
        receipt = cancellation_tail_receipt()
        self.assertEqual(receipt["failure_upper_strictly_below_power_of_two"], -3649)
        self.assertLess(comb((1 << 21) - 1, 1024) << 3649, 65535 ** 1024)


if __name__ == "__main__":
    unittest.main()

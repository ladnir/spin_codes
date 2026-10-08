"""Exact support-law checks and the directed sparse quiet-codeword gate."""

from fractions import Fraction as F
from itertools import combinations
from math import comb
import unittest

import hd_blockwise_shared_cuts as proof


class SharedCutTests(unittest.TestCase):
    def test_exact_regional_generating_function(self):
        for rows, width, threshold in ((8, 4, 3), (12, 4, 3), (8, 2, 2)):
            for active in range(rows+1):
                for mark in (F(1), F(3, 2), F(2)):
                    enumerated = sum((mark**proof.crowded_count(subset, rows, width, threshold)
                                      for subset in combinations(range(rows), active)), F(0))/comb(rows, active)
                    self.assertEqual(enumerated, proof.regional_moment_exact(
                        rows, active, mark, width, threshold))

    def test_coefficient_mark_is_an_upper_bound(self):
        rows, width, threshold = 12, 4, 3
        for active in range(1, rows+1):
            for mark in (F(3, 2), F(2)):
                x = F(active, 2*rows)
                polynomial = proof.packet_polynomial(mark, width, threshold)
                g = sum((coefficient*x**j for j, coefficient in enumerate(polynomial)), F(0))
                bound = g**(rows//width)/(x**active*comb(rows, active))
                self.assertGreaterEqual(bound, proof.regional_moment_exact(
                    rows, active, mark, width, threshold))

    def test_actual_small_field_quiet_words(self):
        receipt = proof.small_quiet_witnesses()
        self.assertGreater(receipt['nonzero_quiet_words_checked'], 0)
        self.assertTrue(all(active >= 3 for active in receipt['active_row_histogram']))
        self.assertTrue(receipt['terminal_free'])

    def test_fixed_marks_cover_every_claimed_class(self):
        for active in range(5, 257):
            x, mark = proof.fixed_marks(active)
            self.assertEqual(x, F(active, 16384))
            self.assertGreaterEqual(mark, 2)
        with self.assertRaises(ValueError):
            proof.fixed_marks(4)
        with self.assertRaises(ValueError):
            proof.fixed_marks(257)

    def test_directed_all_sparse_quiet_classes(self):
        receipt = proof.sparse_quiet_receipt()
        self.assertEqual(receipt['class_count'], 252)
        self.assertEqual(receipt['worst_class']['active_rows'], 8)
        self.assertIn('-2236.709505', receipt['quiet_event_log2_upper'])
        self.assertTrue(receipt['below_2_to_minus_40'])
        self.assertTrue(receipt['output_weight_unrestricted'])
        self.assertTrue(receipt['feedthrough_independent'])
        self.assertFalse(receipt['whole_code_certificate'])

    def test_joint_tile_support_and_feedback_rank(self):
        for offset in (0, 8, 24, 56):
            for w in range(9):
                for u in range(9):
                    source = tuple((j+3) % 8 for j in range(w))
                    target = tuple((j+5) % 8 for j in range(u))
                    receipt = proof.tile_support_rank(source, target, offset)
                    dimension = max(0, w+u-8)
                    self.assertEqual(receipt['support_space_dimension'], dimension)
                    self.assertEqual(receipt['feedback_rank_on_support_space'], min(4, dimension))
                    self.assertEqual(receipt['quiet_tile_dimension'], max(0, w+u-12))

    def test_quiet_packet_branch_bound_is_attained(self):
        receipt = proof.quiet_tile_branch_witness()
        self.assertEqual(receipt['input_weight'], 5)
        self.assertEqual(receipt['output_weight'], 8)
        self.assertEqual(receipt['combined_weight'], 13)
        self.assertTrue(receipt['feedback_zero'])

    def test_base_field_descent(self):
        receipt = proof.small_field_descent()
        self.assertEqual(receipt['quiet_dimension'], 2)
        self.assertEqual(receipt['nonzero_extension_words_checked'], 255)
        self.assertEqual(receipt['base_minimum'], 2)
        self.assertEqual(receipt['extension_minimum'], 2)

    def test_GF256_quiet_first_moment_still_fails(self):
        receipt = proof.base_field_quiet_moment()
        self.assertEqual(receipt['counting_field'], 256)
        self.assertEqual(receipt['quiet_dimension_lower'], 29496)
        self.assertEqual(receipt['feedback_equations'], 52424)
        self.assertEqual(receipt['output_weight_upper'], 838848)
        self.assertIn('4110.59915438587', receipt['quiet_projective_moment_log2_lower'])
        self.assertTrue(receipt['moment_exceeds_one'])
        self.assertTrue(receipt['substantial_failure_probability_not_established'])


if __name__ == '__main__':
    unittest.main()

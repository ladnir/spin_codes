"""Small exact type/count checks and the directed high-type quiet class."""

from fractions import Fraction as F
from itertools import product
import unittest

import hd_blockwise_shared_quiet_compose as proof


class QuietComposeTests(unittest.TestCase):
    def test_feedback_weight_enumerator(self):
        from hd_blockwise_inner import dot
        check = ((1, 1, 1, 1), (0, 1, 2, 3))
        histogram = [0]*5
        for vector in product(range(4), repeat=4):
            if not any(dot(row, vector, 4) for row in check):
                histogram[sum(value != 0 for value in vector)] += 1
        self.assertEqual(tuple(histogram), proof.quiet_weight_polynomial(4, 2, 4))
        large = proof.quiet_weight_polynomial()
        self.assertEqual(sum(large), 256**60)
        self.assertEqual(large[1:5], (0, 0, 0, 0))

    def test_most_concentrated_capped_type(self):
        for rows in range(1, 9):
            for maximum in range(1, rows+1):
                alphabet = 4
                if maximum*alphabet < rows:
                    continue
                proposed = proof.orbit_size(proof.worst_capped_counts(rows, maximum, alphabet))
                for counts in product(range(maximum+1), repeat=alphabet):
                    if sum(counts) == rows:
                        self.assertLessEqual(proposed, proof.orbit_size(counts))

    def test_exact_information_set_membership_bound(self):
        receipt = proof.small_information_set_check()
        self.assertEqual(receipt['fixed_inverse_words_checked'], 4096)
        self.assertGreater(receipt['words_with_membership_probability_one'], 0)
        self.assertGreater(receipt['words_with_zero_membership_probability'], 0)

    def test_dense_support_relaxation_witness(self):
        for translation in (8, 64):
            receipt = proof.dense_input_sparse_output_witness(translation)
            self.assertEqual(receipt['input_weight'], 64)
            self.assertEqual(receipt['output_weight'], 8)
            self.assertEqual(receipt['quiet_support_space_dimension'], 4)
            self.assertGreater(receipt['full_output_support_vectors_per_packet'], 0)
            self.assertTrue(receipt['outer_value_constraints_not_imposed'])
            self.assertFalse(receipt['actual_code_counterexample'])

    def test_directed_high_type_class(self):
        receipt = proof.high_type_receipt()
        self.assertIn('7693823.1287572', receipt['quiet_projective_ball_log2_upper'])
        self.assertEqual(receipt['sufficient_top128_information_bits'], 7693904)
        self.assertIn('-25534.654265', receipt['bounded_multiplicity_class_log2_upper'])
        self.assertTrue(receipt['bounded_multiplicity_class_below_2_to_minus_40'])
        self.assertTrue(receipt['covers_all_active_row_counts'])
        self.assertFalse(receipt['whole_code_certificate'])

    def test_crude_multiplicity48_is_not_a_proof_gate(self):
        receipt = proof.high_type_receipt(maximum=48, failure_bits=40)
        self.assertEqual(receipt['sufficient_top128_information_bits'], 7693864)
        self.assertIn('5352.411435', receipt['bounded_multiplicity_class_log2_upper'])
        self.assertFalse(receipt['bounded_multiplicity_class_below_2_to_minus_40'])


if __name__ == '__main__':
    unittest.main()

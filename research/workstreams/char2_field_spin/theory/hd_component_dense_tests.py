"""Exact small checks and directed large receipts for the dense proof branch."""

from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest

import hd_component_dense as proof
from hd_affine_proof import mul


class DenseProofTests(unittest.TestCase):
    def test_native_placement_order_is_inverse_of_store_address(self):
        rows = 64
        for group in (1, 2, 4):
            for stride in proof.STRIDES:
                order = proof.placement_order(rows, stride, group)
                self.assertEqual(set(order), set(range(rows)))
                inv = pow(stride, -1, rows//group)
                for position, row in enumerate(order):
                    placed = group*((inv*(row//group)) % (rows//group))+row % group
                    self.assertEqual(placed, position)

    def test_unit_basis_differences_have_only_permitted_source_rows(self):
        rows, permitted = 64, 38
        mask = proof.fixed_mask(rows, permitted)
        for group in (1, 4):
            selected = proof.unit_positions(mask, (1, 3, 5, 7), group)
            for c, stride in enumerate((1, 3, 5, 7)):
                order = proof.placement_order(rows, stride, group)
                for position in selected[c]:
                    self.assertLess(position+1, rows)
                    self.assertTrue(mask[order[position]])
                    self.assertTrue(mask[order[position+1]])

    def test_uniform_regional_containment_probability(self):
        for rows in range(2, 9):
            for window in range(1, rows+1):
                for occupied in range(window+1):
                    masks = list(combinations(range(rows), occupied))
                    favorable = sum(all(position < window for position in mask) for mask in masks)
                    self.assertEqual(F(favorable, len(masks)), F(comb(window, occupied), comb(rows, occupied)))

    def test_actual_gf4_codewords_with_both_post_laws(self):
        for group, three_product in product((1, 4), (False, True)):
            receipt = proof.small_actual_code_check(group, three_product)
            self.assertGreaterEqual(receipt['actual_dimension'], receipt['dimension_lower'])
            self.assertGreater(receipt['checked_basis_words'], 0)
            self.assertEqual(receipt['output_weight_upper'], 100)

    def test_information_set_weight_ball(self):
        # An injective [4,2] GF4 code with its first two entries systematic.
        words = [(a, b, a ^ b, mul(2, a) ^ b) for a, b in product(range(4), repeat=2)]
        for cutoff in (0, 1, 2):
            actual = sum(0 < sum(x != 0 for x in word) <= cutoff for word in words)
            upper = sum(comb(2, weight)*3**weight for weight in range(1, cutoff+1))
            self.assertLessEqual(actual, upper)

    def test_dense_geometric_ratio_bound_exactly(self):
        rows, b, rho, Q, cutoff, start = 16, 4, 2, 3, 6, 10
        first_dimension = b*start
        ratio = F(rows-start, start+1)*F(1, Q**rho)*F(first_dimension+1, first_dimension-cutoff+1)**b
        terms = [F(comb(rows, a)*comb(b*a, cutoff), Q**(rho*a)) for a in range(start, rows+1)]
        self.assertLess(ratio, 1)
        for left, right in zip(terms, terms[1:]):
            self.assertLessEqual(right, ratio*left)
        self.assertLessEqual(sum(terms), terms[0]/(1-ratio))

    def test_entire_k20_dense_tail(self):
        receipt = proof.dense_tail_upper()
        self.assertEqual(receipt['number_of_classes'], 2826)
        self.assertTrue(receipt['entire_tail_below_2_to_minus_80'])
        self.assertIn('-322.5829159528', receipt['entire_tail_log2_upper'])
        self.assertFalse(receipt['whole_distance_certificate'])

    def test_both_k20_actual_first_moment_obstructions(self):
        for group, variables, dimension in ((1, 754354, 125106), (4, 756089, 126841)):
            receipt = proof.obstruction_receipt(group)
            self.assertEqual(receipt['unit_variables'], variables)
            self.assertEqual(receipt['actual_code_dimension_lower'], dimension)
            self.assertLess(receipt['output_weight_upper'], receipt['forbidden_cutoff'])
            self.assertTrue(receipt['rules_out_uncapped_first_moment'])
            self.assertFalse(receipt['rejects_40_bit_distance'])


if __name__ == '__main__':
    unittest.main()

"""Integer identities for the fixed-input-set joint-weight kernel."""
from fractions import Fraction as F
from math import comb
import unittest

import numpy as np

import hd_two_sum_input_weight as target
from hd_two_sum_closed_kernel import uniform_mask_count
from hd_two_sum_gf16_fractional_state import fractional_dense_weights


class InputWeightTests(unittest.TestCase):
    def test_ambient_matches_labelled_kernel(self):
        for group in (8, 16):
            for arity in range(group+1):
                active = (1 << arity)-1
                for weight in range(group+1):
                    for overlap in range(max(0, arity+weight-group), min(arity, weight)+1):
                        mask = (1 << overlap)-1 | (((1 << (weight-overlap))-1) << arity)
                        self.assertEqual(target.ambient_count(arity, overlap, weight, group),
                                         uniform_mask_count(active, mask, group))

    def test_positive_coefficients_normalize_by_input_weight(self):
        for group in (8, 16):
            for arity in range(group+1):
                totals = [0]*(arity+1)
                for weight in range(group+1):
                    for overlap in range(max(0, arity+weight-group), min(arity, weight)+1):
                        multiplicity = comb(arity, overlap)*comb(group-arity, weight-overlap)
                        for k, count in enumerate(target.input_weight_counts(arity, overlap, weight, group)):
                            self.assertGreaterEqual(count, 0)
                            totals[k] += multiplicity*count
                self.assertEqual(totals, [comb(arity, k)*255**k for k in range(arity+1)])

    def test_marker_one_recovers_existing_dense_holder(self):
        for alpha in (F(1), F(1, 5)):
            actual = target.type_weights(alpha, 0.)
            expected = fractional_dense_weights(alpha)
            for arity in range(3, 17):
                np.testing.assert_allclose(actual[arity], expected[arity], atol=2e-13, rtol=0)

    def test_zero_output_requires_zero_input(self):
        for arity in range(17):
            self.assertEqual(target.input_weight_counts(arity, 0, 0), (1,)+(0,)*arity)
            self.assertEqual(target.marked_count(arity, 0, 0, 256), 1)

    def test_involution_joint_symmetry(self):
        # H^2=I for these affine labels, so exchanging exact input/output
        # supports preserves the number of assignments for fixed supports.
        for arity in range(17):
            for weight in range(17):
                for overlap in range(max(0, arity+weight-16), min(arity, weight)+1):
                    self.assertEqual(target.exact_support_count(arity, overlap, weight),
                                     target.exact_support_count(weight, overlap, arity))

    def test_marked_mass_identity(self):
        for arity in (1, 3, 9, 16):
            for marker in (F(1), F(3, 2), F(4)):
                total = sum(comb(arity, overlap)*comb(16-arity, weight-overlap)
                            *target.marked_count(arity, overlap, weight, marker)
                            for weight in range(17)
                            for overlap in range(max(0, arity+weight-16), min(arity, weight)+1))
                self.assertEqual(total, (1+255*marker)**arity)


if __name__ == '__main__':
    unittest.main()

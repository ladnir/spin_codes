"""Light exact controls for the one-group dense specialization."""
from fractions import Fraction as F
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

from hd_two_sum_gf16_single_dense import WITNESSES, single_source_tables, arity_tail, run
from hd_two_sum_gf16_small_cover import WITNESSES as OLD_WITNESSES
from hd_two_sum_gf16_arity_dense import weight_vectors


class SingleDenseTests(unittest.TestCase):
    def test_exact_same_three_frozen_witnesses(self):
        self.assertEqual(WITNESSES, OLD_WITNESSES)

    def test_single_source_requires_only_empty_or_singleton_geometry(self):
        raw = {(0, 0): np.array([[0., -np.inf], [-2., -3.]]),
               (1, 0): np.array([[-1., -4.], [-5., -6.]])}
        actual = single_source_tables(raw)
        self.assertEqual(actual.shape, (16, 2, 2, 2))
        for region in actual:
            np.testing.assert_array_equal(region, np.array([raw[0, 0], raw[1, 0]]))

    def test_scalar_tail_is_exact_binomial_positive_power(self):
        # One empty mask and sixteen equally weighted singleton masks;
        # supplied weight vectors have already summed those16 choices.
        vectors = np.full((1, 17, 2), -np.inf)
        vectors[0, 0] = np.log2([2., 1.])
        vectors[0, 1] = np.log2([16.*3, 16.*F(3, 2)])
        weights = np.zeros(17)
        bound, _ = arity_tail(vectors, -1., weights)
        exact = sum(comb(16, j)*2**(16-j)*48**j for j in range(9, 17))
        self.assertAlmostEqual(bound, log2(exact), places=11)

    def test_zero_sparse_dp_identity_counts_masks_once(self):
        raw = {(m, 0): np.array([[0., -np.inf], [-np.inf, 0.]]) for m in (0, 1)}
        vectors = weight_vectors(single_source_tables(raw), -1., sparse_groups=0)
        expected = np.array([[comb(16, w), comb(16, w)/2] for w in range(17)])
        np.testing.assert_allclose(np.exp2(vectors[0]), expected, rtol=2e-14)

    def test_selection_and_scope_do_not_multiply_by_witness_count(self):
        tail = lambda *args: (0., [dict(log2_mass=0., nonempty_dense_columns=9)])
        with patch('hd_two_sum_gf16_single_dense.pair_regions_screen', return_value={}) as bank, \
                patch('hd_two_sum_gf16_single_dense.single_source_tables', return_value=None), \
                patch('hd_two_sum_gf16_single_dense.weight_vectors', return_value=None), \
                patch('hd_two_sum_gf16_single_dense.arity_tail', side_effect=tail):
            result = run()
        self.assertEqual(bank.call_count, 3)
        expected = min(log2(8192)-log2(255)-734003*log2(float(z))-float(p) for _, z, p in WITNESSES)
        self.assertAlmostEqual(result['failure_log2_screen'], expected+log2(14), places=11)
        self.assertEqual(sum(row['fixed_row_sets'] for row in result['entries']), 65399)
        self.assertEqual(result['cohort_counts'], [0, 0, 1])
        self.assertEqual(result['source_supports'], 8192)
        self.assertEqual(result['sparse_normalizer_exponent'], 0)
        self.assertFalse(result['directed'])

    def test_budget_rejects_expansion_and_times_out_before_bank(self):
        for invalid in (True, 0, 21, float('inf')):
            with self.assertRaises(ValueError):
                run(max_seconds=invalid)
        with patch('hd_two_sum_gf16_single_dense.pair_regions_screen') as bank:
            result = run(max_seconds=1e-12)
        bank.assert_not_called()
        self.assertFalse(result['complete_selected_cohort'])
        self.assertNotIn('failure_log2_screen', result)


if __name__ == '__main__':
    unittest.main()

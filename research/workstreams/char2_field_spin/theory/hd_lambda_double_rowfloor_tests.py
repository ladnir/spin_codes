"""Tiny algebra and mocked controls only; no regional or proof-gate run."""
from itertools import product
from math import log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_double_rowfloor_control as target
from hd_lambda_double_middle_tests import kernel, apply


class OriginalRowFloorTests(unittest.TestCase):
    def test_all_scheduled_quiet_supports_and_original_row_mapping(self):
        h = target.base.fiber.outer.matrix(2, 4)
        for j in range(16):
            data = target.base.fiber.column_fiber_counts(target.SOURCE, j, q=4)
            self.assertTrue(data['all_zero_fiber_coordinates_protected'])
            mapped = tuple(data['mapped_source'])
            outside = tuple(i for i in range(16) if i not in mapped)
            basis = kernel([[h[i][k] for k in mapped] for i in outside], 6)
            self.assertEqual(len(basis), int(j < 8))
            for vector in basis:
                x = tuple(vector[mapped.index(i)] if i in mapped else 0 for i in range(16))
                y = apply(h, x)
                self.assertEqual(y, x)
                raw = sum(bool(value) << i for i, value in enumerate(y))
                physical = target.base.fiber._permute_mask(raw, data['output_permutation'])
                original = {i for i in target.SOURCE if x[data['input_permutation'][i]]}
                self.assertEqual(original, set(target.SOURCE) if j in (0, 1, 6, 7) else {2, 3, 4, 5})
                self.assertEqual(target.original_score(physical, data), int(0 in original))
                self.assertEqual(target.original_score(physical, data, 1), int(1 in original))
            self.assertEqual(target.original_score(0, data), 0)
            self.assertEqual(data['zero_fiber_counts'][0], 1)

    def test_refuses_unprotected_or_foreign_source(self):
        data = target.base.fiber.column_fiber_counts(target.SOURCE, 0, q=4)
        with self.assertRaises(ValueError):
            target.original_score(0, dict(data, all_zero_fiber_coordinates_protected=False))
        with self.assertRaises(ValueError):
            target.original_score(0, dict(data, original_source=tuple(range(7))))

    def test_degree_halving_keeps_exact_noncommuting_tail(self):
        # Three chronological matrices, two moment replicas; sparse even degrees.
        low = np.log2(np.arange(1, 3*2*4+1).reshape(3, 2, 2, 2))
        doubled = np.full((3, 3, 2, 2), -np.inf)
        doubled[:, ::2] = low
        for floor in range(7):
            self.assertAlmostEqual(float(target.base.positive_deficit_tail(low, 2, floor)),
                float(target.base.positive_deficit_tail(doubled, 2, 2*floor)), places=12)

    def test_score_inflation_preserves_tail_under_positive_products(self):
        # Two-column toy: actual NQ row0 values may be0 or1, both mapped to1.
        branches = ((0, 0, 3), (0, 1, 5), (1, 1, 7))
        for floor in range(3):
            true_tail = inflated_tail = 0
            for left, right in product(branches, repeat=2):
                mass = left[2]*right[2]
                true_tail += mass*(left[0]+right[0] >= floor)
                inflated_tail += mass*(left[1]+right[1] >= floor)
            self.assertLessEqual(true_tail, inflated_tail)

    def test_raw_quiet_mass_keeps_nonzero_score_zero(self):
        # Fake table gives moment norm1 to every output mask. Only tiny algebra runs.
        table = np.zeros((65536, 2, 2))
        bound = dict(selected=np.full((2, 2), -np.inf), direct_terms=0,
                     factor_2r=np.full((2, 2), -np.inf), direct_union=np.full((2, 2), -np.inf))
        with patch.object(target.base, 'nonquiet_norm', return_value=bound):
            columns, _, _ = target.prepare_columns(table, table)
        for j in range(16):
            if j in (0, 1, 6, 7):
                self.assertTrue(np.allclose(columns[j, 0], 0.))
                self.assertTrue(np.allclose(columns[j, 1], log2(255)))
            elif j < 8:
                self.assertTrue(np.allclose(columns[j, 0], 8.))
                self.assertTrue(np.all(np.isneginf(columns[j, 1])))
            else:
                self.assertTrue(np.allclose(columns[j, 0], 0.))
                self.assertTrue(np.all(np.isneginf(columns[j, 1])))

    def test_frozen_baseline_validation_and_mock_normalization(self):
        baseline, old = target.load_baseline()
        columns = np.zeros((16, 2, 2, 2))
        with patch.object(target.base, 'load_measure', return_value=(baseline, baseline['common_measure'], baseline['common_prefix_log2'])), \
             patch.object(target.base, 'moment_table', return_value=None), \
             patch.object(target, 'prepare_columns', return_value=(columns, np.zeros((16, 2, 2)), [])), \
             patch.object(target.base, 'positive_deficit_tail', return_value=-1000.) as tail:
            data = target.run()
        self.assertTrue(data['complete'])
        self.assertEqual(tail.call_args.args[1:], (18, 162))
        expected_root = -1000/18-48*log2(255)
        self.assertAlmostEqual(data['result']['normalized_fixed_A_root_log2'], expected_root)
        self.assertAlmostEqual(data['result']['selected_complete_fixed_A_bound_log2'],
            baseline['common_prefix_log2']+18*min(expected_root, old['normalized_fixed_A_root_log2']))
        self.assertEqual(data['RS_factor_whole_word'], '255^(-864)')
        self.assertEqual(data['retained_score_deficit'], 126)
        self.assertFalse(data['whole_code_certificate'])
        self.assertEqual(data['regional_provider_runs'], 0)


if __name__ == '__main__':
    unittest.main()

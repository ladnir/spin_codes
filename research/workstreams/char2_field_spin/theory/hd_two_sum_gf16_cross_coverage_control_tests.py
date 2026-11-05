"""Small path/mask enumerations and mocked budget/scope checks only."""
from itertools import product
from math import log2
from unittest.mock import patch
import unittest

import numpy as np

import hd_two_sum_gf16_cross_coverage_control as gate


MEASURE = dict(log2_probability=-1., log2_atom00=-2., log2_atom_cross=-1.,
               log2_scale01=-3., log2_scale10=-4., log2_scale11=-2.)


def brute(weights, groups, measure, arity, marker):
    size = len(weights)-1
    matrices = [np.full((2, 2), -np.inf), np.full((2, 2), -np.inf)]
    p = 2**measure['log2_probability']
    for path in product((0, 1), repeat=size+1):
        edges = list(zip(path, path[1:]))
        for masks in product(range(1 << size), repeat=groups):
            mass = sum(weights[mask.bit_count()]+marker*min(mask.bit_count(), arity) for mask in masks)
            covered = True
            for i, edge in enumerate(edges):
                count = sum(mask >> i & 1 for mask in masks)
                if edge == (0, 0):
                    mass += log2(p)+count*measure['log2_atom00'] if count else 0.
                elif edge == (1, 1):
                    mass += measure['log2_scale11']
                else:
                    mass += measure['log2_scale01' if edge == (0, 1) else 'log2_scale10']
                    mass += count*measure['log2_atom_cross']
                    covered &= count > 0
            matrices[0][path[0], path[-1]] = np.logaddexp2(matrices[0][path[0], path[-1]], mass)
            if covered:
                matrices[1][path[0], path[-1]] = np.logaddexp2(matrices[1][path[0], path[-1]], mass)
    return matrices


class CoverageControlTests(unittest.TestCase):
    def test_paired_columns_against_labelled_state_and_mask_words(self):
        weights = np.log2([1., 2., 3., 4.])
        expected_u, expected_c = brute(weights, 2, MEASURE, 2, .7)
        actual_c, actual_u, costs = gate.paired_columns(weights, 2, MEASURE, arity=2, log_marker=.7)
        np.testing.assert_allclose(actual_u, expected_u, atol=2e-13)
        np.testing.assert_allclose(actual_c, expected_c, atol=2e-13)
        old, _ = gate.base.column_matrix(weights, 2, MEASURE, arity=2, log_marker=.7)
        np.testing.assert_allclose(old, actual_u, atol=2e-13)
        self.assertGreater(costs['distinct_coverage_geometries'], 0)

    def test_exact_zero_cross_when_every_group_mask_is_empty(self):
        covered, _, _ = gate.paired_columns(np.array([0., -np.inf]), 2, MEASURE, arity=1, log_marker=0.)
        self.assertEqual(covered[0, 1], -np.inf)
        self.assertEqual(covered[1, 0], -np.inf)

    def test_marker_uses_sum_of_final_states_then_adjacent_coefficients(self):
        word = np.array([[[0., 2., 1.], [0., 3., 1.]],
                         [[0., 0., 0.], [0., 0., 0.]]])
        initial = np.logaddexp2.reduce(word[0], axis=0)
        self.assertAlmostEqual(gate.base.adjacent_floor_marker(initial, 1), log2(12)-2.)

    def test_timeout_keeps_partial_metadata_without_new_provider(self):
        reference = dict(prefix_log2=0., controls=[{'mode': 'polynomial', 'log2_screen': 0.}])
        with patch.object(gate, 'read_inputs', return_value=(reference, {'source_count': '1'}, MEASURE)), \
             patch.object(gate.base, 'raw_holder_weights', return_value=np.zeros(17)), \
             patch.object(gate.base, 'column_matrix', side_effect=TimeoutError('toy')):
            result = gate.run()
        self.assertFalse(result['complete_requested_controls'])
        self.assertEqual(result['raw_regional_banks_built'], 0)
        self.assertEqual(result['fixed_marker_count'], 1)
        self.assertEqual(result['controls'], [])

    def test_resource_guards_before_loading(self):
        with patch.object(gate, 'read_inputs', side_effect=AssertionError('must not read')):
            for invalid in (True, 0, 31, float('nan')):
                with self.assertRaises(ValueError):
                    gate.run(max_seconds=invalid)


if __name__ == '__main__':
    unittest.main()

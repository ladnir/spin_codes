"""All21 FLOAT wrapper scope/count tests; no production operator builds."""
from copy import deepcopy
from math import comb
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_affine_dense_two_types as gate


class AffineDenseTwoTypesTests(unittest.TestCase):
    def test_exact_ordered_assignment_partition(self):
        rows = gate.type_partition()
        represented = set()
        for row in rows:
            a, b = row['dense_arities']
            self.assertEqual(row['dense_type_multiplicity'], 1 if a == b else 2)
            self.assertEqual(row['ordered_fixed_dense_rowset_assignments'],
                             (1 if a == b else 2)*comb(8, a)*comb(8, b))
            for first in row['fixed_row_sets_left']:
                for second in row['fixed_row_sets_right']:
                    self.assertNotIn((first, second), represented)
                    represented.add((first, second))
                    if a != b:
                        self.assertNotIn((second, first), represented)
                        represented.add((second, first))
        self.assertEqual(represented, {(a, b) for a in range(256) for b in range(256)
                                     if a.bit_count() >= 3 and b.bit_count() >= 3})

    def test_pinned_literal_and_boundary(self):
        source = gate.frozen_witness()
        self.assertEqual(source['winning_literal']['mark'], str(gate.Z))
        self.assertEqual(source['winning_literal']['log2_potential'], '-40')
        self.assertTrue(source['winning_witness_hits_boundary'])
        with patch.object(gate, 'FIT_SHA256', '0'*64), self.assertRaises(ValueError):
            gate.frozen_witness()

    def test_one_shared_bank_and_multiplicity_after_tail(self):
        source = deepcopy(gate.frozen_witness())
        source['winning_literal']['failure_log2_screen'] = 0.
        # Put every unmultiplied cell at zero, isolating the outside factor.
        from math import log2
        source_count = gate.matching_joint_counts(8192, 17)[8, 0]
        prefix = (log2(source_count)+log2(comb(17, 2))+15*log2(28)-log2(255)
                  -gate.LIMIT*log2(float(gate.Z))-float(gate.LOGP)-480*log2(255))
        with patch.object(gate, 'frozen_witness', return_value=source), \
             patch.object(gate, 'affine_class_moments', return_value=((1,), 1)), \
             patch.object(gate, 'class_log_weights', return_value=np.zeros(2)), \
             patch.object(gate, 'pair_regions_screen', return_value=None) as raw, \
             patch.object(gate, 'ordered_region_tables', return_value=None), \
             patch.object(gate, 'occupancy_vectors', return_value=np.zeros((1, 3, 2))) as vectors, \
             patch.object(gate, 'branch_caps', return_value=np.zeros((1, 3))), \
             patch.object(gate, 'column_log_coefficients', return_value=np.zeros(1)), \
             patch.object(gate, 'positive_bivariate_tail', return_value=(-prefix, {})) as tail:
            result = gate.run()
        self.assertEqual(raw.call_count, 1)
        self.assertEqual(vectors.call_count, 1)
        self.assertEqual(tail.call_count, 21)
        self.assertEqual(result['ordered_fixed_dense_rowset_assignments'], 219**2)
        self.assertAlmostEqual(result['failure_log2_screen'], log2(36), places=12)
        for row in result['entries']:
            self.assertEqual(row['failure_log2_screen'], log2(row['dense_type_multiplicity']))
        self.assertTrue(result['complete_selected_cohort_and_source'])
        self.assertFalse(result['complete_source_coverage'])
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['source_colorings'], comb(17, 2))
        self.assertEqual(result['sparse_pair_assignments'], 28**15)
        self.assertEqual((result['sparse_support_floor'], result['dense_nonempty_column_floor'],
                          result['sparse_normalizer_exponent']), (510, 34, 480))

    def test_incomplete_timeout_is_not_coverage(self):
        with patch.object(gate, 'affine_class_moments', side_effect=TimeoutError('mocked')):
            result = gate.run()
        self.assertEqual(result['status'], 'compute-budget-exceeded')
        self.assertFalse(result['complete_dense_rowset_coverage'])
        self.assertFalse(result['complete_selected_cohort_and_source'])
        self.assertNotIn('failure_log2_screen', result)
        self.assertEqual(result['evaluated_arity_pairs'], 0)
        with self.assertRaises(ValueError):
            gate.run(max_seconds=121)


if __name__ == '__main__':
    unittest.main()

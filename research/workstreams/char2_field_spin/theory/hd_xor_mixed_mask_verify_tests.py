"""Receipt scope hardening and independent rounded-cell union checks."""
from copy import deepcopy
from fractions import Fraction as F
import json
import unittest

from hd_xor_mixed_mask_verify import DIRECTORY, checked_rows, rounded_union_512, verify


class MixedMaskVerifyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt = json.loads((DIRECTORY/'hd_xor_mixed_mask_8_0_mu16_receipt256.json').read_text())

    def test_all219_dual_ceilings_and_new512_sum(self):
        result = verify()
        self.assertEqual(result['fixed_row_sets'], 219)
        self.assertTrue(result['exact_partition_checked'])
        self.assertTrue(result['all_per_class_ceilings_agree'])
        self.assertTrue(result['aggregate_ceilings_agree'])
        self.assertEqual(result['independent_sum_precision_bits'], 512)
        self.assertTrue(result['scoped_class_complete'])
        self.assertFalse(result['complete_mixed_class'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertLess(F(result['rounded_cell_union_log2_upper']), -100)
        self.assertEqual(sum(row['fixed_row_sets'] for row in result['arity_summaries']), 219)

    def test_exact_power_of_two_sum(self):
        result = rounded_union_512({mask: F(-50) for mask in (7, 11, 13, 14)})
        self.assertGreaterEqual(result, -48)
        self.assertLessEqual(result, -48+F(1, 10**12))

    def test_duplicate_wrong_arity_wrong_literal_and_wrong_scope_rejected(self):
        alterations = (
            lambda data: data['entries'].__setitem__(0, data['entries'][1]),
            lambda data: data['entries'][0].__setitem__('active_constituents', 8),
            lambda data: data['entries'][0].__setitem__('mark', '99/100'),
            lambda data: data.__setitem__('cohort_counts', [1, 15, 1]),
            lambda data: data.__setitem__('mu_range', [1, 16]),
            lambda data: data.__setitem__('two_row_geometries', 28**16),
            lambda data: data.__setitem__('projective_normalizations', True),
            lambda data: data.__setitem__('source_supports', data['source_supports']-1),
        )
        for alter in alterations:
            candidate = deepcopy(self.receipt)
            alter(candidate)
            with self.assertRaises(ValueError):
                checked_rows(candidate, 256)

    def test_partial_receipt_cannot_be_promoted_to_complete(self):
        partial = json.loads((DIRECTORY/'hd_xor_mixed_mask_A195_receipt256.json').read_text())
        with self.assertRaises(ValueError):
            checked_rows(partial, 256)


if __name__ == '__main__':
    unittest.main()

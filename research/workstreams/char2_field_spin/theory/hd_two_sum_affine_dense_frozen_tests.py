"""Deterministic arity selection is outside setup randomness, not a new union."""
from copy import deepcopy
import unittest

from hd_two_sum_affine_dense import ARITY_CLASSES, CONSTRUCTION
from hd_two_sum_affine_dense_frozen import combine


def bank(label, values):
    return dict(construction=CONSTRUCTION, cohort_counts=[0, 16, 1],
        source_matching_counts=[8, 0], labels=list(range(8, 16)), field_order=256,
        forbidden_weight=734003, total_columns=32, sparse_pair_assignments=28**16,
        source_colorings=17, sparse_support_floor=544, dense_nonempty_column_floor=17,
        sparse_normalizer_exponent=512, old_xor_mask_counts_used=False, directed=False,
        witness_label=label, mark='99/100', log2_potential='-60',
        entries=[dict(active_constituents=a, active_row_masks=list(masks),
                      failure_log2_screen=value) for (a, masks), value in zip(ARITY_CLASSES, values)])


class AffineDenseFrozenTests(unittest.TestCase):
    def test_deterministic_class_choice_no_extra_bank_factor(self):
        first = bank('first', [0, 1, 2, 3, 4, 5])
        second = bank('second', [5, 4, 3, 2, 1, 0])
        result = combine(first, second)
        self.assertEqual([row['failure_log2_screen'] for row in result['selected']], [0, 1, 2, 2, 1, 0])
        from math import log2
        self.assertAlmostEqual(result['failure_log2_screen'], log2(14))
        self.assertFalse(result['setup_dependent_selection'])
        self.assertFalse(result['extra_witness_bank_factor'])
        self.assertEqual(sum(len(row['active_row_masks']) for row in result['selected']), 219)

    def test_wrong_kernel_and_changed_partition_fail_closed(self):
        initial = bank('first', [-10]*6)
        for key, value in (('labels', list(range(8))), ('old_xor_mask_counts_used', True),
                           ('source_matching_counts', [0, 8]), ('sparse_pair_assignments', 28)):
            changed = deepcopy(initial)
            changed[key] = value
            with self.assertRaises(ValueError):
                combine(initial, changed)
        changed = deepcopy(initial)
        changed['entries'][0]['active_row_masks'].pop()
        with self.assertRaises(ValueError):
            combine(initial, changed)


if __name__ == '__main__':
    unittest.main()

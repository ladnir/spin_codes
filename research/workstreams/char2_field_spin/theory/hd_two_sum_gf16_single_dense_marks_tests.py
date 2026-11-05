from fractions import Fraction as F
from math import log2
import unittest
from unittest.mock import patch

from hd_two_sum_gf16_single_dense_marks import MARKS, run
from hd_two_sum_gf16_single_dense_matrix import MARKS as OLD_MARKS


class SingleDenseMarksTests(unittest.TestCase):
    def test_only_the_four_authorized_rational_marks(self):
        self.assertEqual(MARKS, (F('1999/2000'), F('.9996'), F('.9997'), F('.9998')))
        self.assertEqual(tuple(z for _, z in OLD_MARKS), (F('499451537/500000000'), F('1199/1200')))

    def test_all14_arities_select_before_union_without_bank_factor(self):
        ns = 'hd_two_sum_gf16_single_dense_marks.'
        with patch(ns+'pair_regions_screen', return_value={}) as banks, \
                patch(ns+'single_source_tables', return_value=None), \
                patch(ns+'weight_matrices', return_value=None), \
                patch(ns+'matrix_branches', return_value=None), \
                patch(ns+'matrix_tail', return_value=(0., [dict(nonempty_dense_columns=9, log2_mass=0.)])):
            result = run()
        self.assertEqual(banks.call_count, 4)
        self.assertEqual([call.args for call in banks.call_args_list], [(float(z), 1) for z in MARKS])
        expected = log2(8192)-log2(255)-734003*log2(float(MARKS[-1]))+log2(14)
        self.assertAlmostEqual(result['failure_log2_screen'], expected, places=11)
        self.assertEqual(len(result['entries']), 14)
        self.assertEqual(sum(row['fixed_row_sets'] for row in result['entries']), 65399)
        self.assertEqual(result['endpoint_potential_factors'], 0)
        self.assertFalse(result['optimizer_used'])
        self.assertFalse(result['directed'])

    def test_budget_is_explicit_and_cannot_expand(self):
        for invalid in (True, 0, 21, float('nan')):
            with self.assertRaises(ValueError):
                run(max_seconds=invalid)
        with patch('hd_two_sum_gf16_single_dense_marks.pair_regions_screen') as bank:
            result = run(max_seconds=1e-12)
        bank.assert_not_called()
        self.assertFalse(result['complete_selected_cohort'])
        self.assertNotIn('failure_log2_screen', result)


if __name__ == '__main__':
    unittest.main()

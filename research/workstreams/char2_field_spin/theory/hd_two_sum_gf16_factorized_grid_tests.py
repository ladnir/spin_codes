"""Mocked fixed-grid scope and once-only regional-bank checks."""
import unittest
from unittest.mock import patch

import hd_two_sum_gf16_factorized_grid as grid


class FactorizedGridTests(unittest.TestCase):
    def test_nine_fixed_pairs_three_banks_no_setup_selection(self):
        def evaluate(regions, weights, mark, logp, check):
            check()
            value = float(mark)+float(logp)
            return dict(mark=str(mark), log2_potential=str(logp), failure_log2_screen=value,
                envelope_contacts=[0, 32], entries=[dict(active_groups=r, failure_log2_screen=value+r)
                                                  for r in range(18, 33)])
        with patch.object(grid, 'forced_nonzero_regions', return_value={}) as raw, \
             patch.object(grid, 'dense_arity_logs', return_value=({}, {})), \
             patch.object(grid, 'evaluate', side_effect=evaluate) as caps:
            result = grid.run()
        self.assertEqual(raw.call_count, 3)
        self.assertEqual(caps.call_count, 9)
        self.assertEqual({(row['mark'], row['log2_potential']) for row in result['witnesses']},
                         {(str(z), str(p)) for z in grid.MARKS for p in grid.LOG_POTENTIALS})
        self.assertEqual(result['completed_witnesses'], 9)
        self.assertTrue(result['complete_requested_band'])
        self.assertFalse(result['directed'])
        self.assertEqual(result['z_p_fits'], 0)
        self.assertEqual(result['witness_bank_penalty'], 1)

    def test_timeout_is_not_complete_and_first_failure_is_pinned(self):
        with patch.object(grid, 'dense_arity_logs', side_effect=TimeoutError('mocked')):
            result = grid.run()
        self.assertFalse(result['complete_requested_band'])
        self.assertEqual(result['completed_witnesses'], 0)
        self.assertNotIn('failure_log2_screen', result)
        with patch.object(grid, 'FIRST_SHA256', '0'*64), self.assertRaises(ValueError):
            grid.source_scope()


if __name__ == '__main__':
    unittest.main()

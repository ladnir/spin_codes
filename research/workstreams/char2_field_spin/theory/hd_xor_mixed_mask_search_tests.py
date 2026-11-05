"""Budget, literal, and scope checks without repeating the proof screen."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from hd_xor_mixed_mask_search import BOUNDS, bounded_search, evaluate_at


class MixedMaskSearchTests(unittest.TestCase):
    def test_search_is_bounded_to_one_mask_and_one_all_A_followup(self):
        calls = []

        def fake(mark, potential, active_masks=(255,), progress=None):
            masks = tuple(active_masks)
            z, p = float(F(mark)), float(F(potential))
            calls.append((F(mark), F(potential), masks))
            value = 100*((z-.9935)/.002)**2+((p+50)/5)**2
            return dict(failure_log2_screen=value, elapsed_seconds=0.,
                        entries=[dict(largest_dense_degree=32, largest_sparse_degree_at_peak=1000)])

        with patch('hd_xor_mixed_mask_search.evaluate_at', side_effect=fake):
            result = bounded_search()
        trials = result['proposals']
        self.assertLessEqual(len(trials), 36)
        self.assertGreater(len(trials), 1)
        self.assertEqual(len(calls), len(trials)+1)
        self.assertTrue(all(masks == (255,) for _, _, masks in calls[:-1]))
        self.assertEqual(calls[-1][2], tuple(mask for mask in range(256) if mask.bit_count() >= 3))
        for trial, (z, p, _) in zip(trials, calls):
            self.assertEqual(F(trial['mark']), z)
            self.assertEqual(F(trial['log2_potential']), p)
            for value, (low, high) in zip((float(z), float(p)), BOUNDS):
                self.assertLessEqual(low, value)
                self.assertLessEqual(value, high)
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['complete_mixed_class'])
        self.assertTrue(result['exact_positive_bivariate_tail'])
        self.assertEqual(result['sparse_support_floor'], 544)
        self.assertEqual(result['dense_nonempty_column_floor'], 17)

    def test_invalid_witness_or_dense_scope_is_rejected_before_work(self):
        for mark, potential, masks in ((1, -1, (255,)), (F(1, 2), 1, (255,)),
                                      (F(1, 2), -1, (255, 255)), (F(1, 2), -1, (3,))):
            with self.assertRaises(ValueError):
                evaluate_at(mark, potential, masks)


if __name__ == '__main__':
    unittest.main()

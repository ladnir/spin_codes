"""Bounded fit bookkeeping, with mocked objectives and no operator run."""
from fractions import Fraction as F
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from hd_two_sum_affine_dense_two_fit import fit, Z, LOGP


class AffineDenseTwoFitTests(unittest.TestCase):
    def test_budget_exact_literals_and_one_fresh_replay(self):
        def objective(z, p, mass):
            return dict(mark=str(F(z)), log2_potential=str(F(p)),
                failure_log2_screen=float(z)*1000+float(p), elapsed_seconds=0.,
                largest_dense_degree=34, largest_sparse_degree=510)
        def optimizer(fn, start, **kwargs):
            self.assertEqual(kwargs['options']['maxfev'], 23)
            for point in ((0., 0.), (-3., 0.), (-3., -4.)):
                fn(point)
            return SimpleNamespace(success=False, message='mocked bounded calls')
        namespace = 'hd_two_sum_affine_dense_two_fit.'
        with patch(namespace+'evaluate', side_effect=objective) as evaluated, \
                patch(namespace+'minimize', side_effect=optimizer):
            result = fit()
        self.assertEqual(evaluated.call_count, 4)
        self.assertEqual(result['optimizer_evaluations'], 3)
        self.assertEqual(result['fresh_literal_replays'], 1)
        self.assertEqual(result['maximum_evaluations'], 24)
        self.assertEqual(result['total_evaluations'], 4)
        self.assertEqual(result['dense_arities'], [8, 8])
        self.assertFalse(result['directed'])
        self.assertFalse(result['complete_dense_rowset_coverage'])
        winner = result['winning_literal']
        self.assertIsInstance(F(winner['mark']), F)
        self.assertEqual(evaluated.call_args.args[:2], (F(winner['mark']), F(winner['log2_potential'])))


if __name__ == '__main__':
    unittest.main()

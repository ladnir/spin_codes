from fractions import Fraction as F
import unittest
from unittest.mock import patch

from hd_shifted_swap_search import exact_witness, proposal
from hd_shifted_swap_screen import CONSTRUCTION


class SwapSearchTests(unittest.TestCase):
    def test_finite_objective_budget_and_explicit_float_scope(self):
        calls = []
        def fake_bound(profile, z, alpha, odds, **options):
            calls.append(options['method'])
            return dict(profile=list(profile), mark=z, alpha=alpha,
                        marker_log_odds=[list(row) for row in odds],
                        log2_upper_screen=(z-.5)**2+(alpha-.1)**2,
                        directed=False, whole_code_certificate=False)
        with patch('hd_shifted_swap_search.bound', side_effect=fake_bound):
            result = proposal((1024, 400, 400), max_evaluations=25, grid=64)
        self.assertLessEqual(result['objective_evaluations'], 25)
        self.assertLessEqual(calls.count('global'), 25)
        self.assertEqual(calls.count('alias'), 1)
        self.assertEqual(result['construction'], CONSTRUCTION)
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])

    def test_rational_witness_and_work_validation(self):
        result = exact_witness(dict(mark=.5, alpha=.25,
                                    marker_log_odds=[[-3., -4.], [-5., -6.]]))
        self.assertEqual(F(result['z']), F(1, 2))
        self.assertEqual(F(result['alpha']), F(1, 4))
        self.assertEqual(F(result['markers'][0][0]), F(1, 8))
        self.assertTrue(result['pair_swaps'])
        with self.assertRaises(ValueError):
            proposal((1, 0, 0), max_evaluations=0)


if __name__ == '__main__':
    unittest.main()

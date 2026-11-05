from copy import deepcopy
from fractions import Fraction as F
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from hd_shifted_potential_retune import retune_potential


class RetuneTests(unittest.TestCase):
    def test_preserves_bank_and_source_but_exactifies_one_common_potential(self):
        original = dict(z='1/2', alpha='1/10', markers=['1 2'], potential='1/2',
                        alpha_and_potential='unchanged-from-source')
        screen = SimpleNamespace(source_witness=original, log_ratio=-1.,
                                 regional=lambda weight: ('same-region', weight))
        observed = []
        def objective(pairs, weight, region, ratio):
            observed.append((pairs, weight, region))
            return (ratio+3.)**2
        with patch('hd_shifted_potential_retune.fixed_screen_from_artifact', return_value=screen), \
                patch('hd_shifted_potential_retune.run_log_screen_float', side_effect=objective):
            result, stats = retune_potential({}, 640)
        self.assertEqual(original['potential'], '1/2')
        self.assertAlmostEqual(result.log_ratio, -3., places=6)
        self.assertEqual(F(result.source_witness['potential']), F(2.**result.log_ratio))
        self.assertNotIn('alpha_and_potential', result.source_witness)
        self.assertTrue(result.source_witness['alpha_preserved_from_source'])
        for key in ('z', 'alpha', 'markers'):
            self.assertEqual(result.source_witness[key], original[key])
        self.assertTrue(all(row == (4096, 640, ('same-region', 640)) for row in observed))
        self.assertLessEqual(stats['selected_support_sum'], stats['initial_support_sum'])
        self.assertFalse(stats['directed'])
        self.assertEqual(stats['new_marker_fits'], 0)

    def test_original_candidate_is_retained_when_optimizer_is_worse(self):
        screen = SimpleNamespace(source_witness={'potential': '1/8'}, log_ratio=-3.,
                                 regional=lambda weight: {})
        with patch('hd_shifted_potential_retune.fixed_screen_from_artifact', return_value=screen), \
                patch('hd_shifted_potential_retune.run_log_screen_float',
                      side_effect=lambda p, w, r, v: (v+3)**2), \
                patch('hd_shifted_potential_retune.minimize_scalar',
                      return_value=SimpleNamespace(x=-2., nfev=1)):
            result, _ = retune_potential({}, 640)
        self.assertEqual(result.source_witness['potential'], '1/8')

    def test_invalid_inputs_fail_before_build(self):
        with patch('hd_shifted_potential_retune.fixed_screen_from_artifact') as build:
            for weight in (True, 0, 8193, 1.5):
                with self.assertRaises(ValueError):
                    retune_potential({}, weight)
            for bounds in ((0., 0.), (-513., 0.), (0., float('inf'))):
                with self.assertRaises(ValueError):
                    retune_potential({}, 640, bounds=bounds)
            build.assert_not_called()


if __name__ == '__main__':
    unittest.main()

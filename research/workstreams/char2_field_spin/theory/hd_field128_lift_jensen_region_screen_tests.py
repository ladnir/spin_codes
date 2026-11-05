"""Tiny Jensen, ordered-state, normalization, and scope tests; no gate."""

from fractions import Fraction as F
from itertools import product
import unittest

import numpy as np

import hd_field128_lift_jensen_region_screen as jensen


old = jensen.old


class JensenRegionTests(unittest.TestCase):
    def test_conditional_jensen_and_independent_persistent_regions(self):
        # Fixed noncommuting raw operators; two routes with nonuniform exact probabilities.
        a = np.array([[.5, .2], [.1, .7]])
        b = np.array([[.4, .3], [.2, .8]])
        routes = (a @ b, b @ a)
        probabilities = (F(1, 3), F(2, 3))
        alpha = .5
        target = sum(float(p) * matrix**alpha for p, matrix in zip(probabilities, routes))
        mean = sum(float(p) * matrix for p, matrix in zip(probabilities, routes))
        cap = mean**alpha
        self.assertTrue(np.all(target <= cap + 1e-15))
        self.assertTrue(np.any(target < cap - 1e-8))
        self.assertFalse(np.allclose(routes[0], routes[1]))
        # Expand the actual scalar path sum first. Regional state is not reset.
        actual = sum(float(probabilities[i] * probabilities[k])
                     * np.sum((routes[i] @ routes[k])[0])**alpha for i, k in product(range(2), repeat=2))
        self.assertLessEqual(actual, np.sum((target @ target)[0]) + 1e-14)
        self.assertLessEqual(actual, np.sum((cap @ cap)[0]) + 1e-14)
        # Both regional upper bounds dominate the same target, permitting entrywise min.
        alternative = target + np.array([[.03, .001], [.001, .03]])
        mixed = np.minimum(cap, alternative)
        self.assertTrue(np.all(target <= mixed + 1e-15))
        self.assertLessEqual(actual, np.sum((mixed @ mixed)[0]) + 1e-14)

    def test_normalizer_is_inside_alpha(self):
        alpha = .25
        mean = old.regional_screen(.5, 1., 2, q=256, paired=False, max_seconds=10)
        actual = jensen.regional_screen(.5, alpha, 2, q=256, mode='jensen', max_seconds=10)
        for key in mean:
            np.testing.assert_allclose(actual[key], alpha * np.asarray(mean[key]), rtol=0, atol=1e-13)
        # Applying alpha to an unnormalized numerator and dividing afterward is too small.
        normalizer_log = old.log_multinomial((old.PAIRS - 1, 1, 0))
        normalized = np.asarray(mean[1, 0])
        wrong = alpha * (normalized + normalizer_log) - normalizer_log
        self.assertTrue(np.all(wrong < actual[1, 0]))

    def test_default_min8_and_alpha1_equality(self):
        for alpha in (.25, 1.):
            selected, pure = jensen.regional_tables(.5, alpha, 2, q=256, max_seconds=10)
            baseline = jensen.eight.regional_screen(.5, alpha, 2, q=256, max_seconds=10)
            for key in selected:
                np.testing.assert_allclose(selected[key], np.minimum(pure[key], baseline[key]), rtol=0, atol=1e-12)
                if alpha == 1:
                    np.testing.assert_allclose(pure[key], baseline[key], rtol=0, atol=1e-12)

    def test_scope_metadata_on_tiny_selected_weight(self):
        result = jensen.weight_screen(1, .5, .25, q=256, profile=(0, 0), max_seconds=10)
        self.assertEqual(result['combination_mode'], 'min8')
        self.assertLessEqual(result['log2_upper_screen'], result['pure_jensen_log2_upper_screen'] + 1e-10)
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['covers_all_supports_at_selected_weight'])
        self.assertTrue(result['state_persists_across_regions'])
        self.assertTrue(result['normalizer_inside_fractional_power'])
        pure = jensen.weight_screen(1, .5, .25, q=256, mode='jensen', max_seconds=10)
        self.assertEqual(pure['log2_upper_screen'], pure['pure_jensen_log2_upper_screen'])
        self.assertTrue(pure['covers_all_supports_at_selected_weight'])

    def test_guards_and_import_closure(self):
        for kwargs in (dict(mode='unknown'), dict(max_v_full=3)):
            with self.assertRaises(ValueError):
                jensen.regional_screen(.5, .1, 2, **kwargs)
        with self.assertRaises(ValueError):
            jensen.regional_screen(.5, .1, 65)
        with self.assertRaises(TimeoutError):
            jensen.regional_screen(.5, .1, 2, max_seconds=1e-12)
        pins = jensen.source_pins()
        for name in ('hd_field128_lift_jensen_region_screen.py', 'hd_field128_lift_eight_packet_screen.py',
                     'hd_field128_lift_four_packet_screen.py', 'hd_field128_lift_screen.py'):
            self.assertIn(name, pins)


if __name__ == '__main__':
    unittest.main()

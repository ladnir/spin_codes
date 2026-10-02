import copy
from contextlib import ExitStack, redirect_stdout
from fractions import Fraction as Q
import io
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from flint import arb, arb_mat, ctx
import region_cache as cache_module

rc = cache_module.rc


class TinyModel:
    def __init__(self):
        self.data = {'windows': 32, 'authenticated_marker': object()}
        self.tilt = Q(3, 16)
        self.features = [Q(1, 4), Q(3, 4)]
        self.active = [1, 1]
        self.coefficients = [Q(1, 2), Q(1, 3)]
        self.q_min = 1
        self.threshold = 2
        self.variance_shuffle = True
        self.variance_bins = 1
        self.regional_count = True
        self.loss_cache = {}

    def weights(self, cell, tilt):
        return [1-cell[0], cell[1]/tilt], None

    def family(self, tilt):
        return self.coefficients, None, None

    def comparison_loss(self, cell, tilt, dual=None):
        return Q(1)


def witness(direct=False):
    part = dict(interval=['0', '1/4'], dual=['0', '0', '0'])
    return dict(tilt='3/16', parameters=['1/100', '0', '0'], variance_dual=['0', '0'],
        variance_partition=[part], regional_direct_counts=direct,
        regional_count_parts=[dict(part, mgf_witnesses=[dict(tilt='0', dual=['0', '0', '0'])])])


class RegionCacheTests(unittest.TestCase):
    def setUp(self):
        ctx.prec = 192
        self.model = TinyModel()
        self.cache = cache_module.RegionalCache()
        self.polynomial = [arb_mat([[1, arb(1)/8], [0, arb(1)/2]]),
                           arb_mat([[arb(1)/4, arb(1)/2], [arb(1)/8, arb(1)/4]]),
                           arb_mat([[arb(1)/8, arb(1)/4], [arb(1)/2, arb(1)/8]])]
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.stack.enter_context(patch.object(rc.sc, 'G', 2))
        self.stack.enter_context(patch.object(rc.sc, 'REGIONS', 3))
        self.stack.enter_context(patch.object(rc.sc, 'PACKETS', 6))
        self.real_placement = rc.placement
        self.local = self.stack.enter_context(patch.object(rc, 'local_operators', return_value=object()))
        self.place = self.stack.enter_context(patch.object(rc, 'placement', side_effect=lambda _: [arb_mat(m) for m in self.polynomial]))

    def test_exact_outward_equivalence_and_reuse(self):
        for direct in (False, True):
            candidate = witness(direct)
            for cell in ((Q(1, 4), Q(1, 2)), (Q(1, 3), Q(5, 12))):
                with redirect_stdout(io.StringIO()):
                    expected, checked = rc.outward(self.model, cell, candidate)
                actual, actual_checked = self.cache.outward(self.model, cell, candidate)
                self.assertEqual(actual, expected)
                self.assertEqual(actual_checked, checked)
        self.assertEqual(self.cache.misses, 1)
        self.assertEqual(self.cache.hits, 3)
        self.assertEqual(self.place.call_count, 5)  # Four uncached controls, one cached build.

    def test_mixture_changes_recompute_outer_bound(self):
        cell, candidate = (Q(1, 4), Q(1, 2)), witness(True)
        first, _ = self.cache.outward(self.model, cell, candidate)
        self.model.coefficients = [Q(1, 8), Q(1, 9)]
        actual, _ = self.cache.outward(self.model, cell, candidate)
        with redirect_stdout(io.StringIO()):
            expected, _ = rc.outward(self.model, cell, candidate)
        self.assertEqual(actual, expected)
        self.assertNotEqual(actual, first)
        self.assertEqual(self.cache.misses, 1)

    def test_precision_lambda_data_and_each_selector_invalidate(self):
        candidate = witness()
        self.cache.polynomial(self.model, candidate)
        self.assertTrue(self.cache.is_warm(self.model, candidate))
        cases = [dict(regional_exact_zero=True), dict(regional_joint_return_through=1),
            dict(regional_lazy_density_through=1),
            dict(regional_feedback_classes_from=1, regional_feedback_classes_through=2),
            dict(regional_feedback_classes_from=2, regional_feedback_classes_through=2),
            dict(regional_feedback_classes_from=2, regional_feedback_classes_through=3),
            dict(regional_feedback_classes_from=2, regional_feedback_classes_through=3,
                 regional_feedback_uniform_classes=True),
            dict(regional_feedback_classes_from=2, regional_feedback_classes_through=3,
                 regional_feedback_uniform_classes=True, regional_feedback_uniform_replace=True)]
        for changes in cases:
            changed = dict(candidate, **changes)
            self.assertFalse(self.cache.is_warm(self.model, changed))
            self.cache.polynomial(self.model, changed)
            self.assertTrue(self.cache.is_warm(self.model, changed))
        changed = dict(candidate, parameters=['1/99', '0', '0'])
        self.assertFalse(self.cache.is_warm(self.model, changed))
        self.cache.polynomial(self.model, changed)
        ctx.prec = 256
        self.assertFalse(self.cache.is_warm(self.model, candidate))
        self.cache.polynomial(self.model, candidate)
        self.model.data = dict(self.model.data)
        self.assertFalse(self.cache.is_warm(self.model, candidate))
        self.cache.polynomial(self.model, candidate)
        self.assertEqual(self.cache.misses, 12)
        self.assertEqual(self.cache.hits, 0)

    def test_warmth_check_is_read_only_and_clear_makes_it_cold(self):
        candidate = witness()
        original = copy.deepcopy(candidate)
        self.assertFalse(self.cache.is_warm(self.model, candidate))
        self.assertEqual(self.place.call_count, 0)
        self.assertEqual((self.cache.hits, self.cache.misses), (0, 0))
        self.cache.polynomial(self.model, candidate)
        saved = self.cache._saved
        for _ in range(3):
            self.assertTrue(self.cache.is_warm(self.model, candidate))
        # Outer duals and count selectors do not alter the local operator.
        changed = dict(candidate, parameters=['1/100', '2', '3'], regional_direct_counts=True)
        self.assertTrue(self.cache.is_warm(self.model, changed))
        self.assertIs(self.cache._saved, saved)
        self.assertEqual(candidate, original)
        self.assertEqual(self.place.call_count, 1)
        self.assertEqual((self.cache.hits, self.cache.misses), (0, 1))
        self.cache.clear()
        self.assertFalse(self.cache.is_warm(self.model, candidate))

    def test_public_results_do_not_mutate_retained_polynomial(self):
        first = self.cache.polynomial(self.model, witness())
        first[0][0, 0] = 999
        second = self.cache.polynomial(self.model, witness())
        self.assertEqual(second[0][0, 0], 1)
        self.assertEqual(self.cache.misses, 1)
        self.cache.clear()
        self.cache.polynomial(self.model, witness())
        self.assertEqual(self.cache.misses, 2)

    def test_unknown_and_malformed_local_selectors_rejected(self):
        for changes in (dict(regional_new_refinement=True), dict(regional_exact_zero=1),
                        dict(regional_joint_return_through=True), dict(regional_lazy_density_through=33),
                        dict(regional_feedback_classes_from=1), dict(regional_feedback_uniform_replace=True)):
            self.assertRaises(ValueError, self.cache.polynomial, self.model, dict(witness(), **changes))
            self.assertRaises(ValueError, self.cache.is_warm, self.model, dict(witness(), **changes))
        self.assertEqual(self.place.call_count, 0)

    def test_warmth_checks_function_identity_and_malformed_lambda(self):
        candidate = witness()
        self.cache.polynomial(self.model, candidate)
        with patch.object(rc, 'placement'):
            self.assertFalse(self.cache.is_warm(self.model, candidate))
        with patch.object(rc, 'local_operators'):
            self.assertFalse(self.cache.is_warm(self.model, candidate))
        self.assertTrue(self.cache.is_warm(self.model, candidate))
        for params in (['0', '0', '0'], ['1'], ['-1', '0', '0']):
            self.assertRaises(ValueError, self.cache.is_warm, self.model,
                              dict(candidate, parameters=params))
        self.assertEqual((self.cache.hits, self.cache.misses), (0, 1))

    def test_construction_precision_change_and_bad_geometry_rejected(self):
        def changed(_):
            ctx.prec += 1
            return self.polynomial
        self.place.side_effect = changed
        self.assertRaises(ArithmeticError, self.cache.polynomial, self.model, witness())
        self.place.side_effect = lambda _: self.polynomial[:2]
        self.assertRaises(ArithmeticError, self.cache.polynomial, self.model, witness())

    def test_bound_adapter_dispatch_and_validation(self):
        cell = (Q(1, 4), Q(1, 2))
        candidate = witness()
        actual = self.cache.bound(self.model, cell, candidate)
        with redirect_stdout(io.StringIO()):
            expected, _ = rc.outward(self.model, cell, candidate)
        self.assertEqual(actual, expected)
        nonregional = dict(candidate)
        nonregional.pop('regional_count_parts')
        with patch.object(self.model, 'outward', create=True, return_value=arb(7)) as base:
            self.assertEqual(self.cache.bound(self.model, cell, nonregional), 7)
            base.assert_called_once_with(cell, nonregional)
        for changed in (dict(candidate, parameters=['1/100', '0', '-1']),
                        dict(candidate, parameters=['0', '0', '0']),
                        dict(candidate, parameters=['1']),
                        dict(candidate, tilt='1/4'),
                        dict(candidate, regional_new_refinement=True)):
            self.assertRaises(ValueError, self.cache.bound, self.model, cell, changed)
        changed = dict(candidate)
        changed.pop('variance_dual')
        self.assertRaises(ValueError, self.cache.bound, self.model, cell, changed)
        self.model.regional_count = False
        self.assertRaises(ValueError, self.cache.bound, self.model, cell, candidate)

    def test_actual_small_ordered_placement(self):
        # Exercise the real noncommuting polynomial-placement implementation
        # on two one-packet epochs; no production census or proof search runs.
        self.local.return_value = self.polynomial[:2]
        self.place.side_effect = lambda local: self.real_placement(local, epochs=2)
        candidate, cell = witness(True), (Q(1, 4), Q(1, 2))
        with redirect_stdout(io.StringIO()):
            expected, _ = rc.outward(self.model, cell, candidate)
        self.assertEqual(self.cache.bound(self.model, cell, candidate), expected)
        self.assertEqual(self.cache.bound(self.model, cell, candidate), expected)
        self.assertEqual(self.cache.misses, 1)
        self.assertEqual(self.cache.hits, 1)


if __name__ == '__main__':
    unittest.main()

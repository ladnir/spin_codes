"""Mock-only contract tests for the original-weight dense mark-search adapter."""

from contextlib import ExitStack
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import hd_field128_gamma_midpoint_dense_mark_search as dense


class DenseSearchTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name)/'proposal.json'
        self.weight = 5120
        self.source = dict(
            construction=dense.provider.gamma.construction(30), state_dimension=30,
            target_distance='2/5', requested_weights=[self.weight],
            witness=dict(z=str(F(0.1)), alpha='1/8',
                         markers=[['2/3', '5/7'], ['11/13', '17/19']], grid=128))
        self.save()
        self.calls = []
        self.checkpoints = []

    def save(self):
        self.path.write_text(json.dumps(self.source), encoding='utf-8')

    def mocks(self, *, points=1, fail_at=None, failure=TimeoutError,
              wrong_shape=False, wrong_profiles=False, count_failure=None):
        owner = self
        class Counts:
            def __init__(self, pairs, weight, **options):
                owner.calls.append(('counts', pairs, weight, options['max_singletons']))
                if count_failure:
                    raise count_failure('mock count resource stop')
            def diagnostics(self):
                return dict(mock=True)
        class Coefficients:
            def __init__(self, z, alpha, **options):
                previous = sum(call[0] == 'coeff' for call in owner.calls)
                owner.calls.append(('coeff', F(z), F(alpha), options['state']))
                if previous == fail_at:
                    raise failure('mock bank resource stop')
                self.composed = np.array([z, alpha])
        class Bank:
            def __init__(self, z, alpha, odds, potential, **options):
                owner.calls.append(('bank', tuple(odds), potential, options['floor_policy']))
                self.alpha, self.density, self.constant = alpha, 3., -7.
                self.regional_diagnostics = dict(mock=True, selected_positive_fallback_entries=0)
            def regional(self, weight):
                owner.calls.append(('regional', weight))
                count = len(dense.checked_domain(weight))
                return np.zeros((count-int(wrong_shape), 2, 2))
            def metadata(self):
                return dict(construction=owner.source['construction'])
        def support(pairs, weight, regions, **options):
            owner.calls.append(('support', pairs, weight, regions.shape, options['counts']))
            return dict(log2_support_sum_screen=-123., largest_profile=[weight-4096, weight//2],
                        profile_pairs_evaluated=len(regions)**2-int(wrong_profiles))
        def minimize(objective, initial, **options):
            self.assertEqual(options['method'], 'Nelder-Mead')
            objective(initial.copy())
            if points >= 2:
                # Repeated exact marks must use the cache, not rebuild a bank.
                objective(initial.copy())
                objective(initial+np.array([.25, 0.]))
            return SimpleNamespace(success=False, message='mock evaluation limit', nfev=points)
        stack = ExitStack()
        stack.enter_context(patch.object(dense, 'source_pins', return_value={}))
        stack.enter_context(patch.object(dense.mode, 'MidpointCounts', Counts))
        stack.enter_context(patch.object(dense.mode, 'support_sum', side_effect=support))
        stack.enter_context(patch.object(dense.provider, 'AliasCoefficients', Coefficients))
        stack.enter_context(patch.object(dense.provider, 'AliasBankScreen', Bank))
        stack.enter_context(patch.object(dense, 'minimize', side_effect=minimize))
        return stack

    def run_mock(self, **options):
        with self.mocks(**options):
            return dense.run_search(self.path, self.weight, max_evaluations=3,
                                    checkpoint=lambda row: self.checkpoints.append(row) or 'mock-link')

    def test_original_weight_complete_domain_prefix_and_literals(self):
        result = self.run_mock(points=2)
        row = result['history'][0]
        full = self.checkpoints[0]
        self.assertEqual(result['completed_trials'], 2)
        self.assertEqual(result['attempted_banks'], 2)
        self.assertEqual(result['counts_built'], 1)
        self.assertEqual(self.calls[0], ('counts', 4096, 5120, 64))
        self.assertTrue(all(call[1] == 5120 for call in self.calls if call[0] == 'regional'))
        self.assertTrue(all(call[1:4] == (4096, 5120, (1537, 2, 2))
                            for call in self.calls if call[0] == 'support'))
        self.assertEqual(full['feasible_full_range'], [1024, 2560])
        self.assertEqual(full['canonical_count_weight'], 3072)
        self.assertEqual(full['count_full_shift'], 1024)
        self.assertEqual(full['profile_pairs_evaluated'], 1537**2)
        self.assertEqual(full['log2_outer_prefix'], -7.+.125*3.*5120)
        self.assertEqual(full['log2_upper_screen'], full['log2_outer_prefix']-123.)
        self.assertEqual(full['witness']['markers'], self.source['witness']['markers'])
        self.assertEqual(full['witness']['z'], self.source['witness']['z'])
        self.assertEqual(full['witness']['alpha'], '1/8')
        self.assertEqual(row['checkpoint'], 'mock-link')
        self.assertFalse(full['directed'])
        self.assertFalse(full['selected_weight_certificate'])
        self.assertFalse(full['whole_code_certificate'])
        self.assertFalse(result['common_potential_used'])
        self.assertTrue(result['count_complement_only'])
        banks = [call for call in self.calls if call[0] == 'bank']
        self.assertEqual(banks[0], banks[1])
        self.assertEqual(banks[0][2:], (0., 'positive-fallback'))

    def test_domain_and_resource_guards(self):
        for weight in (4095, 8192, True, 5120.):
            with self.assertRaises(ValueError):
                dense.checked_domain(weight)
        self.assertEqual(list(dense.checked_domain(8191)), [4095])
        self.assertEqual(int(dense.checked_domain(4096)[0]), 0)
        for kwargs in (dict(max_evaluations=0), dict(max_evaluations=41),
                       dict(max_evaluations=True), dict(max_seconds=0),
                       dict(max_seconds=121), dict(max_seconds=float('nan'))):
            with self.assertRaises(ValueError):
                dense.run_search(self.path, self.weight, **kwargs)

    def test_source_weight_and_exact_initial_mark_guards(self):
        for change in ('weight', 'mark', 'marker', 'grid'):
            with self.subTest(change=change):
                original = json.loads(json.dumps(self.source))
                if change == 'weight':
                    self.source['requested_weights'] = [3072]
                elif change == 'mark':
                    self.source['witness']['z'] = '1/3'
                elif change == 'marker':
                    self.source['witness']['markers'][0][0] = '0'
                else:
                    self.source['witness']['grid'] = 17
                self.save()
                with self.mocks(), self.assertRaises(ValueError):
                    dense.run_search(self.path, self.weight)
                self.assertFalse(self.calls)
                self.source = original
                self.save()

    def test_timeout_and_memory_preserve_only_completed_trials(self):
        for failure in (TimeoutError, MemoryError):
            with self.subTest(failure=failure.__name__):
                self.calls.clear()
                self.checkpoints.clear()
                result = self.run_mock(points=2, fail_at=1, failure=failure)
                self.assertEqual(result['completed_trials'], 1)
                self.assertEqual(result['attempted_banks'], 2)
                self.assertEqual(len(self.checkpoints), 1)
                self.assertEqual(result['optimizer']['resource_stop'], failure.__name__)
                self.assertFalse(result['optimizer']['success'])
                self.assertFalse(result['failures'][0]['complete_selected_weight_screen'])
                self.assertEqual(result['best'], self.checkpoints[0])

    def test_resource_stop_before_first_complete_trial(self):
        for options in (dict(count_failure=MemoryError), dict(fail_at=0, failure=TimeoutError)):
            with self.subTest(options=options):
                self.calls.clear()
                result = self.run_mock(**options)
                self.assertEqual(result['completed_trials'], 0)
                self.assertIsNone(result['best'])
                self.assertFalse(result['whole_code_certificate'])
                self.assertFalse(result['optimizer']['success'])

    def test_wrong_domain_or_incomplete_profiles_fail_closed(self):
        for options in (dict(wrong_shape=True), dict(wrong_profiles=True)):
            with self.subTest(options=options):
                self.calls.clear()
                result = self.run_mock(**options)
                self.assertEqual(result['completed_trials'], 0)
                self.assertIsNone(result['best'])
                self.assertEqual(result['failures'][0]['error_type'], 'ArithmeticError')
                self.assertFalse(result['failures'][0]['complete_selected_weight_screen'])
                self.assertFalse(self.checkpoints)


if __name__ == '__main__':
    unittest.main()

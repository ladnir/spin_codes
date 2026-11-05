"""Tiny polynomial/offset tests and mocked orchestration; no coefficient gate."""

from contextlib import ExitStack
from fractions import Fraction as F
import io
import json
from math import log2
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import hd_field128_gamma_midpoint_dense_calibrate as mod


class Budget:
    def check(self):
        pass


def synthetic_source(path):
    markers = [[str(F(i+2, 128)), str(F(i+2, 64))] for i in range(79)]
    data = dict(construction=mod.provider.gamma.construction(30), state_dimension=30,
                target_distance='2/5', active_groups=4096, source_sha256={},
                witness=dict(z='1/2', alpha='1/4', markers=markers, grid=128))
    path.write_text(json.dumps(data), encoding='utf-8')
    return data


class MockRuntime:
    def __init__(self, source, *, timeout_fit=None, mutate=False):
        self.source, self.timeout_fit, self.mutate = source, timeout_fit, mutate
        self.coefficient_calls = self.count_calls = self.bank_calls = self.fit_calls = self.sum_calls = 0
        self.coeff = SimpleNamespace(composed=np.array([1., 2., 3.]))

    def coefficients(self, z, alpha, **kwargs):
        self.coefficient_calls += 1
        assert kwargs['state'] == 30 and (z, alpha) == (.5, .25)
        return self.coeff

    def bank(self, z, alpha, odds, ratio, **kwargs):
        self.bank_calls += 1
        assert kwargs['_coefficients'] is self.coeff and ratio == 0
        return SimpleNamespace(constant=5., alpha=.25, density=2., score=-10.*self.bank_calls,
                               metadata=lambda: {'synthetic': True, 'markers': len(odds)})

    def regions(self, bank, weight, budget):
        n = len(mod.provider.alias.full_counts(4096, weight))
        return np.full((n, 2, 2), bank.score), np.full((n, 2, 2), True)

    def counts(self, pairs, weight, **kwargs):
        self.count_calls += 1
        assert (pairs, weight) == (4096, 5120)
        return object()

    def fit(self, numerator, pairs, weight, full, entry, odds, budget, **kwargs):
        self.fit_calls += 1
        if self.fit_calls == self.timeout_fit:
            raise TimeoutError('synthetic fit budget')
        return (F(100000+4*full+entry), F(1)), dict(entry=entry, full=full)

    def support(self, pairs, weight, regions, **kwargs):
        self.sum_calls += 1
        assert (pairs, weight) == (4096, 5120) and regions.shape == (1537, 2, 2)
        if self.mutate:
            self.source.write_text(self.source.read_text()+' ', encoding='utf-8')
        return dict(profile_pairs_evaluated=len(regions)**2, log2_support_sum_screen=float(np.max(regions)),
                    largest_profile=dict(even_full=1024, odd_full=1024))

    def patched(self):
        stack = ExitStack()
        for target, replacement in ((mod.provider, ('AliasCoefficients', self.coefficients)),
                                    (mod.provider, ('AliasBankScreen', self.bank)),
                                    (mod.mode, ('MidpointCounts', self.counts)),
                                    (mod.mode, ('support_sum', self.support)),
                                    (mod, ('regional_with_selection', self.regions)),
                                    (mod, ('fit_entry', self.fit)),
                                    (mod, ('source_pins', lambda: {}))):
            stack.enter_context(patch.object(target, *replacement))
        stack.enter_context(patch('sys.stdout', new=io.StringIO()))
        return stack


class DenseCalibrationTests(unittest.TestCase):
    def test_actual_dense_geometry(self):
        self.assertEqual(mod.target_full_counts(4096, 5120), (1024, 1025, 1028, 1040, 1408, 1792, 2176, 2560))
        self.assertEqual(mod.histogram(4096, 5120, 1024), (0, 3072, 1024))
        self.assertEqual(mod.histogram(4096, 5120, 2560), (1536, 0, 2560))
        with self.assertRaises(ValueError):
            mod.histogram(4096, 5120, 1023)

    def test_four_boundary_entries_positive_polynomial_oracle(self):
        # Two degree4096 monomials per entry; normalized optimum u-v=i+1.
        def numerator(odds):
            u, v = odds
            return np.array([np.logaddexp2(3073*u+1023*v-d, 3071*u+1025*v+d)
                             for d in (1., 2., 3., 4.)])
        proposals = []
        for entry in range(4):
            pair, record = mod.fit_entry(numerator, 4096, 5120, 1024, entry, [(0., 0.)], Budget())
            proposals.append(pair)
            self.assertLess(abs(log2(float(pair[0]))-log2(float(pair[1]))-(entry+1)), .001)
            self.assertEqual(pair[1], F(2)**64)
            self.assertEqual(tuple(map(F, record['positive_marker_literals'])), pair)
            self.assertEqual(tuple(F(2.**float(F(x))) for x in record['fitted_log_literals']), pair)
            self.assertLessEqual(record['objective_calls'], 128)
            self.assertEqual(record['entry'], entry)
        self.assertEqual(len(set(proposals)), 4)

    def test_interior_and_zero_single_finite_seeds(self):
        def numerator(odds):
            u, v = odds
            return np.full(4, 8*np.logaddexp2(0., np.logaddexp2(u, v)))
        pair, record = mod.fit_entry(numerator, 8, 10, 3, 0, [], Budget())
        self.assertAlmostEqual(log2(float(pair[0])), 2., places=2)
        self.assertAlmostEqual(log2(float(pair[1])), log2(3.), places=2)
        pair, record = mod.fit_entry(numerator, 8, 10, 5, 0, [], Budget())
        self.assertEqual(pair[0], F(1, 2**64))
        self.assertAlmostEqual(log2(float(pair[1])), log2(F(5, 3)), places=3)

    def test_regional_selection_uses_nonzero_full_offset(self):
        tables = []
        for offset in (0., -3.):
            logs = np.arange(64, dtype=float).reshape((4, 4, 4))+offset
            reliable = np.ones((4, 4, 4), dtype=bool)
            reliable[1, 2, 2] = False
            tables.append(SimpleNamespace(logs=logs, reliable=reliable, positive_numerator=np.full(4, -20.)))
        odds = [(1., 2.), (0., 1.)]
        expected, flags = [], []
        for f in (2, 3):
            s = 6-2*f
            values, chosen = [], []
            for entry in range(4):
                choices = [(float(table.logs[entry, s, f]) if table.reliable[entry, s, f]
                            else float(table.positive_numerator[entry]))-s*u-f*v
                           for (u, v), table in zip(odds, tables)]
                index = int(np.argmin(choices))
                values.append(choices[index])
                chosen.append(tables[index].reliable[entry, s, f])
            expected.append(values)
            flags.append(chosen)
        expected = np.asarray(expected).reshape((2, 2, 2))
        bank = SimpleNamespace(grid=4, odds=odds, tables=tables, regional=lambda r: expected.copy())
        actual, reliable = mod.regional_with_selection(bank, 6, Budget(), pairs=4)
        np.testing.assert_array_equal(actual, expected)
        np.testing.assert_array_equal(reliable, np.asarray(flags).reshape((2, 2, 2)))
        self.assertFalse(reliable[0, 0, 1])

    def test_mock_complete_fresh_baseline_and_every_round(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'source.json'
            saved = synthetic_source(path)
            runtime = MockRuntime(path)
            with runtime.patched():
                result = mod.run(path, 5120)
            self.assertEqual(result['status'], 'complete')
            self.assertEqual((runtime.coefficient_calls, runtime.count_calls, runtime.fit_calls, runtime.sum_calls), (1, 1, 32, 9))
            self.assertEqual(result['marker_count'], 111)
            self.assertEqual(result['witness']['markers'][:79], saved['witness']['markers'])
            self.assertEqual(result['feasible_full_range'], [1024, 2560])
            self.assertEqual(result['completed_target_count'], 8)
            self.assertIsNone(result['pending_round'])
            self.assertFalse(result['directed'])
            self.assertFalse(result['old_numerical_bounds_reused'])
            for row in result['rounds']:
                self.assertEqual(row['log2_outer_prefix'], 2565.)
                self.assertEqual(row['profile_pairs_evaluated'], 1537**2)
                self.assertEqual(row['selected_reliable_entries'], 6148)

    def test_timeout_keeps_only_completed_round_score(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'source.json'
            synthetic_source(path)
            runtime = MockRuntime(path, timeout_fit=6)
            with runtime.patched():
                result = mod.run(path, 5120)
            self.assertEqual(result['status'], 'timeout')
            self.assertEqual(len(result['rounds']), 2)
            self.assertEqual(result['last_complete_round']['round'], 1)
            self.assertEqual(result['last_complete_round']['marker_count'], 83)
            self.assertEqual(result['pending_round']['round'], 2)
            self.assertTrue(result['witness_may_include_incomplete_round'])
            self.assertNotIn('log2_upper_screen', result)

    def test_mutated_source_and_bad_scope_fail_closed(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'source.json'
            synthetic_source(path)
            runtime = MockRuntime(path, mutate=True)
            with runtime.patched(), self.assertRaisesRegex(RuntimeError, 'artifact changed'):
                mod.run(path, 5120)
            for weight, seconds in ((4096, 90), (8192, 90), (5120, 91)):
                with self.assertRaises(ValueError):
                    mod.run(path, weight, max_seconds=seconds)


if __name__ == '__main__':
    unittest.main()

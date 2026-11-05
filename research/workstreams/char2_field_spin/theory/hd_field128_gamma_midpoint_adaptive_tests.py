"""Tiny original-weight geometry and mocked adaptive calibration; no bank."""

from contextlib import ExitStack
from fractions import Fraction as F
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import hd_field128_gamma_midpoint_adaptive as mod


class Budget:
    def check(self):
        pass


def source(path):
    data = dict(construction=mod.provider.gamma.construction(30), state_dimension=30, target_distance='2/5',
                active_groups=1280, source_sha256={}, complete_selected_weight_screen=True,
                witness=dict(z='1/2', alpha='1/4', grid=128,
                             markers=[[str(F(i+2, 128)), str(F(i+2, 64))] for i in range(111)]))
    path.write_text(json.dumps(data), encoding='utf-8')
    return data


class Runtime:
    def __init__(self, path, *, timeout_sum=None, duplicate=False, mutate=False):
        self.path, self.timeout_sum, self.duplicate, self.mutate = path, timeout_sum, duplicate, mutate
        self.coefficient_calls = self.count_calls = self.bank_calls = self.fit_calls = self.sum_calls = 0
        self.targets = []
        self.coeff = SimpleNamespace(composed=np.array([1., 2., 3.]))

    def coefficients(self, z, alpha, **kwargs):
        self.coefficient_calls += 1
        assert (z, alpha) == (.5, .25) and kwargs['state'] == 30
        return self.coeff

    def bank(self, z, alpha, odds, ratio, **kwargs):
        self.bank_calls += 1
        assert kwargs['_coefficients'] is self.coeff
        assert len(odds) == (111 if self.bank_calls == 1 else 1)
        return SimpleNamespace(score=-10.*self.bank_calls, constant=5., alpha=.25, density=2.,
                               metadata=lambda: dict(synthetic=True, marker_count=len(odds)))

    def regions(self, bank, weight, budget):
        assert weight == 1280
        return np.full((641, 2, 2), bank.score), np.full((641, 2, 2), True)

    def winners(self, bank, weight, budget):
        regions, flags = self.regions(bank, weight, budget)
        return regions, flags, np.zeros(regions.shape, dtype=int)

    def counts(self, pairs, weight, **kwargs):
        self.count_calls += 1
        assert (pairs, weight) == (4096, 1280)
        return object()

    def support(self, pairs, weight, regions, **kwargs):
        self.sum_calls += 1
        assert (pairs, weight) == (4096, 1280) and regions.shape == (641, 2, 2)
        if self.sum_calls == self.timeout_sum:
            raise TimeoutError('synthetic all-profile budget')
        if self.mutate:
            self.path.write_text(self.path.read_text()+' ', encoding='utf-8')
        # The selected target must follow the current sum, not source metadata.
        f = 100 if self.sum_calls == 1 else 200+self.sum_calls
        return dict(profile_pairs_evaluated=641**2, log2_support_sum_screen=float(np.max(regions)),
                    largest_profile=dict(even_full=f, odd_full=f))

    def fit(self, numerator, pairs, weight, full, entry, odds, budget, **kwargs):
        self.fit_calls += 1
        self.targets.append((full, entry))
        positive = (F(2, 128), F(2, 64)) if self.duplicate else (F(100000+self.fit_calls), F(1))
        return positive, dict(entry=entry, histogram=dict(full=full), positive_marker_literals=list(map(str, positive)))

    def patched(self):
        stack = ExitStack()
        for target, name, replacement in ((mod.provider, 'AliasCoefficients', self.coefficients),
                    (mod.provider, 'AliasBankScreen', self.bank), (mod.mode, 'MidpointCounts', self.counts),
                    (mod.mode, 'support_sum', self.support), (mod, 'fit_entry', self.fit),
                    (mod.frozen, 'regional_with_selection', self.regions), (mod, 'regional_winners', self.winners),
                    (mod, 'source_pins', lambda: {})):
            stack.enter_context(patch.object(target, name, replacement))
        stack.enter_context(patch('sys.stdout', new=io.StringIO()))
        return stack


class AdaptiveTests(unittest.TestCase):
    def test_target_entries_and_both_parities(self):
        profile = dict(even_full=100, odd_full=201)
        attempts = {}
        for f in (100, 201):
            for entry in range(4):
                self.assertEqual(mod.choose_target(profile, attempts, pairs=4096, weight=1280), (f, entry))
                attempts[f, entry] = 1
        self.assertEqual(mod.choose_target(profile, attempts, pairs=4096, weight=1280), (100, 0))
        with self.assertRaises(ValueError):
            mod.choose_target(dict(even_full=1023, odd_full=201), {}, pairs=4096, weight=1280)

    def test_dense_winners_match_independent_scalar_min(self):
        odds, tables = [(1., 2.), (0., 1.)], []
        for offset in (0., -3.):
            logs = np.arange(64, dtype=float).reshape((4, 4, 4))+offset
            reliable = np.ones((4, 4, 4), dtype=bool)
            reliable[1, 2, 2] = False
            tables.append(SimpleNamespace(logs=logs, reliable=reliable, positive_numerator=np.full(4, -20.)))
        expected, flags, winners = [], [], []
        lf = mod.provider.alias.LOG_FACTORIAL
        for f in (2, 3):
            s, e = 6-2*f, 4-6+f
            h = lf[4]-lf[e]-lf[s]-lf[f]
            values, chosen_flags, chosen_indices = [], [], []
            for entry in range(4):
                choices = [(float(t.logs[entry, s, f]) if t.reliable[entry, s, f]
                            else float(t.positive_numerator[entry]))-h-s*u-f*v
                           for (u, v), t in zip(odds, tables)]
                i = int(np.argmin(choices))
                values.append(choices[i])
                chosen_flags.append(tables[i].reliable[entry, s, f])
                chosen_indices.append(i)
            expected.append(values)
            flags.append(chosen_flags)
            winners.append(chosen_indices)
        expected = np.asarray(expected).reshape((2, 2, 2))
        bank = SimpleNamespace(grid=4, odds=odds, tables=tables, regional=lambda r: expected.copy())
        actual, reliable, selected = mod.regional_winners(bank, 6, Budget(), pairs=4)
        np.testing.assert_array_equal(actual, expected)
        np.testing.assert_array_equal(reliable, np.asarray(flags).reshape((2, 2, 2)))
        np.testing.assert_array_equal(selected, np.asarray(winners).reshape((2, 2, 2)))
        bank.regional = lambda r: expected+1.
        with self.assertRaisesRegex(ArithmeticError, 'reconstruction'):
            mod.regional_winners(bank, 6, Budget(), pairs=4)

    def test_every_single_addition_gets_complete_original_sum(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'source.json'
            saved = source(path)
            runtime = Runtime(path)
            with runtime.patched():
                result = mod.run(path, 1280)
            self.assertEqual(result['status'], 'marker-limit')
            self.assertEqual((runtime.coefficient_calls, runtime.count_calls, runtime.fit_calls, runtime.sum_calls), (1, 1, 17, 18))
            self.assertEqual(result['marker_count'], 128)
            self.assertEqual(result['witness']['markers'][:111], saved['witness']['markers'])
            self.assertEqual(result['witness']['z'], '1/2')
            self.assertEqual(result['witness']['alpha'], '1/4')
            self.assertEqual(runtime.targets[0], (100, 0))
            self.assertEqual(runtime.targets[1], (202, 0))
            self.assertEqual(result['possible_winner_only_count'], 1)
            self.assertFalse(result['pruning_performed'])
            self.assertEqual(result['winning_marker_indices_at_last_complete'], [127])
            for step, row in enumerate(result['history']):
                self.assertEqual(row['marker_count'], 111+step)
                self.assertEqual(row['profile_pairs_evaluated'], 641**2)
                self.assertEqual(row['log2_outer_prefix'], 645.)
                self.assertEqual(row['selected_reliable_entries'], 2564)

    def test_timeout_exports_only_last_scored_marker_prefix(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'source.json'
            saved = source(path)
            runtime = Runtime(path, timeout_sum=2)
            with runtime.patched():
                result = mod.run(path, 1280)
            self.assertEqual(result['status'], 'timeout')
            self.assertEqual((result['marker_count'], result['applied_marker_count']), (111, 112))
            self.assertEqual(result['witness']['markers'], saved['witness']['markers'])
            self.assertEqual(len(result['unscored_applied_markers']), 1)
            self.assertEqual(len(result['history']), 1)
            self.assertEqual(result['pending_step']['stage'], 'all-profile-sum')
            self.assertEqual(result['winning_marker_indices_at_last_complete'], [0])
            self.assertNotIn('log2_upper_screen', result)

    def test_duplicate_is_not_added_and_fit_budget_still_applies(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'source.json'
            source(path)
            runtime = Runtime(path, duplicate=True)
            with runtime.patched():
                result = mod.run(path, 1280, max_fits=5)
            self.assertEqual((result['marker_count'], runtime.fit_calls, runtime.bank_calls, runtime.sum_calls), (111, 5, 1, 1))
            self.assertEqual(runtime.targets, [(100, 0), (100, 1), (100, 2), (100, 3), (100, 0)])
            self.assertTrue(all(f['addition_status'] == 'duplicate-not-added' for f in result['fits']))

    def test_full_zero_fit_axis_and_budget(self):
        points = []
        def numerator(odds):
            points.append(odds)
            # At r=6 and f=0 the objective becomes (u+2)^2.
            return np.full(4, (odds[0]+2)**2+6*odds[0])
        positive, record = mod.fit_entry(numerator, 16, 6, 0, 2, [(0., 0.)], Budget(), max_calls=64)
        self.assertLessEqual(len(points), 64)
        self.assertEqual(record['objective_calls'], len(points))
        self.assertTrue(all(v == -64. for u, v in points))
        self.assertAlmostEqual(float(positive[0]), .25, places=6)
        self.assertEqual(positive[1], F(1, 2**64))
        self.assertEqual(record['histogram'], dict(empty=10, single=6, full=0))
        with self.assertRaises(ArithmeticError):
            mod.fit_entry(lambda p: np.full(4, np.nan), 16, 6, 0, 0, [], Budget())
        with patch.object(mod.frozen, 'fit_entry', return_value='delegated') as inherited:
            self.assertEqual(mod.fit_entry(numerator, 16, 6, 2, 1, [], Budget()), 'delegated')
            inherited.assert_called_once()

    def test_source_mutation_and_domain_guards(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'source.json'
            source(path)
            runtime = Runtime(path, mutate=True)
            with runtime.patched(), self.assertRaisesRegex(RuntimeError, 'Source witness changed'):
                mod.run(path, 1280, max_fits=0)
            for kwargs in ({'weight': 8192}, {'weight': 1280, 'max_seconds': 91}, {'weight': 1280, 'max_fits': 33}):
                with self.assertRaises(ValueError):
                    mod.run(path, **kwargs)


if __name__ == '__main__':
    unittest.main()

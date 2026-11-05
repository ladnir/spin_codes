"""Mock-only refined-objective search checks; no numerical distance gates."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import hd_field128_gamma_all_pair_search as search


def source():
    return dict(construction='source-literals-only', active_groups=256,
                field_order=1 << 128, log2_upper_screen=-999999., log2_potential_ratio=-8.,
                marker_calibrations=[dict(weight=256)],
                witness=dict(z='1/100000000000000000000', alpha='1/160', grid=64,
                             log2_potential_ratio_exact='-1',
                             marker_log_odds_exact=[['-7/4', '-64'], ['-3/2', '-1/3']]))


class FakeCoefficients:
    created = []

    def __init__(self, z, alpha, *, state, _budget):
        _budget.check()
        self.z, self.alpha, self.state = z, alpha, state
        self.created.append(self)


class FakeBank:
    created = []

    def __init__(self, z, alpha, odds, ratio, *, state, distance, grid, floor_policy,
                 max_table_bytes, _coefficients, _budget):
        _budget.check()
        self.z, self.alpha, self.odds, self.log_ratio = z, alpha, odds, ratio
        self.state, self.distance, self.grid = state, distance, grid
        self.coefficients, self.calls = _coefficients, []
        self.created.append(self)

    def witness(self):
        return dict(z=str(F(self.z)), alpha=str(F(self.alpha)), grid=self.grid,
                    log2_potential_ratio_exact=str(F(self.log_ratio)),
                    marker_log_odds_exact=[[str(F(x)) for x in pair] for pair in self.odds])

    def evaluate(self, *args, **kwargs):
        raise AssertionError('The common-potential objective must never be evaluated')


def selected(bank, r, *, support_method, optimize_potential, _budget):
    _budget.check()
    if support_method != 'exact-profiles' or optimize_potential:
        raise AssertionError('Only the refined no-potential objective is allowed')
    bank.calls.append((r, support_method, optimize_potential))
    return dict(active_groups=r, construction=search.provider.gamma.construction(bank.state),
                state_dimension=bank.state, support_method=support_method,
                log2_upper_screen=bank.alpha*1000, log2_potential_ratio=bank.log_ratio,
                largest_profile=dict(even_full=0, odd_full=r//2, log2_term=2.),
                regional_alias_diagnostics=dict(selected_positive_fallback_entries=0),
                common_potential_used=False, directed=False, whole_code_certificate=False, arithmetic='FLOAT')


def one_point(objective, initial, **kwargs):
    objective(initial)
    return SimpleNamespace(success=True, message='mock complete', nfev=1)


def two_points(objective, initial, **kwargs):
    objective(initial)
    objective([initial[0], initial[1]-1])
    return SimpleNamespace(success=True, message='mock complete', nfev=2)


class RefinedSearchTests(unittest.TestCase):
    def setUp(self):
        FakeCoefficients.created, FakeBank.created = [], []

    def mocked(self):
        return patch.multiple(search.provider, AliasCoefficients=FakeCoefficients,
                              AliasBankScreen=FakeBank, selected_screen=selected)

    def test_initial_marks_exact_source_bank_and_refined_objective_only(self):
        rows = []
        with self.mocked(), patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {}, r=256, state=30, checkpoint=rows.append)
        self.assertEqual(len(FakeBank.created), 1)
        bank = FakeBank.created[0]
        self.assertIs(bank.coefficients, FakeCoefficients.created[0])
        self.assertEqual(bank.z, float(F(source()['witness']['z'])))
        self.assertEqual(bank.alpha, float(F(source()['witness']['alpha'])))
        self.assertEqual(bank.odds, ((-7/4, -64.), (-3/2, -1/3)))
        self.assertEqual(bank.log_ratio, -8.)
        self.assertEqual(bank.calls, [(256, 'exact-profiles', False)])
        self.assertEqual(result['best']['construction'], search.provider.gamma.construction(30))
        self.assertEqual(result['best']['log2_upper_screen'], 6.25)
        self.assertEqual(result['completed_potential_function_calls'], 0)
        self.assertFalse(result['common_potential_used'])
        self.assertFalse(result['old_numerical_bounds_reused'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['marker_bank_exact_source'], source()['witness']['marker_log_odds_exact'])
        self.assertEqual(result['best']['witness']['log2_potential_ratio_exact'], '-8')
        self.assertEqual(rows, [result['best']])
        for data in (result, result['best']):
            replay = search.saved.prepare_proposal(json.loads(json.dumps(data)), {}, r=256, state=30, distance=F(2, 5))
            self.assertEqual(replay['witness']['marker_log_odds_exact'], result['best']['witness']['marker_log_odds_exact'])

    def test_cache_preserves_one_bank_per_distinct_mark(self):
        def repeat(objective, initial, **kwargs):
            first = objective(initial)
            self.assertEqual(first, objective(initial))
            objective([initial[0], initial[1]-1])
            self.assertEqual(kwargs['options']['maxfev'], 3)
            return SimpleNamespace(success=True, message='mock complete', nfev=3)
        with self.mocked(), patch.object(search, 'minimize', side_effect=repeat):
            result = search.run_search(source(), {}, r=256, max_evaluations=3)
        self.assertEqual(result['mark_evaluations'], 2)
        self.assertEqual(result['completed_mark_evaluations'], 2)
        self.assertEqual(len(FakeCoefficients.created), 2)
        self.assertEqual([bank.log_ratio for bank in FakeBank.created], [-8., -8.])
        self.assertEqual(result['best']['evaluation_index'], 1)

    def test_checkpoint_precedes_next_bank_and_timeout_retains_previous(self):
        rows = []
        def timeout(bank, *args, **kwargs):
            if len(FakeBank.created) == 2:
                self.assertEqual(len(rows), 1)
                raise TimeoutError('mock next-mark budget')
            return selected(bank, *args, **kwargs)
        with self.mocked(), patch.object(search.provider, 'selected_screen', side_effect=timeout), \
             patch.object(search, 'minimize', side_effect=two_points):
            result = search.run_search(source(), {}, r=256, checkpoint=rows.append)
        self.assertTrue(result['resource_limit_reached'])
        self.assertEqual(result['best']['evaluation_index'], 0)
        self.assertEqual(result['mark_evaluations'], 2)
        self.assertEqual(result['completed_mark_evaluations'], 1)
        self.assertEqual(result['failures'][0]['error_type'], 'TimeoutError')
        self.assertEqual(len(rows), 1)

    def test_arithmetic_failure_no_result_or_later_recovery(self):
        def failure(bank, *args, **kwargs):
            if len(FakeBank.created) == 1:
                raise ArithmeticError('mock alias issue')
            return selected(bank, *args, **kwargs)
        with self.mocked(), patch.object(search.provider, 'selected_screen', side_effect=failure), \
             patch.object(search, 'minimize', side_effect=two_points):
            result = search.run_search(source(), {}, r=256)
        self.assertEqual(result['best']['evaluation_index'], 1)
        self.assertEqual(result['completed_mark_evaluations'], 1)
        self.assertEqual(result['failures'][0]['error_type'], 'ArithmeticError')
        FakeBank.created.clear()
        with self.mocked(), patch.object(search.provider, 'selected_screen', side_effect=failure), \
             patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {}, r=256)
        self.assertIsNone(result['best'])
        self.assertEqual(result['completed_mark_evaluations'], 0)

    def test_scope_calibration_resource_guards_and_limits(self):
        for override in (dict(state=32), dict(max_evaluations=0), dict(max_seconds=3601), dict(grid=32),
                         dict(max_table_bytes=1), dict(alpha='1/1099511627776')):
            with self.assertRaises(ValueError):
                search.run_search(source(), {}, r=256, **override)
        with self.assertRaises(ValueError):
            search.run_search(source(), {}, r=257, marker_weight=257)
        with self.mocked(), patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {}, r=256, max_seconds=1e-12)
        self.assertIsNone(result['best'])
        self.assertEqual(result['mark_evaluations'], 0)
        self.assertTrue(result['resource_limit_reached'])
        with self.mocked(), patch.object(search, 'minimize', side_effect=two_points):
            result = search.run_search(source(), {}, r=256, max_evaluations=1)
        self.assertEqual(result['mark_evaluations'], 1)
        self.assertEqual(result['completed_mark_evaluations'], 1)
        self.assertEqual(len(FakeBank.created), 1)
        self.assertTrue(result['resource_limit_reached'])

    def test_actual_optimizer_only_mock_banks_and_maxfev(self):
        with self.mocked():
            result = search.run_search(source(), {}, r=256, max_evaluations=4)
        self.assertLessEqual(result['mark_evaluations'], 4)
        self.assertEqual(result['optimizer']['function_calls'], 4)
        self.assertEqual(result['completed_mark_evaluations'], len(FakeBank.created))

    def test_exclusive_outputs_pin_closure_and_source_change(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'evaluation_000_exact_profiles.json'
            record = search.saved.write_new(path, {'scope': 'mock only'})
            self.assertEqual(record['sha256'], sha256(path.read_bytes()).hexdigest())
            with self.assertRaises(FileExistsError):
                search.saved.write_new(path, {})
        pins = search.source_pins()
        for name in ('hd_field128_gamma_all_pair_search.py', 'hd_field128_gamma_all_pair_screen.py',
                     'hd_field128_gamma_profile_diagnostic.py', 'hd_field128_full_rs_gamma_alias_search.py',
                     'hd_field128_full_rs_gamma_screen.py', 'hd_shifted_pairs_count.py'):
            self.assertIn(name, pins)
        rows = []
        with self.mocked(), patch.object(search, 'minimize', side_effect=one_point), \
             patch.object(search, 'source_pins', side_effect=[{'file': 'old'}, {'file': 'changed'}]):
            with self.assertRaises(RuntimeError):
                search.run_search(source(), {}, r=256, checkpoint=rows.append)
        self.assertEqual(rows, [])


if __name__ == '__main__':
    unittest.main()

"""Mock-only search/resource/provenance tests; no distance gates or timings."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import hd_field128_full_rs_state_alias_search as search


def source():
    return dict(construction=search.provider.cover.base.CONSTRUCTION,
                field_order=1 << 128, pair_swaps=True, full_field_random_coefficients=True,
                center=4096, selected_weight=4096, marker_calibrations=[dict(kind='mock-calibrated', target_full='0')],
                state_potential_ratio='1/4', log2_upper_screen=-999999,
                witness=dict(z='1/100000000000000000000', alpha='1/160', grid=64,
                             markers=[['3/4', '5/8'], ['1/2', '1/3']]))


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
        self.z, self.alpha, self.odds, self.log_ratio = z, alpha, odds, ratio
        self.state, self.distance, self.grid = state, distance, grid
        self.coefficients, self.calls = _coefficients, []
        self.created.append(self)

    def evaluate(self, r, *, method, optimize_potential, _budget, max_evaluations=None):
        _budget.check()
        self.calls.append(optimize_potential)
        return dict(active_groups=r, construction=search.provider.local.construction(self.state),
                    state_dimension=self.state, inner_width=64,
                    log2_upper_screen=self.alpha*1000-(50 if optimize_potential else 0),
                    log2_potential_ratio=-123. if optimize_potential else self.log_ratio,
                    potential_evaluations=3 if optimize_potential else 0,
                    regional_alias_diagnostics=dict(entries=8196, selected_positive_fallback_entries=0),
                    covers_all_matching_profiles=True, covers_all_even_odd_splits=True,
                    directed=False, whole_code_certificate=False, arithmetic='FLOAT')

    def witness(self):
        return dict(z=str(F(self.z)), alpha=str(F(self.alpha)), log2_potential_ratio=self.log_ratio,
                    log2_potential_ratio_exact=str(F(self.log_ratio)), marker_log_odds=self.odds)


def one_point(objective, initial, **kwargs):
    objective(initial)
    return SimpleNamespace(success=True, message='mock complete', nfev=1)


class AliasSearchTests(unittest.TestCase):
    def setUp(self):
        FakeCoefficients.created, FakeBank.created = [], []

    def mocked(self):
        return patch.multiple(search.provider, AliasCoefficients=FakeCoefficients, AliasBankScreen=FakeBank)

    def test_source_markers_calibration_and_construction_are_preserved(self):
        artifact = source()
        prepared = search.prepare_proposal(artifact, {'sha256': 'mock'}, r=4096, state=32, distance=F(9, 20))
        self.assertEqual(prepared['witness']['markers'], artifact['witness']['markers'])
        self.assertEqual(prepared['marker_calibrations'], artifact['marker_calibrations'])
        self.assertEqual(prepared['construction'], search.provider.local.construction(32))
        self.assertFalse(prepared['old_numerical_bounds_reused'])
        self.assertNotIn('log2_upper_screen', prepared)
        with self.assertRaises(ValueError):
            search.prepare_proposal(artifact, {}, r=4095, state=32, distance=F(9, 20))

    def test_potential_reuses_bank_and_saved_witness_contains_fitted_ratio(self):
        rows = []
        with self.mocked(), patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {'sha256': 'mock'}, r=4096, checkpoint=rows.append)
        self.assertEqual(len(FakeCoefficients.created), 1)
        self.assertEqual(len(FakeBank.created), 1)
        self.assertIs(FakeBank.created[0].coefficients, FakeCoefficients.created[0])
        self.assertEqual(FakeBank.created[0].calls, [False, True])
        self.assertEqual([row['phase'] for row in rows], ['baseline', 'potential-fit'])
        self.assertEqual(result['best']['log2_potential_ratio'], -123.)
        self.assertEqual(result['best']['witness']['log2_potential_ratio'], -123.)
        self.assertEqual(result['best']['witness']['log2_potential_ratio_exact'], '-123')
        self.assertEqual(result['initial_log2_potential_ratio'], -2.)
        self.assertEqual(result['marker_bank_exact_source'], source()['witness']['markers'])
        self.assertEqual(result['completed_potential_function_calls'], 3)
        self.assertFalse(result['whole_code_certificate'])

    def test_repeated_mark_does_not_rebuild_coefficients_or_bank(self):
        def repeat(objective, initial, **kwargs):
            first = objective(initial)
            self.assertEqual(first, objective(initial))
            objective([initial[0], initial[1]-1])
            return SimpleNamespace(success=True, message='mock complete', nfev=3)
        with self.mocked(), patch.object(search, 'minimize', side_effect=repeat):
            result = search.run_search(source(), {}, r=4096, max_evaluations=3)
        self.assertEqual(result['mark_evaluations'], 2)
        self.assertEqual(len(FakeBank.created), 2)
        self.assertEqual(FakeBank.created[1].log_ratio, -123.)
        self.assertLess(FakeBank.created[1].alpha, FakeBank.created[0].alpha)
        self.assertEqual(result['best']['evaluation_index'], 1)

    def test_timeout_retains_immutable_baseline(self):
        evaluate = FakeBank.evaluate
        def timeout(bank, *args, **kwargs):
            if kwargs['optimize_potential']:
                raise TimeoutError('mock potential budget')
            return evaluate(bank, *args, **kwargs)
        rows = []
        with self.mocked(), patch.object(FakeBank, 'evaluate', new=timeout), \
             patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {}, r=4096, checkpoint=rows.append)
        self.assertTrue(result['resource_limit_reached'])
        self.assertEqual(result['best']['phase'], 'baseline')
        self.assertEqual(len(rows), 1)
        self.assertTrue(result['failures'][0]['completed_baseline_retained'])

    def test_worse_fitted_potential_does_not_replace_baseline(self):
        evaluate = FakeBank.evaluate
        def worse(bank, *args, **kwargs):
            row = evaluate(bank, *args, **kwargs)
            if kwargs['optimize_potential']:
                row['log2_upper_screen'] = 100000.
            return row
        with self.mocked(), patch.object(FakeBank, 'evaluate', new=worse), \
             patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {}, r=4096)
        self.assertEqual(result['best']['phase'], 'baseline')
        self.assertEqual(result['best']['witness']['log2_potential_ratio'], -2.)

    def test_resource_validation_and_timeout_before_first_completed_mark(self):
        for override in (dict(max_evaluations=0), dict(potential_evaluations=257), dict(max_seconds=3601),
                         dict(grid=32), dict(max_table_bytes=1), dict(log_ratio=4097), dict(alpha='1/1099511627776')):
            with self.assertRaises(ValueError):
                search.run_search(source(), {}, r=4096, **override)
        with self.mocked(), patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {}, r=4096, max_seconds=1e-12)
        self.assertIsNone(result['best'])
        self.assertEqual(result['mark_evaluations'], 0)
        self.assertTrue(result['resource_limit_reached'])

    def test_exact_checkpoints_are_exclusive_and_hash_pinned(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'evaluation_000_baseline.json'
            receipt = search.write_new(path, {'directed': False, 'scope': 'mock only'})
            self.assertEqual(receipt['sha256'], sha256(path.read_bytes()).hexdigest())
            self.assertEqual(json.loads(path.read_bytes())['scope'], 'mock only')
            with self.assertRaises(FileExistsError):
                search.write_new(path, {})
        pins = search.source_pins()
        for name in ('hd_field128_full_rs_state_alias_search.py', 'hd_field128_full_rs_state_alias_screen.py',
                     'hd_field128_full_rs_state_screen.py', 'hd_field128_lift_alias_screen.py'):
            self.assertIn(name, pins)


if __name__ == '__main__':
    unittest.main()

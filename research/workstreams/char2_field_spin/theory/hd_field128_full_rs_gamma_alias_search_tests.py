"""Mock-only search, calibration, and checkpoint tests; no numerical gates."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import hd_field128_full_rs_gamma_alias_search as search


def source():
    return dict(construction='old-construction-witness-only', active_groups=256,
                field_order=1 << 128, log2_upper_screen=-999999., log2_potential_ratio=-8.,
                marker_calibrations=[dict(kind='mock-only', weight=256)],
                witness=dict(z='1/100000000000000000000', alpha='1/160', grid=64,
                             log2_potential_ratio=-1., log2_potential_ratio_exact='-1',
                             marker_log_odds=[[-999., -999.]],
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
        self.z, self.alpha, self.odds, self.log_ratio = z, alpha, odds, ratio
        self.state, self.distance, self.grid = state, distance, grid
        self.coefficients, self.calls = _coefficients, []
        self.created.append(self)

    def evaluate(self, r, *, method, optimize_potential, _budget, max_evaluations=None):
        _budget.check()
        self.calls.append(optimize_potential)
        return dict(active_groups=r, construction=search.provider.construction(self.state),
                    state_dimension=self.state, inner_width=64, gamma_extension_degree=16,
                    log2_upper_screen=self.alpha*1000-(50 if optimize_potential else 0),
                    log2_potential_ratio=-123. if optimize_potential else self.log_ratio,
                    potential_evaluations=3 if optimize_potential else 0,
                    regional_alias_diagnostics=dict(entries=8196, selected_positive_fallback_entries=0),
                    covers_all_matching_profiles=True, covers_all_even_odd_splits=True,
                    directed=False, whole_code_certificate=False, arithmetic='FLOAT')

    def witness(self):
        return dict(z=str(F(self.z)), alpha=str(F(self.alpha)), grid=self.grid,
                    log2_potential_ratio=self.log_ratio, log2_potential_ratio_exact=str(F(self.log_ratio)),
                    marker_log_odds=self.odds,
                    marker_log_odds_exact=[[str(F(x)) for x in pair] for pair in self.odds])


def one_point(objective, initial, **kwargs):
    objective(initial)
    return SimpleNamespace(success=True, message='mock complete', nfev=1)


class GammaAliasSearchTests(unittest.TestCase):
    def setUp(self):
        FakeCoefficients.created, FakeBank.created = [], []

    def mocked(self):
        return patch.multiple(search.provider, AliasCoefficients=FakeCoefficients, AliasBankScreen=FakeBank)

    def prepare(self, artifact, **kwargs):
        return search.prepare_proposal(artifact, {'sha256': 'mock'}, r=256, state=24, distance=F(2, 5), **kwargs)

    def test_direct_row_and_summary_preserve_exact_markers_not_bounds(self):
        artifact = source()
        for state in (16, 24, 28, 30):
            for data in (artifact, dict(selected_weight=256, best=artifact)):
                result = search.prepare_proposal(data, {}, r=256, state=state, distance=F(2, 5))
                self.assertEqual(result['construction'], search.provider.construction(state))
                self.assertEqual(result['marker_literals'], artifact['witness']['marker_log_odds_exact'])
                self.assertEqual(result['witness']['z'], artifact['witness']['z'])
                self.assertEqual(result['witness']['alpha'], '1/160')
                self.assertEqual(result['witness']['log2_potential_ratio_exact'], '-8')
                self.assertEqual(result['marker_calibrations'], artifact['marker_calibrations'])
                self.assertEqual(result['source_role'], 'exact witness literals only; source bounds are discarded')
                self.assertNotIn('log2_upper_screen', result)
                self.assertNotIn('log2_upper_screen', result['witness'])
                self.assertFalse(result['old_numerical_bounds_reused'])
        prepared = self.prepare(artifact)
        self.assertEqual(prepared['witness']['marker_log_odds_exact'][0], ['-7/4', '-64'])
        self.assertEqual(float(F(prepared['witness']['marker_log_odds_exact'][1][1])), -1/3)

    def test_positive_exact_markers_and_potential_are_supported(self):
        artifact = source()
        del artifact['log2_potential_ratio']
        for key in ('marker_log_odds_exact', 'log2_potential_ratio', 'log2_potential_ratio_exact'):
            del artifact['witness'][key]
        artifact['witness']['markers'] = [['3/4', '5/8'], ['1/2', '1/3']]
        artifact['state_potential_ratio'] = '1/4'
        result = self.prepare(artifact)
        self.assertEqual(result['marker_representation'], 'positive-exact')
        self.assertEqual(result['marker_literals'], artifact['witness']['markers'])
        self.assertEqual(result['witness']['log2_potential_ratio_exact'], '-2')
        self.assertEqual(float(F(result['witness']['marker_log_odds_exact'][1][0])), -1.)

    def test_calibration_must_be_explicit_and_cannot_be_overridden(self):
        missing = source()
        del missing['active_groups']
        with self.assertRaises(ValueError):
            self.prepare(missing)
        self.assertEqual(self.prepare(missing, marker_weight=256)['source_weight_calibration'],
                         {'explicit_marker_weight': 256})
        cases = [dict(source(), selected_weight=255), dict(selected_weight=255, best=source()),
                 dict(source(), active_groups=True), dict(best=None), dict(source(), field_order=256)]
        for artifact in cases:
            with self.assertRaises(ValueError):
                self.prepare(artifact, marker_weight=256)
        wrong = source()
        wrong['witness']['center'] = 255
        with self.assertRaises(ValueError):
            self.prepare(wrong)
        with self.assertRaises(ValueError):
            search.prepare_proposal(source(), {}, r=256, state=32, distance=F(2, 5))

    def test_missing_inexact_duplicate_or_invalid_markers_are_rejected(self):
        for replacement in (None, [[-1., '-2']], [['-1', '-2'], ['-1', '-2']], [['4097', '0']], [['0']]):
            artifact = source()
            if replacement is None:
                del artifact['witness']['marker_log_odds_exact']
            else:
                artifact['witness']['marker_log_odds_exact'] = replacement
            with self.assertRaises(ValueError):
                self.prepare(artifact)

    def test_one_bank_per_mark_fitted_witness_and_roundtrip(self):
        rows = []
        with self.mocked(), patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {}, r=256, checkpoint=rows.append)
        self.assertEqual(len(FakeCoefficients.created), 1)
        self.assertEqual(len(FakeBank.created), 1)
        self.assertIs(FakeBank.created[0].coefficients, FakeCoefficients.created[0])
        self.assertEqual(FakeBank.created[0].calls, [False, True])
        self.assertEqual(FakeBank.created[0].z, float(F(source()['witness']['z'])))
        self.assertEqual(FakeBank.created[0].alpha, float(F(source()['witness']['alpha'])))
        self.assertEqual(FakeBank.created[0].odds, ((-7/4, -64.), (-3/2, -1/3)))
        self.assertEqual([row['phase'] for row in rows], ['baseline', 'potential-fit'])
        self.assertEqual(result['best']['witness']['log2_potential_ratio_exact'], '-123')
        self.assertEqual(result['completed_potential_function_calls'], 3)
        self.assertEqual(result['marker_bank_exact_source'], source()['witness']['marker_log_odds_exact'])
        self.assertFalse(result['whole_code_certificate'])
        for data in (result, result['best'], rows[0]):
            replay = self.prepare(json.loads(json.dumps(data)))
            expected = data.get('best', data)['witness']
            self.assertEqual(replay['witness']['marker_log_odds_exact'], expected['marker_log_odds_exact'])
            self.assertEqual(replay['witness']['log2_potential_ratio_exact'], expected['log2_potential_ratio_exact'])

    def test_repeated_marks_are_cached_and_next_mark_uses_best_potential(self):
        def repeat(objective, initial, **kwargs):
            first = objective(initial)
            self.assertEqual(first, objective(initial))
            objective([initial[0], initial[1]-1])
            return SimpleNamespace(success=True, message='mock complete', nfev=3)
        with self.mocked(), patch.object(search, 'minimize', side_effect=repeat):
            result = search.run_search(source(), {}, r=256, max_evaluations=3)
        self.assertEqual(result['mark_evaluations'], 2)
        self.assertEqual(len(FakeBank.created), 2)
        self.assertEqual(FakeBank.created[1].log_ratio, -123.)
        self.assertEqual(result['best']['evaluation_index'], 1)

    def test_checkpoint_precedes_fitting_and_timeout_keeps_baseline(self):
        evaluate, rows = FakeBank.evaluate, []
        def timeout(bank, *args, **kwargs):
            if kwargs['optimize_potential']:
                self.assertEqual(len(rows), 1)
                raise TimeoutError('mock potential budget')
            return evaluate(bank, *args, **kwargs)
        with self.mocked(), patch.object(FakeBank, 'evaluate', new=timeout), \
             patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {}, r=256, checkpoint=rows.append)
        self.assertTrue(result['resource_limit_reached'])
        self.assertEqual(result['best']['phase'], 'baseline')
        self.assertTrue(result['failures'][0]['completed_baseline_retained'])
        self.assertEqual(result['best']['marker_literals'], source()['witness']['marker_log_odds_exact'])

    def test_worse_fitting_keeps_better_baseline(self):
        evaluate = FakeBank.evaluate
        def worse(bank, *args, **kwargs):
            row = evaluate(bank, *args, **kwargs)
            if kwargs['optimize_potential']:
                row['log2_upper_screen'] = 100000.
            return row
        with self.mocked(), patch.object(FakeBank, 'evaluate', new=worse), \
             patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {}, r=256)
        self.assertEqual(result['best']['phase'], 'baseline')
        self.assertEqual(result['best']['witness']['log2_potential_ratio_exact'], '-8')

    def test_resource_validation_and_time_and_evaluation_limits(self):
        for override in (dict(max_evaluations=0), dict(potential_evaluations=257), dict(max_seconds=3601),
                         dict(grid=32), dict(max_table_bytes=1), dict(log_ratio=4097),
                         dict(alpha='1/1099511627776')):
            with self.assertRaises(ValueError):
                search.run_search(source(), {}, r=256, **override)
        with self.mocked(), patch.object(search, 'minimize', side_effect=one_point):
            result = search.run_search(source(), {}, r=256, max_seconds=1e-12)
        self.assertIsNone(result['best'])
        self.assertEqual(result['mark_evaluations'], 0)
        self.assertTrue(result['resource_limit_reached'])
        def exceed(objective, initial, **kwargs):
            objective(initial)
            objective([initial[0], initial[1]-1])
        with self.mocked(), patch.object(search, 'minimize', side_effect=exceed):
            result = search.run_search(source(), {}, r=256, max_evaluations=1)
        self.assertEqual(result['mark_evaluations'], 1)
        self.assertEqual(len(FakeBank.created), 1)
        self.assertTrue(result['resource_limit_reached'])
        self.assertIsNotNone(result['best'])

    def test_exclusive_checkpoints_recursive_pins_and_source_change_guard(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'evaluation_000_baseline.json'
            receipt = search.write_new(path, {'scope': 'mock only'})
            self.assertEqual(receipt['sha256'], sha256(path.read_bytes()).hexdigest())
            with self.assertRaises(FileExistsError):
                search.write_new(path, {})
        pins = search.source_pins()
        for name in ('hd_field128_full_rs_gamma_alias_search.py', 'hd_field128_full_rs_gamma_screen.py',
                     'hd_field128_full_rs_state_alias_screen.py', 'hd_field128_full_rs_state_screen.py',
                     'hd_field128_lift_alias_screen.py'):
            self.assertIn(name, pins)
        rows = []
        with self.mocked(), patch.object(search, 'minimize', side_effect=one_point), \
             patch.object(search, 'source_pins', side_effect=[{'file': 'old'}, {'file': 'changed'}]):
            with self.assertRaises(RuntimeError):
                search.run_search(source(), {}, r=256, checkpoint=rows.append)
        self.assertEqual(rows, [])


if __name__ == '__main__':
    unittest.main()

"""Tiny independent ordinary-matrix controls; no construction-size run."""
from itertools import product
import io
import json
from pathlib import Path
import unittest
from unittest.mock import Mock, patch

import numpy as np
import hd_field128_translated_full_quartet_v54 as pilot


def logs(values):
    with np.errstate(divide='ignore'):
        return np.log2(values)


def brute(matrices, weights, alpha):
    """Independent ordinary positive products; no log-product helper."""
    two, four = np.zeros((2, 2)), np.zeros((2, 2))
    for i, j in product(range(len(matrices)), repeat=2):
        two += weights[i]*weights[j]*(matrices[i]@matrices[j])**alpha
    for i, j, k, ell in product(range(len(matrices)), repeat=4):
        four += (weights[i]*weights[j]*weights[k]*weights[ell]
                 *(((matrices[i]@matrices[j])@matrices[k])@matrices[ell])**alpha)
    return two, four


class Checks(unittest.TestCase):
    def toy(self):
        matrices = np.array([[[1., 2.], [3., 4.]], [[.5, 3.], [1., .25]],
                             [[2., .5], [.125, 1.5]]])
        raw = logs(np.array([[a@b for b in matrices] for a in matrices]))
        return matrices, raw

    def test_independent_ordered_sums_and_weights_outside_alpha(self):
        matrices, raw = self.toy()
        weights = np.array([[2., 3., 5.], [.25, 7., 1.]])
        for alpha in (.125, .5, 1.):
            actual = pilot.full_quartets(raw, logs(weights), alpha, chunk=2)
            for i, row in enumerate(weights):
                two, four = brute(matrices, row, alpha)
                np.testing.assert_allclose(np.exp2(actual['two'][i]), two, rtol=3.e-14)
                np.testing.assert_allclose(np.exp2(actual['four'][i]), four, rtol=4.e-14)
            self.assertEqual(actual['ordered_quartets_per_marker'], 81)
            self.assertEqual(actual['raw_quartets_composed'], 81)
            self.assertEqual(actual['pair_blocks'], 5)

    def test_alpha_one_is_exact_matrix_power(self):
        matrices, raw = self.toy()
        weights = np.array([2., 3., 5.])
        actual = pilot.full_quartets(raw, logs(weights), 1.)
        single = np.sum(weights[:, None, None]*matrices, axis=0)
        np.testing.assert_allclose(np.exp2(actual['four'][0]), np.linalg.matrix_power(single, 4), rtol=3.e-14)
        np.testing.assert_allclose(actual['four'], actual['old_four'], atol=2.e-14, rtol=0)

    def test_chunk_independence_and_complete_pair_blocks(self):
        _, raw = self.toy()
        weights = logs([[2., 3., 5.], [7., 1., .125]])
        expected = pilot.full_quartets(raw, weights, .2, chunk=1)
        for chunk in (2, 4, 8, 32):
            actual = pilot.full_quartets(raw, weights, .2, chunk=chunk)
            np.testing.assert_allclose(actual['four'], expected['four'], atol=2.e-14, rtol=0)
            self.assertEqual(actual['ordered_quartets_per_marker'], 3**4)

    def test_count_formula_without_construction_size_arithmetic(self):
        self.assertEqual(pilot.quartet_count(65), 17_850_625)
        self.assertEqual(pilot.quartet_count(16), 65_536)
        for bad in (0, 66, True, 3.5):
            with self.assertRaises(ValueError):
                pilot.quartet_count(bad)

    def test_exact_zero_logs_and_zero_marker_weights(self):
        matrices = np.array([np.eye(2), np.array([[0., 2.], [0., 0.]])])
        raw = logs(np.array([[a@b for b in matrices] for a in matrices]))
        for weights in (np.array([1., 3.]), np.array([1., 0.]), np.array([0., 0.])):
            actual = pilot.full_quartets(raw, logs(weights), .5, chunk=1)
            two, four = brute(matrices, weights, .5)
            np.testing.assert_allclose(np.exp2(actual['two'][0]), two, atol=0, rtol=2.e-14)
            np.testing.assert_allclose(np.exp2(actual['four'][0]), four, atol=0, rtol=3.e-14)
            self.assertFalse(np.any(np.isnan(actual['four'])))
        zero = pilot.full_quartets(np.full((2, 2, 2, 2), -np.inf), [0., 0.], .125)
        self.assertTrue(np.all(np.isneginf(zero['four'])))

    def test_boundary_floor_and_input_guards(self):
        _, raw = self.toy()
        actual = pilot.full_quartets(raw, [0., 1., 2.], .125)
        self.assertTrue(np.all(actual['four'] <= actual['old_four']+1.e-13))
        self.assertTrue(np.all(actual['four'] >= actual['optimistic_floor']-1.e-13))
        for kwargs in (dict(alpha=0), dict(alpha=1.1), dict(alpha=.5, chunk=0),
                       dict(alpha=.5, chunk=33), dict(alpha=float('nan'))):
            with self.assertRaises(ValueError):
                pilot.full_quartets(raw, [0., 1., 2.], **kwargs)
        for invalid in (np.full_like(raw, np.inf), np.full_like(raw, np.nan), raw[:, :2]):
            with self.assertRaises(ValueError):
                pilot.full_quartets(invalid, [0., 1., 2.], .5)
        with self.assertRaises(ValueError):
            pilot.full_quartets(raw, np.zeros((5, 3)), .5)
        with self.assertRaises(TimeoutError):
            pilot.full_quartets(raw, [0., 1., 2.], .5, check=Mock(side_effect=TimeoutError))

    def fixture(self):
        markers = [['1', '1'], ['2', '1'], ['1', '2'], ['2', '2']]
        common = dict(active_groups=4096, target_distance='9/20', conditional_on_principal_setup=True,
                      directed=False, whole_code_certificate=False, setup_risk_added=False)
        v52 = dict(common, schema='reversed-witness-exchange-v52-FLOAT', all_profiles_retained=True,
            final=dict(profiles_evaluated=2049**2,
                       largest_profile=dict(even_full=1613, odd_full=1727, witness=0)),
            witness_specs=[dict(z='1/4', alpha='1/128', markers=markers)])
        v53 = dict(common, schema='translated-quartet-diagnostic-v53-FLOAT', selected_profile=[1613, 1727],
            witness_index=0, z='1/4', alpha='1/128',
            input=dict(path=pilot.INPUT52, sha256=pilot.INPUT52_SHA),
            markers=[dict(marker_index=i, marker=pair) for i, pair in enumerate(markers)])
        return v52, v53

    def test_fixed_input_scope_and_literal_controls(self):
        a, b = self.fixture()
        self.assertEqual(pilot.validate_inputs(a, b)[1], [0, 1, 2, 3])
        for item, key, value in ((a, 'active_groups', 4095), (b, 'target_distance', '2/5'),
                                 (b, 'witness_index', 1), (b, 'z', '1/2'),
                                 (b, 'conditional_on_principal_setup', False)):
            old = item[key]
            item[key] = value
            with self.assertRaises(ValueError):
                pilot.validate_inputs(a, b)
            item[key] = old
        b['markers'][0]['marker'] = ['3', '1']
        with self.assertRaises(ValueError):
            pilot.validate_inputs(a, b)
        with self.assertRaises(ArithmeticError):
            pilot.close_control(1., 2., 'mock')

    def test_write_limits_without_files_and_hard_guard_entry(self):
        path = pilot.HERE/'never-written-v54-test.json'
        with patch.object(pilot.shutil, 'disk_usage', return_value=Mock(free=pilot.RESERVE)), \
                patch.object(Path, 'open') as opening:
            with self.assertRaises(RuntimeError):
                pilot.write_result(path, {'ok': True})
            opening.assert_not_called()
        with patch.object(pilot.shutil, 'disk_usage', return_value=Mock(free=pilot.RESERVE+pilot.LIMIT)), \
                patch.object(Path, 'open') as opening:
            with self.assertRaises(RuntimeError):
                pilot.write_result(path, {'oversized': 'x'*pilot.LIMIT})
            opening.assert_not_called()
        buffer, context = io.BytesIO(), Mock()
        context.__enter__ = Mock(return_value=buffer)
        context.__exit__ = Mock(return_value=False)
        with patch.object(pilot.shutil, 'disk_usage', return_value=Mock(free=pilot.RESERVE+pilot.LIMIT)), \
                patch.object(Path, 'open', return_value=context) as opening, \
                patch.object(Path, 'read_bytes', side_effect=lambda: buffer.getvalue()):
            receipt = pilot.write_result(path, {'ok': True})
            opening.assert_called_once_with('xb')
            self.assertEqual(json.loads(buffer.getvalue()), {'ok': True})
            self.assertEqual(receipt['bytes'], len(buffer.getvalue()))
        with patch.object(pilot.guard, 'remote_limits', return_value=context), \
                patch.object(pilot, '_run', return_value={'mock': True}) as run:
            self.assertEqual(pilot.run('output'), {'mock': True})
            run.assert_called_once_with('output')


if __name__ == '__main__':
    unittest.main()

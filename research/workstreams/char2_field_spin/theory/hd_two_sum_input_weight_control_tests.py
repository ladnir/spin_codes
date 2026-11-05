"""Lightweight scope and normalization checks; no regional bank is built."""
from hashlib import sha256
import json
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_input_weight_control as target


class InputWeightControlTests(unittest.TestCase):
    def fixture(self):
        return json.dumps(dict(banks=[dict(active_groups=r, mark=str(mark),
            alpha=str(target.low.ALPHA), curves=[[0.]*(r+1)]*4,
            initial_log2_slopes=[0.]*4) for r, _, mark in target.low.ANCHORS])).encode()

    def test_empty_floor_template_reaches_declared_target(self):
        alpha = float(target.low.ALPHA)
        arities = np.arange(1, 17)
        for r, _, mark in target.low.ANCHORS:
            templates = dict(target.marker_templates(r, mark))
            self.assertEqual(set(templates), {'zero', 'uniform1', 'uniform2',
                                              'empty-floor1', 'empty-floor2'})
            prefix = alpha*(-target.low.fractional.LIMIT*log2(float(mark))-log2(255))
            empty = np.array([log2(comb(16, a))-8*a*alpha*log2(255) for a in arities])
            for name in ('empty-floor1', 'empty-floor2'):
                values = empty-9*arities*alpha*templates[name]
                bound = prefix+log2(comb(8192, r))+r*target.state.logsum(values)
                self.assertLessEqual(bound, -60.+1e-9)

    def test_preparation_keeps_markers_fixed_per_arity(self):
        slopes = np.array([-.1, -.3, .2, -.4])
        markers = np.array([0., 2.]*8)
        def weights(alpha, marker):
            self.assertEqual(alpha, target.low.ALPHA)
            return {a: .01*a+np.arange(17)*(.02*a+marker) for a in range(1, 17)}
        with patch.object(target.joint, 'type_weights', side_effect=weights) as builder:
            prepared = target.prepared_joint(np.zeros((4, 4)), slopes, markers)
        self.assertEqual(builder.call_count, 2)
        for index, row in enumerate(prepared['classes']):
            positions = np.repeat(slopes, row[2])
            expected = [.01*a+np.logaddexp2(0., positions+.02*a+markers[a-1]).sum()
                        for a in range(1, 17)]
            np.testing.assert_allclose(prepared['typed'][index], expected, atol=1e-12, rtol=0)
        with self.assertRaises(ValueError):
            target.prepared_joint(np.zeros((4, 4)), slopes, [-1.]*16)

    def test_complete_mock_run_has_exactly_fifteen_fixed_controls(self):
        raw = self.fixture()
        prepared = {'intercepts': np.zeros(4)}
        with patch.object(target.Path, 'read_bytes', return_value=raw), \
                patch.object(target, 'BANK_SHA', sha256(raw).hexdigest()), \
                patch.object(target.low, 'checked_dependencies', return_value={}), \
                patch.object(target, 'prepared_joint', return_value=prepared) as prepare, \
                patch.object(target, 'unmarked_bound', return_value=123.), \
                patch.object(target.exact, 'evaluate_prepared',
                             return_value={'failure_log2_screen': 99.}) as exact, \
                patch.object(target, 'perf_counter', return_value=0.):
            result = target.run(1)
        self.assertTrue(result['complete_requested_controls'])
        self.assertFalse(result['directed'])
        self.assertFalse(result['high_side_included'])
        self.assertEqual(result['regional_banks_rebuilt'], 0)
        self.assertEqual(prepare.call_count, 15)
        self.assertEqual(exact.call_count, 15)
        self.assertEqual(len(result['witnesses']), 15)
        for index, (r, cutoff, mark) in enumerate(target.low.ANCHORS):
            rows = result['witnesses'][5*index:5*index+5]
            self.assertEqual({row['active_groups'] for row in rows}, {r})
            self.assertEqual({row['maximum_total_arity'] for row in rows}, {cutoff-1})
            self.assertEqual({row['mark'] for row in rows}, {str(mark)})
            for call in exact.call_args_list[5*index:5*index+5]:
                self.assertEqual(call.args[1:], (r, cutoff, mark, target.low.ALPHA))

    def test_timeout_preserves_completed_controls(self):
        raw = self.fixture()
        with patch.object(target.Path, 'read_bytes', return_value=raw), \
                patch.object(target, 'BANK_SHA', sha256(raw).hexdigest()), \
                patch.object(target.low, 'checked_dependencies', return_value={}), \
                patch.object(target, 'prepared_joint', return_value={'intercepts': np.zeros(4)}), \
                patch.object(target, 'unmarked_bound', return_value=123.), \
                patch.object(target.exact, 'evaluate_prepared', side_effect=[
                    {'failure_log2_screen': 99.}, TimeoutError('bounded fixture')]), \
                patch.object(target, 'perf_counter', return_value=0.):
            result = target.run(1)
        self.assertEqual(result['status'], 'compute-budget-exceeded')
        self.assertFalse(result['complete_requested_controls'])
        self.assertEqual(len(result['witnesses']), 1)

    def test_pin_and_anchor_mismatch_fail_before_preparation(self):
        raw = self.fixture()
        with patch.object(target.Path, 'read_bytes', return_value=raw), \
                patch.object(target, 'prepared_joint') as prepare:
            with self.assertRaises(ValueError):
                target.run(1)
            prepare.assert_not_called()
        data = json.loads(raw)
        data['banks'][0]['alpha'] = '1/4'
        raw = json.dumps(data).encode()
        with patch.object(target.Path, 'read_bytes', return_value=raw), \
                patch.object(target, 'BANK_SHA', sha256(raw).hexdigest()), \
                patch.object(target.low, 'checked_dependencies', return_value={}), \
                patch.object(target, 'prepared_joint') as prepare:
            with self.assertRaises(ValueError):
                target.run(1)
            prepare.assert_not_called()

    def test_invalid_budget_does_not_read_bank(self):
        with patch.object(target.Path, 'read_bytes') as read:
            for budget in (0, -1, 61, True, float('nan'), float('inf')):
                with self.assertRaises(ValueError):
                    target.run(budget)
            read.assert_not_called()


if __name__ == '__main__':
    unittest.main()

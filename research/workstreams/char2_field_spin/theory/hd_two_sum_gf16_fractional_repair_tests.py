"""Small repair-scope controls; no regional provider or production search."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_fractional_repair as target


def baseline(index=0):
    mark, alpha = target.source.CONTROLS[index//5]
    return dict(mark=str(mark), alpha=str(alpha), bank_index=index//5,
        log2_floor_markers=['0']*16, log2_slopes=[1., 2., 3., 4.],
        floor_template='fixture', floor_anchor=None, floor_scale=0,
        failure_log2_screen=35.)


def fake_evaluate(curves, weights, mark, alpha, floors, slopes):
    value = 10+(slopes[3]+1)**2
    return dict(mark=str(mark), alpha=str(alpha), log2_floor_markers=list(map(str, floors)),
        log2_slopes=list(slopes), failure_log2_screen=value,
        entries=[dict(active_groups=r, failure_log2_screen=value) for r in range(18, 33)])


class FractionalRepairTests(unittest.TestCase):
    def test_actual_frozen_pins_scope_and_exact_floor_literals(self):
        data, banks = target.load_frozen()
        self.assertEqual(len(banks), 5)
        self.assertTrue(all(bank.shape == (4, 33) for bank in banks))
        self.assertEqual(len(data['witnesses'][:25]), 25)
        self.assertEqual(data['witnesses'][25]['coordinate_rounds'], 1)
        for row in data['witnesses'][:25]:
            self.assertTrue(all(F(value) >= 0 for value in row['log2_floor_markers']))
        self.assertTrue(np.isneginf(banks[0][1, 0]))

    def test_hash_mismatch_rejected_without_parsing(self):
        raw = b'not a proof artifact'
        with patch.object(Path, 'read_bytes', return_value=raw):
            self.assertEqual(target.read_pinned('fixture', sha256(raw).hexdigest()), raw)
            with self.assertRaises(ValueError):
                target.read_pinned('fixture', '0'*64)

    def test_only_live_coordinate_changes_and_zero_candidate_is_tested(self):
        seen = []
        def evaluate(*args):
            seen.append(list(args[-1]))
            return fake_evaluate(*args)
        candidates = ((0.,), (0.,), (0.,), (-2., -1., 0., 4.))
        with patch.object(target.state, 'slope_candidates', return_value=(np.zeros(4), candidates)), \
                patch.object(target.source, 'evaluate', side_effect=evaluate):
            row = target.repair_live_slope(np.zeros((4, 33)), {}, baseline(), lambda: None)
        self.assertEqual(row['log2_slopes'], [1., 2., 3., -1.])
        self.assertTrue(all(slopes[:3] == [1., 2., 3.] for slopes in seen))
        self.assertEqual({slopes[3] for slopes in seen}, {-2., -1., 0., 4.})
        self.assertEqual(row['live_slope_evaluations'], 4)
        self.assertEqual(row['original_failure_log2_screen'], 35.)
        self.assertEqual(row['log2_floor_markers'], ['0']*16)

    def test_stale_baseline_rejected_and_candidate_budget_enforced(self):
        candidates = ((0.,), (0.,), (0.,), (-1., 0.))
        bad = dict(baseline(), failure_log2_screen=0.)
        with patch.object(target.state, 'slope_candidates', return_value=(np.zeros(4), candidates)), \
                patch.object(target.source, 'evaluate', side_effect=fake_evaluate):
            with self.assertRaises(ArithmeticError):
                target.repair_live_slope(np.zeros((4, 33)), {}, bad, lambda: None)
        candidates = ((0.,), (0.,), (0.,), tuple(range(17)))
        with patch.object(target.state, 'slope_candidates', return_value=(np.zeros(4), candidates)), \
                patch.object(target.source, 'evaluate') as evaluate:
            with self.assertRaises(ArithmeticError):
                target.repair_live_slope(np.zeros((4, 33)), {}, baseline(), lambda: None)
        evaluate.assert_not_called()

    def test_twenty_five_repairs_then_only_one_global_coordinate_round(self):
        frozen = dict(witnesses=[baseline(i) for i in range(25)]+[{'not_a_baseline': True}])
        counter = []
        def repair(curves, weights, original, check):
            check()
            index = len(counter)
            counter.append(index)
            return dict(original, failure_log2_screen=-float(index),
                live_slope_evaluations=4,
                entries=[dict(active_groups=r, failure_log2_screen=-float(index)) for r in range(18, 33)])
        def refine(curves, weights, original, check):
            check()
            self.assertEqual(original['failure_log2_screen'], -24.)
            return dict(original, slope_evaluations=1)
        with patch.object(target, 'load_frozen', return_value=(frozen, [np.zeros((4, 33))]*5)), \
                patch.object(target.source, 'fractional_dense_weights', return_value={}) as weights, \
                patch.object(target.source, 'fractional_regions') as provider, \
                patch.object(target, 'repair_live_slope', side_effect=repair) as repair_call, \
                patch.object(target.source, 'refine_one_round', side_effect=refine) as refine_call:
            result = target.run()
        self.assertEqual(weights.call_count, 5)
        self.assertEqual(repair_call.call_count, 25)
        self.assertEqual(refine_call.call_count, 1)
        provider.assert_not_called()
        self.assertEqual(result['best_repaired_baseline_index'], 24)
        self.assertEqual(result['repaired_baseline_count'], 25)
        self.assertEqual(len(result['witnesses']), 26)
        self.assertTrue(result['complete_requested_band'])
        self.assertEqual(result['regional_banks_rebuilt'], 0)
        self.assertFalse(result['directed'])
        self.assertFalse(result['extra_witness_factor'])

    def test_guard_and_timeout_do_not_load_or_rebuild_banks(self):
        for kwargs in ({'max_seconds': 31}, {'max_seconds': True}, {'max_seconds': float('nan')}, {'refine': 1}):
            with patch.object(target, 'load_frozen') as load, self.assertRaises(ValueError):
                target.run(**kwargs)
            load.assert_not_called()
        with patch.object(target, 'load_frozen') as load:
            result = target.run(max_seconds=1e-12)
        load.assert_not_called()
        self.assertFalse(result['complete_requested_band'])
        self.assertEqual(result['status'], 'compute-budget-exceeded')


if __name__ == '__main__':
    unittest.main()

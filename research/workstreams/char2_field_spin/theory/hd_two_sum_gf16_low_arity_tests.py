"""Exact toy low-arity restrictions and bounded orchestration, without banks."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_low_arity as target


class LowArityTests(unittest.TestCase):
    def test_exact_whole_word_arity_marker_bound(self):
        masses = {1: F(3), 2: F(5), 3: F(7)}
        for groups in (1, 2, 3):
            for bound in range(groups, 3*groups+1):
                exact = sum((np.prod([masses[a] for a in types])
                             for types in product(masses, repeat=groups) if sum(types) <= bound), F(0))
                for t in (F(1), F(1, 2), F(1, 16)):
                    majorant = t**(-bound)*sum((value*t**a for a, value in masses.items()), F(0))**groups
                    self.assertLessEqual(exact, majorant)
                if bound == 3*groups:
                    self.assertEqual(exact, sum(masses.values())**groups)

    def test_column_tilt_is_outside_alpha_and_applies_once_per_whole_type(self):
        typed = np.zeros((2, 16))
        for log_t in target.LOG_T:
            tilted = target.tilt_type_logs(typed, log_t)
            np.testing.assert_array_equal(16*tilted[0], np.arange(1, 17)*log_t)
        with self.assertRaises(ValueError):
            target.tilt_type_logs(typed, 1)

    def test_low_type_union_under_noncommuting_path_kernel(self):
        classes = target.state.path_classes(3)
        path_a = np.array([sum(n*v for n, v in zip(row[2], (0., 1., 2., 0.))) for row in classes])
        typed = np.full((len(classes), 16), -np.inf)
        for index, (_, _, counts, _) in enumerate(classes):
            typed[index, :3] = np.log2([1+counts[0], 2+counts[1], 1+counts[3]])
        exact = 0
        groups, columns, bound = 2, 2, 4
        for types in product(range(3), repeat=groups):
            if sum(a+1 for a in types) > bound:
                continue
            matrix = np.zeros((2, 2), dtype=object)
            for index, (initial, final, _, count) in enumerate(classes):
                matrix[initial, final] += count*round(2**(path_a[index]+sum(typed[index, a] for a in types)))
            exact += sum(np.linalg.matrix_power(matrix, columns)[0])
        for log_t in (0, -1, -4):
            kernel = target.state.kernel_from_types(classes, path_a,
                target.tilt_type_logs(typed, log_t, columns), groups, columns)
            power = target.state.matrix_power_log(kernel, columns)
            upper = 2**(-bound*log_t+float(target.state.logsum(power[0])))
            self.assertGreaterEqual(upper, exact*(1-1e-13))

    def test_peak_tail_initialization_avoids_deep_empty_endpoint(self):
        curve = np.array([-500., -500., -20., -22., -26., -30., -40., -60., -80.])
        curves = np.stack((curve, curve-1, curve-2, curve-3))
        initial, candidates, points = target.peak_tail_slopes(curves)
        self.assertTrue(np.all(initial < 0))
        self.assertTrue(all(row['nonempty_peak'] == 2 for row in points))
        self.assertTrue(all(0. in row and len(row) <= 7 for row in candidates))
        self.assertTrue(all(0 not in row['secant_positions'] for row in points))
        intercepts, slack = target.state.affine_envelope(curves, initial)
        self.assertTrue(np.all(slack >= -1e-12))
        self.assertTrue(np.all(np.isfinite(intercepts)))

    def test_actual_anchor_floors_and_frozen_scope(self):
        for r, a0, mark in target.ANCHORS:
            self.assertEqual(target.checked_low_scope(r, a0, mark), a0-1)
            templates = target.floor_templates(r, mark)
            self.assertEqual(len(templates), 3)
            for row in templates[1:]:
                self.assertIn(f'anchor{r}', row['name'])
                self.assertTrue(all(isinstance(x, F) and x >= 0 for x in row['log2_floor_markers']))
                empty = target.fractional.all_empty_log(mark, target.ALPHA, row['log2_floor_markers'], r)
                self.assertLessEqual(empty, -60+1e-8)
        for args in ((18, 248, F(691, 700)), (19, 247, F(691, 700)), (18, 247, F(99, 100))):
            with self.assertRaises(ValueError):
                target.checked_low_scope(*args)

    def test_outer_arity_correction_not_fractionally_scaled(self):
        prepared = dict(classes=(), path_a=np.zeros(0), typed=np.zeros((0, 16)),
            floors=[F(0)]*16, slopes=np.zeros(4), intercepts=np.zeros(4), slack=np.zeros((4, 19)))
        power = np.array([[0., -np.inf], [-np.inf, 0.]])
        with patch.object(target.state, 'kernel_from_types', return_value=power), \
                patch.object(target.state, 'matrix_power_log', return_value=power):
            got = target.evaluate_prepared(prepared, *target.ANCHORS[0], -4)
        r, a0, mark = target.ANCHORS[0]
        expected = .2*(-target.fractional.LIMIT*log2(float(mark))-log2(255))+log2(comb(8192, r))+4*(a0-1)
        self.assertAlmostEqual(got['failure_log2_screen'], expected)
        self.assertEqual(got['outside_total_arity_log_factor'], 984)
        self.assertEqual(got['total_arity_marker'], '1/16')

    def test_three_banks_forty_five_baselines_and_only_optional_worst_refinement(self):
        curves = np.tile(-np.arange(33), (4, 1)).astype(float)
        calls = []
        def evaluate(prepared, r, a0, mark, log_t):
            calls.append((r, log_t))
            return dict(active_groups=r, high_arity_cutoff=a0, low_total_active_constituents_max=a0-1,
                mark=str(mark), alpha=str(target.ALPHA), log2_total_arity_marker=log_t,
                log2_slopes=[-1.]*4, failure_log2_screen=float(r-log_t))
        def refine(curves, weights, anchor, baseline, templates, candidates, check):
            check()
            self.assertEqual(anchor[0], 32)
            return dict(baseline, joint_evaluations=1)
        with patch.object(target, 'checked_dependencies', return_value={}), \
                patch.object(target.fractional, 'fractional_dense_weights', return_value={}), \
                patch.object(target.fractional, 'fractional_regions') as bank, \
                patch.object(target.state, 'regional_curves', return_value=curves), \
                patch.object(target, 'prepared_envelope', return_value={}), \
                patch.object(target, 'evaluate_prepared', side_effect=evaluate), \
                patch.object(target, 'one_joint_round', side_effect=refine) as round_call:
            result = target.run(refine_worst=True)
        self.assertEqual(bank.call_count, 3)
        self.assertTrue(all(call.args[2] == 32 for call in bank.call_args_list))
        self.assertEqual(len(calls), 45)
        self.assertEqual(round_call.call_count, 1)
        self.assertEqual(result['refined_anchor_index'], 2)
        self.assertTrue(result['complete_requested_anchors'])
        self.assertFalse(result['continuous_weight_band_covered'])
        self.assertFalse(result['high_side_included'])
        self.assertFalse(result['directed'])

    def test_guard_and_timeout_do_not_build_a_bank(self):
        for kwargs in ({'max_seconds': 31}, {'max_seconds': True}, {'max_seconds': float('nan')}, {'refine_worst': 1}):
            with patch.object(target.fractional, 'fractional_regions') as bank, self.assertRaises(ValueError):
                target.run(**kwargs)
            bank.assert_not_called()
        with patch.object(target.fractional, 'fractional_regions') as bank:
            result = target.run(max_seconds=1e-12)
        bank.assert_not_called()
        self.assertFalse(result['complete_requested_anchors'])
        self.assertEqual(result['status'], 'compute-budget-exceeded')


if __name__ == '__main__':
    unittest.main()

"""Small exact/combinatorial checks; no production packet bank is built."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_fractional_state as target
import hd_two_sum_gf16_state_envelope as old
from hd_two_sum_closed_kernel import arity_power_sum, uniform_mask_count
from hd_two_sum_gf16_factorized import dense_arity_logs


class FractionalStateTests(unittest.TestCase):
    def test_closed_moments_match_every_fixed_rowset_in_toy_group(self):
        for alpha in (F(1), F(1, 16), F(3, 10)):
            for arity in range(3, 9):
                active_sets = [mask for mask in range(256) if mask.bit_count() == arity]
                for weight in range(9):
                    mask = (1 << weight)-1
                    counts = [uniform_mask_count(active, mask, 8) for active in active_sets]
                    logs = [16*float(alpha)*log2(n) for n in counts if n]
                    expected = float(np.logaddexp2.reduce(logs)) if logs else -np.inf
                    got = target.fractional_power_sum_log(arity, weight, alpha, group=8)
                    if np.isneginf(expected):
                        self.assertTrue(np.isneginf(got))
                    else:
                        self.assertAlmostEqual(got, expected, places=11)

    def test_integer_fractional_moments_and_overlap_edge_arities(self):
        for alpha in (F(1), F(1, 2), F(1, 4), F(1, 8), F(1, 16)):
            exponent = int(16*alpha)
            for arity in range(3, 17):
                for weight in range(17):
                    value = arity_power_sum(arity, weight, 16, exponent)
                    got = target.fractional_power_sum_log(arity, weight, alpha)
                    if value:
                        self.assertAlmostEqual(got, log2(value), places=11)
                    else:
                        self.assertTrue(np.isneginf(got))
        # A_g=all coordinates gives N_A(T)=255^|T| exactly.
        weights = target.fractional_dense_weights(F(1, 4))
        np.testing.assert_allclose(weights[16], (np.arange(17)-8)*log2(255)/4, atol=1e-12)
        self.assertAlmostEqual(weights[3][0], log2(comb(16, 3))/16-3*log2(255)/8)

    def test_alpha1_dense_weights_match_old_but_fractional_weights_are_not_old_power(self):
        expected, _ = dense_arity_logs()
        got = target.fractional_dense_weights(F(1))
        for arity in range(3, 17):
            np.testing.assert_allclose(got[arity], expected[arity], atol=1e-12)
        half = target.fractional_dense_weights(F(1, 2))
        self.assertAlmostEqual(half[3][0]-.5*got[3][0], log2(comb(16, 3))/32)
        self.assertGreater(half[3][0]-.5*got[3][0], .1)

    def test_packet_power_precedes_geometry_average_and_zero_stays_zero(self):
        left, right = (1, 0, 0, 0), (0, 1, 0, 0)
        table = {left: np.array([F(1, 16), F(1, 4), 0, 1], dtype=float),
                 right: np.array([F(9, 16), 1, F(1, 4), 0], dtype=float)}
        geometry = [(1, 0, 3, left), (1, 0, 1, right)]
        marker = .25
        got = target.fractional_packet_coefficients(table, geometry, F(1, 2), 2, marker)
        expected = marker*(1.5*np.sqrt(table[left])+.5*np.sqrt(table[right]))
        np.testing.assert_allclose(got[:, 1, 0], expected)
        self.assertNotAlmostEqual(got[0, 1, 0], marker*np.sqrt(1.5*table[left][0]+.5*table[right][0]))
        alpha1 = target.fractional_packet_coefficients(table, geometry, 1, 2, marker)
        np.testing.assert_allclose(alpha1[:, 1, 0], marker*(1.5*table[left]+.5*table[right]))
        only = target.fractional_packet_coefficients(table, geometry[:1], F(1, 2), 2, marker)
        self.assertEqual(only[2, 1, 0], 0.)

    def test_route_normalizer_is_unpowered(self):
        maximum, pairs, marker = 2, 7, .25
        powered = np.zeros((4, 3, 2))
        expected = np.array([.125, .25, .5, 1.])
        powered[:, 1, 0] = comb(pairs, 1)*marker*expected
        powered[:, 0, 1] = comb(pairs, 1)*marker*expected
        result = target.normalize_regional_coefficients(powered, maximum, marker, pairs)
        np.testing.assert_allclose(2**result[1, 0].ravel(), expected)
        np.testing.assert_allclose(2**result[2, 1].ravel(), expected)
        self.assertTrue(np.all(np.isneginf(result[0, 0])))

    def test_fractional_sparse_and_dense_type_masses_against_labelled_masks(self):
        q, group, columns, dimension, minimum = 7, 4, 5, 2, 3
        alpha = F(1, 2)
        floors = np.array([2., 3., 5., 7.])
        marks = np.array([2., 3., 5., 7.])
        dense = {3: np.array([1., 2., 3., 4., 5.]), 4: np.array([2., 1., 4., 3., 6.])}
        got = 2**target.fractional_type_weights(np.log2(marks),
            {a: np.log2(w) for a, w in dense.items()}, np.log2(floors), alpha,
            columns=columns, dimension=dimension, minimum=minimum, q=q)
        full = np.prod(marks)
        af = float(alpha)
        one = (group/(q-1)**(dimension*af)/floors[0]**minimum)**(1/columns)
        one *= 1+(q-1)**af*floors[0]*full
        two = (comb(group, 2)/(q-1)**(2*dimension*af)/floors[1]**(2*minimum))**(1/columns)
        two *= (1+2*(q-1)**af*floors[1]*full
                +floors[1]**2*(((q-1)*(q-1-group))**af*full
                            +(q-1)**af*sum(full/b for b in marks)))
        expected = [one, two]
        for arity in (3, 4):
            floor = floors[arity-1]
            value = sum(dense[arity][sum(mask)]*np.prod(marks**mask)
                        *(floor if any(mask) else 1) for mask in product(range(2), repeat=group))
            expected.append(value/floor**(minimum/columns))
        np.testing.assert_allclose(got, expected, rtol=2e-14)

    def test_alpha1_entire_evaluation_matches_state_helper_at_identical_slopes(self):
        m = np.arange(33)
        curves = np.array([-m, -4.-m/2, -8.-m/3, -2.-m/4])
        curves[1:3, 0] = -np.inf
        weights = target.fractional_dense_weights(1)
        slopes, _ = old.slope_candidates(curves)
        expected = old.evaluate(curves, weights, F(99, 100), 0, slopes)
        got = target.evaluate(curves, weights, F(99, 100), 1, 0, slopes)
        np.testing.assert_allclose([r['failure_log2_screen'] for r in got['entries']],
                                   [r['failure_log2_screen'] for r in expected['entries']], atol=1e-10)
        self.assertAlmostEqual(got['failure_log2_screen'], expected['failure_log2_screen'], places=9)

    def test_type_specific_templates_suppress_all_empty_mass_at_their_anchor(self):
        for mark, alpha in target.CONTROLS:
            rows = target.floor_templates(mark, alpha)
            self.assertEqual(len(rows), 5)
            self.assertTrue(all(value == 0 for value in rows[0]['log2_floor_markers']))
            for row in rows[1:]:
                values = row['log2_floor_markers']
                self.assertEqual(len(values), 16)
                self.assertTrue(all(isinstance(value, F) and value >= 0 for value in values))
                self.assertLessEqual(target.all_empty_log(mark, alpha, values, row['anchor']), -60+1e-8)
            self.assertTrue(any(value.denominator != 1 for value in rows[1]['log2_floor_markers']))
        with self.assertRaises(ValueError):
            target.floor_array([0]*15)
        with self.assertRaises(ValueError):
            target.floor_array([-1]+[0]*15)

    def test_fixed_twenty_five_baselines_and_only_one_global_refinement(self):
        curves = np.zeros((4, 33))
        def fake_evaluate(c, w, mark, alpha, floors, slopes):
            value = float(mark)-float(alpha)+sum(map(float, floors))
            return dict(mark=str(mark), alpha=str(alpha), log2_floor_markers=list(map(str, floors)),
                log2_slopes=list(slopes), failure_log2_screen=value,
                entries=[dict(active_groups=r, failure_log2_screen=value) for r in range(18, 33)])
        def fake_refine(c, w, baseline, check):
            check()
            return dict(baseline, coordinate_rounds=1, slope_evaluations=1)
        with patch.object(target, 'fractional_regions') as bank, \
                patch.object(target, 'regional_curves', return_value=curves), \
                patch.object(target, 'fractional_dense_weights', return_value={}), \
                patch.object(target, 'evaluate', side_effect=fake_evaluate) as evaluate, \
                patch.object(target, 'refine_one_round', side_effect=fake_refine) as refine:
            result = target.run()
        self.assertEqual(bank.call_count, 5)
        self.assertEqual(evaluate.call_count, 25)
        self.assertEqual(refine.call_count, 1)
        self.assertEqual(result['best_baseline_index'], 0)
        self.assertEqual(len(result['witnesses']), 26)
        self.assertTrue(result['complete_requested_band'])
        self.assertFalse(result['directed'])
        self.assertFalse(result['route_normalizer_powered'])
        self.assertEqual(result['projective_normalizations'], 1)
        self.assertEqual(result['endpoint_potential_factors'], 0)

    def test_scope_guards_and_timeout_without_hidden_bank(self):
        for alpha in (0, -1, F(3, 2)):
            with self.assertRaises(ValueError):
                target.checked_alpha(alpha)
        for kwargs in ({'max_seconds': 31}, {'max_seconds': True}, {'max_seconds': float('nan')}, {'refine': 1}):
            with patch.object(target, 'fractional_regions') as bank, self.assertRaises(ValueError):
                target.run(**kwargs)
            bank.assert_not_called()
        with patch.object(target, 'fractional_regions') as bank:
            result = target.run(max_seconds=1e-12)
        bank.assert_not_called()
        self.assertFalse(result['complete_requested_band'])
        self.assertEqual(result['status'], 'compute-budget-exceeded')


if __name__ == '__main__':
    unittest.main()

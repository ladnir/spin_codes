"""Tiny matrix/histogram and memory-only guards; no construction-size work."""
from fractions import Fraction as F
from itertools import product
import io
import json
from math import comb, log2
from pathlib import Path
import unittest
from unittest.mock import Mock, patch

import numpy as np
import hd_field128_gamma_four_packet_positive_v35 as pilot


class Checks(unittest.TestCase):
    def toy(self):
        matrices = np.array([[[1., 2.], [3., 4.]], [[.5, 3.], [1., .25]],
                             [[2., .5], [.125, 1.5]]])
        weights = np.array([2., 3., 5.])
        logs = np.log2(matrices)
        return matrices, weights, pilot.compose(logs[:, None], logs[None])

    def test_positive_partition_matches_brute_ordered_quartets(self):
        matrices, weights, raw = self.toy()
        selected = (0, 2)
        for alpha in (.125, .5, 1.):
            actual = pilot.positive_partition(raw, np.log2(weights), alpha, selected)
            c2, rest, c4, full = (np.zeros((2, 2)) for _ in range(4))
            for i, j in product(range(3), repeat=2):
                contribution = weights[i]*weights[j]*(matrices[i]@matrices[j])**alpha
                if i in selected and j in selected:
                    c2 += contribution
                else:
                    rest += contribution
            for i, j, k, ell in product(range(3), repeat=4):
                multiplicity = weights[i]*weights[j]*weights[k]*weights[ell]
                if all(index in selected for index in (i, j, k, ell)):
                    contribution = multiplicity*(matrices[i]@matrices[j]@matrices[k]@matrices[ell])**alpha
                    c4 += contribution
                else:
                    contribution = multiplicity*((matrices[i]@matrices[j])**alpha
                                                 @ (matrices[k]@matrices[ell])**alpha)
                full += contribution
            for key, expected in (('selected_two', c2), ('remainder_two', rest),
                                  ('selected_four', c4), ('four', full), ('two', c2+rest)):
                np.testing.assert_allclose(np.exp2(actual[key]), expected, rtol=2e-14, atol=2e-11)
            self.assertEqual(actual['quartets_evaluated'], 16)

    def test_alpha_one_and_boundary_inequality(self):
        matrices, weights, raw = self.toy()
        self.assertFalse(np.allclose(matrices[0]@matrices[1], matrices[1]@matrices[0]))
        for selected in ((0,), (1, 2), (0, 1, 2)):
            for alpha in (.01, .25, 1.):
                actual = pilot.positive_partition(raw, np.log2(weights), alpha, selected)
                self.assertTrue(np.all(actual['four'] <= actual['old_four']+1e-12))
                self.assertTrue(np.all(actual['four'] >= actual['old_four']-(1-alpha)-1e-12))
                if alpha == 1:
                    np.testing.assert_allclose(actual['four'], actual['old_four'], rtol=0, atol=1e-12)
                if len(selected) == 3:
                    self.assertTrue(np.all(np.isneginf(actual['remainder_two'])))

    def test_histogram_mass_outside_alpha_and_no_swap_multiplier(self):
        x, y, slots = 3, 5, 3
        actual = np.exp2(pilot.occupancy_weights((log2(x), log2(y)), slots))
        expected = np.zeros(2*slots+1)
        for assignment in product(range(3), repeat=slots):
            single, full = assignment.count(1), assignment.count(2)
            expected[single+2*full] += x**single*y**full
        np.testing.assert_allclose(actual, expected, rtol=1e-14)
        self.assertAlmostEqual(actual.sum(), (1+x+y)**slots)
        self.assertEqual(pilot.occupancy_weights((0., 0.), 1).tolist(), [0., 0., 0.])

    def test_original_histogram_normalizer_and_block_power(self):
        two = np.log2(np.array([[.5, .25], [.125, .75]]))
        odds = (log2(3), log2(5))
        # Toy region:4pair slots,4active coordinates,1full,2single,1empty.
        old = pilot.regional(two, 2, 4, 1, odds, pairs=4)
        new = pilot.regional(pilot.compose(two, two), 1, 4, 1, odds, pairs=4)
        expected = np.linalg.matrix_power(np.exp2(two), 2)/(comb(4, 2)*comb(2, 1)*3**2*5)
        np.testing.assert_allclose(np.exp2(old), expected, rtol=1e-14)
        np.testing.assert_allclose(new, old, rtol=0, atol=1e-14)
        # Dense original-weight geometry must not silently replace r by2P-r.
        dense = pilot.regional(two, 2, 6, 2, odds, pairs=4)
        expected_dense = np.linalg.matrix_power(np.exp2(two), 2)/(comb(4, 2)*3**2*5**2)
        np.testing.assert_allclose(np.exp2(dense), expected_dense, rtol=1e-14)
        with self.assertRaises(ValueError):
            pilot.regional(two, 2, 7, 1, odds, pairs=4)

    def test_conditional_state_and_source_projection(self):
        row = dict(pilot.previous.state_metadata(32), active_groups=4096, target_distance='9/20',
                   profile_pairs_evaluated=2049**2, complete_selected_weight_screen=True,
                   all_profiles_retained=True, witness=dict(z='1/4', alpha='1/128', markers=[['2', '3']]),
                   largest_profile=dict(even_full=2026, odd_full=2026, log2_count_upper=688.),
                   log2_outer_prefix=589316., log2_upper_screen=4071.)
        artifact = dict(schema='gamma-state-ladder-conditional-FLOAT-v34', active_groups=4096,
                        target_distance='9/20', source_sha256={'mock': 'pin'},
                        fixed_positive_marker_literals=[['2', '3']], best_by_state={'32': row})
        with patch.object(pilot.previous, 'source_pins', return_value={'mock': 'pin'}):
            self.assertIs(pilot.validate_source(artifact), row)
            row['conditional_on_global_principal_nonsingularity'] = False
            with self.assertRaises(ValueError):
                pilot.validate_source(artifact)
            row['conditional_on_global_principal_nonsingularity'] = True
            row['state_dimension'] = 30
            with self.assertRaises(ValueError):
                pilot.validate_source(artifact)

    def test_selection_is_bounded_and_partition_support_retained(self):
        _, weights, raw = self.toy()
        terms = pilot.pair_terms(raw, np.log2(weights), .125)
        selected = pilot.choose_occupancies(terms, 2)
        self.assertEqual(len(selected), 2)
        actual = pilot.positive_partition(raw, np.log2(weights), .125, selected)
        np.testing.assert_allclose(np.logaddexp2(actual['selected_two'], actual['remainder_two']),
                                   actual['two'], rtol=0, atol=1e-14)
        for bad in ((0, 0), (3,), ()):
            with self.assertRaises(ValueError):
                pilot.positive_partition(raw, np.log2(weights), .125, bad)

    def test_budget_and_numpy_output_without_disk_writes(self):
        path = pilot.HERE/'never-written-v35-test.json'
        with patch.object(pilot.shutil, 'disk_usage', return_value=Mock(free=pilot.RESERVE)), \
                patch.object(Path, 'open') as opening:
            with self.assertRaises(RuntimeError):
                pilot.write_result(path, {'n': np.int64(2)})
            opening.assert_not_called()
        buffer = io.BytesIO()
        context = Mock()
        context.__enter__ = Mock(return_value=buffer)
        context.__exit__ = Mock(return_value=False)
        with patch.object(pilot.shutil, 'disk_usage', return_value=Mock(free=pilot.RESERVE+pilot.OUTPUT_LIMIT)), \
                patch.object(Path, 'open', return_value=context), \
                patch.object(Path, 'read_bytes', side_effect=lambda: buffer.getvalue()):
            result = pilot.write_result(path, {'n': np.int64(2)})
        self.assertEqual(json.loads(buffer.getvalue()), {'n': 2})
        self.assertEqual(result['bytes'], len(buffer.getvalue()))
        with self.assertRaises(ValueError):
            pilot._run('source', 'not-approved', 'out')

    def test_entry_installs_hard_guard(self):
        context = Mock()
        context.__enter__ = Mock(return_value=None)
        context.__exit__ = Mock(return_value=False)
        with patch.object(pilot, 'remote_limits', return_value=context), \
                patch.object(pilot, '_run', return_value={'mock': True}) as run:
            self.assertEqual(pilot.run('s', 'h', 'o'), {'mock': True})
        context.__enter__.assert_called_once()
        run.assert_called_once_with('s', 'h', 'o')


if __name__ == '__main__':
    unittest.main()

from collections import Counter
from itertools import combinations
import json
from math import comb
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from hd_shifted_b128_screen import (geometry, sparse_profile_counts, sparse_regions,
                                    sparse_weight_screen, regional_alias_log_table,
                                    medium_weight_screen, LOG_DENSITY)
from hd_shifted_outer_density import log_density_float, TOEPLITZ_OUTER_MEASURE
from hd_shifted_swap_sparse import regional_screen


class SmallerGroupTests(unittest.TestCase):
    def test_geometry_and_density_are_consistent_without_old_default_change(self):
        g = geometry()
        self.assertEqual(g['outer_groups']*g['input_symbols_per_group'], 1 << 20)
        self.assertEqual(g['regions']*g['packets_per_region']*64, 1 << 21)
        self.assertEqual((g['pairs_per_region'], g['paired_blocks_per_region'], g['regional_pair_repetitions']),
                         (8192, 128, 64))
        self.assertEqual(LOG_DENSITY*2, log_density_float(TOEPLITZ_OUTER_MEASURE))
        self.assertFalse(g['directed'])

    def test_sparse_joint_counts_exhaust_all_small_cyclic_supports(self):
        for groups in (4, 6, 8, 10):
            for weight in (1, 2):
                observed = Counter()
                for positions in combinations(range(groups), weight):
                    support = set(positions)
                    a = sum(i in support and i+1 in support for i in range(0, groups, 2))
                    b = sum(i in support and (i+1) % groups in support for i in range(1, groups, 2))
                    observed[a, b] += 1
                self.assertEqual(dict(observed), sparse_profile_counts(groups, weight))
                self.assertEqual(sum(observed.values()), comb(groups, weight))

    def test_b256_control_matches_existing_sparse_normalization(self):
        actual = sparse_regions(.999, .25, 2, pairs=4096, packet_steps=128)
        expected = regional_screen(.999, .25, 2)
        for key in expected:
            np.testing.assert_allclose(actual[key], expected[key], rtol=0, atol=1.e-12)

    def test_alpha_one_pairing_equals_unpaired_two_hundred_fifty_six_steps(self):
        paired = sparse_regions(.999, 1., 2)
        single = sparse_regions(.999, 1., 2, paired=False)
        for key in paired:
            np.testing.assert_allclose(paired[key], single[key], rtol=0, atol=2.e-12)

    def test_unweighted_identity_coefficients_normalize_to_identity(self):
        coefficients = np.zeros((4, 3, 2))
        for s, f in ((0, 0), (1, 0), (2, 0), (0, 1)):
            coefficients[(0, 3), s, f] = comb(32, s)*comb(32-s, f)
        with patch('hd_shifted_b128_screen.swap_single_screen', return_value=coefficients):
            regions = sparse_regions(.9, 1., 2, paired=False)
        for region in regions.values():
            np.testing.assert_allclose(np.exp2(region), np.eye(2), rtol=0, atol=2.e-10)

    def test_zero_start_and_ordered_pair_power_are_retained(self):
        regions = {(2, 0): ((0., -1.), (-2., -3.)), (2, 1): ((-4., -2.), (-1., -5.))}
        with patch('hd_shifted_b128_screen.sparse_regions', return_value=regions), \
                patch('hd_shifted_b128_screen.power', return_value=((0., -1.), (100., 100.))) as powered:
            result = sparse_weight_screen(2, .999, .1)
        self.assertEqual(powered.call_count, 3)
        self.assertTrue(all(call.args[1] == 64 for call in powered.call_args_list))
        self.assertNotEqual(powered.call_args_list[1].args[0], powered.call_args_list[2].args[0])
        self.assertEqual(sum(row['support_count'] for row in result['profiles']), comb(16384, 2))
        self.assertFalse(result['directed'])

    def test_regional_alias_power_and_scaling_on_constant_polynomial(self):
        matrix = np.array([[.5, .125], [.25, .375]])
        def frequencies(*args):
            return np.broadcast_to(matrix.ravel()[:, None, None], (4, 64, 64)).astype(complex).copy()
        with patch('hd_shifted_b128_screen.swap_coefficients', return_value=np.zeros((4, 33, 33))), \
                patch('hd_shifted_b128_screen.good_step_coefficients', return_value=None), \
                patch('hd_shifted_b128_screen.good_pair_coefficients', return_value=None), \
                patch('hd_shifted_b128_screen.low_degree_frequency', return_value=None), \
                patch('hd_shifted_b128_screen.paired_frequency', side_effect=frequencies):
            for blocks in (64, 128):
                table, diagnostic = regional_alias_log_table(.9, .1, (-2., -3.),
                                                              grid=64, paired_blocks=blocks)
                expected = np.log2(np.linalg.matrix_power(matrix, blocks)).ravel()
                np.testing.assert_allclose(table[:, 0, 0], expected, rtol=0, atol=1.e-9)
                self.assertEqual(diagnostic['paired_blocks'], blocks)

    def test_medium_gate_is_fixed_at_one_weight_and_new_geometry(self):
        from hd_shifted_swap_screen import CONSTRUCTION
        exact = dict(active_groups=1024, z='1/2', alpha='1/10', potential='1/4', markers=['1/2 1/4'])
        artifact = dict(construction=CONSTRUCTION, fixed_witness=dict(construction=CONSTRUCTION,
            pair_swaps=True, orientation_aggregation='independent-fair-swaps', grid=128,
            source_exact_witness=exact))
        with patch('hd_shifted_b128_screen.regional_alias_log_table',
                   return_value=(np.zeros((4, 128, 128)), {})) as atlas, \
                patch('hd_shifted_b128_screen.run_log_screen_float', return_value=-100.) as run:
            result = medium_weight_screen(artifact)
        self.assertEqual(atlas.call_count, 1)
        self.assertEqual(atlas.call_args.kwargs['paired_blocks'], 128)
        self.assertEqual(run.call_args.args[:2], (8192, 2048))
        self.assertEqual(run.call_args.kwargs['repetitions'], 128)
        self.assertEqual(len(run.call_args.args[2]), 1025)
        self.assertEqual(result['exact_frozen_witness'], exact)
        self.assertEqual(result['witness_fits'], 0)
        self.assertFalse(result['directed'])
        exact['active_groups'] = 512
        with self.assertRaises(ValueError):
            medium_weight_screen(artifact)

    def test_optional_frozen_gates_are_float_only_and_did_not_fit(self):
        paths = [Path(__file__).with_name(name) for name in (
            'hd_shifted_b128_sparse_screen.json', 'hd_shifted_b128_medium2048_screen.json')]
        if not all(path.exists() for path in paths):
            self.skipTest('Optional structural screens absent; tests never rebuild the medium bank')
        sparse, medium = [json.loads(path.read_text(encoding='utf-8')) for path in paths]
        for result in (sparse, medium):
            self.assertFalse(result['directed'])
            self.assertFalse(result['whole_code_certificate'])
            self.assertEqual(result['outer_groups'], 16384)
            self.assertEqual(result['regions'], 128)
            json.dumps(result, allow_nan=False)
        self.assertTrue(sparse['all_pass_40_bits'])
        self.assertEqual([row['active_groups'] for row in sparse['results']], [1, 2])
        self.assertEqual(sparse['new_fits'], 0)
        self.assertEqual(medium['active_groups'], 2048)
        self.assertEqual(medium['source_weight'], 1024)
        self.assertEqual(medium['atlases_built'], 22)
        self.assertEqual(medium['witness_fits'], 0)
        self.assertGreater(medium['log2_upper_screen'], 4000.)


if __name__ == '__main__':
    unittest.main()

import json
from contextlib import redirect_stdout
from io import StringIO
from fractions import Fraction as F
from pathlib import Path
import runpy
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from hd_shifted_bulk_search import LOG_DENSITY, LOG_FACTORIAL
from hd_shifted_potential_band_screen import FixedAtlasScreen
from hd_shifted_outer_density import TOEPLITZ_OUTER_MEASURE, log_density_float
from hd_shifted_swap_screen import CONSTRUCTION


class FixedAtlasTests(unittest.TestCase):
    def tables(self, z, alpha, odds, grid):
        return np.full((4, grid, grid), 5. if odds == (1., 2.) else 0.)

    def test_tables_reused_and_minimum_follows_each_normalization(self):
        with patch('hd_shifted_potential_band_screen.alias_table', side_effect=self.tables) as atlas:
            screen = FixedAtlasScreen(.5, .1, ((0, 0), (1, 2), (0, 0)), .25, grid=64)
            first = screen.regional(10)
            screen.regional(11)
        self.assertEqual(atlas.call_count, 2)
        multinomial = LOG_FACTORIAL[4096]-LOG_FACTORIAL[4089]-LOG_FACTORIAL[4]-LOG_FACTORIAL[3]
        np.testing.assert_allclose(first[3], np.full((2, 2), -5.-multinomial))
        self.assertEqual(tuple(first), tuple(range(6)))

    def test_complete_band_has_no_sampling_and_one_common_potential(self):
        ratios, weights = [], []
        with patch('hd_shifted_potential_band_screen.alias_table', side_effect=self.tables):
            screen = FixedAtlasScreen(.5, .1, ((0, 0),), -.125, grid=64)
        def fake_run(pairs, weight, regions, ratio):
            weights.append(weight)
            ratios.append(ratio)
            self.assertEqual(set(regions), set(range(weight//2+1)))
            return -screen.alpha*LOG_DENSITY*weight-screen.constant-(20 if weight == 2 else 100)
        with patch('hd_shifted_potential_band_screen.run_log_screen', side_effect=fake_run):
            result = screen.band(1, 3, max_weights=3)
        self.assertEqual(weights, [1, 2, 3])
        self.assertEqual(ratios, [-.125]*3)
        self.assertEqual(result['failed_weights'], [2])
        self.assertTrue(result['every_weight_evaluated'])
        self.assertFalse(result['all_weights_pass_threshold'])
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['witness_searches'], 0)
        self.assertAlmostEqual(result['union_log2_upper_screen'], -20.)
        self.assertNotIn('outer_measure', result)
        self.assertNotIn('outer_measure', result['fixed_witness'])

    def test_sharper_outer_measure_keeps_atlas_and_changes_only_prefix(self):
        with patch('hd_shifted_potential_band_screen.alias_table', side_effect=self.tables) as atlas:
            legacy = FixedAtlasScreen(.5, .1, ((0, 0),), 0., grid=64)
            sharp = FixedAtlasScreen(.5, .1, ((0, 0),), 0., grid=64,
                                     outer_measure=TOEPLITZ_OUTER_MEASURE)
        np.testing.assert_array_equal(legacy.tables[0], sharp.tables[0])
        self.assertEqual(atlas.call_args_list[0], atlas.call_args_list[1])
        with patch('hd_shifted_potential_band_screen.run_log_screen', return_value=-100.):
            low, high = sharp.evaluate(16), legacy.evaluate(16)
            receipt = sharp.band(16, 16)
        expected = 1.6*(LOG_DENSITY-log_density_float(TOEPLITZ_OUTER_MEASURE))
        self.assertAlmostEqual(high['log2_upper_screen']-low['log2_upper_screen'], expected, delta=2.e-11)
        self.assertEqual(receipt['fixed_witness']['outer_measure'], TOEPLITZ_OUTER_MEASURE)
        self.assertEqual(receipt['outer_measure_bound']['density_formula'], '256^160/255^32')
        with self.assertRaises(ValueError):
            FixedAtlasScreen(.5, .1, ((0, 0),), 0., outer_measure='unproved')

    def test_dense_histogram_start_is_not_assumed_zero(self):
        with patch('hd_shifted_potential_band_screen.alias_table', side_effect=self.tables):
            screen = FixedAtlasScreen(.5, .1, ((0, 0),), 0., grid=64)
        regions = screen.regional(4097)
        self.assertNotIn(0, regions)
        self.assertEqual(min(regions), 1)
        self.assertEqual(max(regions), 2048)
        expected = -(LOG_FACTORIAL[4096]-LOG_FACTORIAL[4095])
        np.testing.assert_allclose(regions[1], np.full((2, 2), expected))

    def test_checkpoint_loads_fixed_literals_without_search(self):
        from hd_shifted_swap_potential import SAVED_RUN_WITNESSES
        with patch('hd_shifted_potential_band_screen.alias_table', side_effect=self.tables) as atlas:
            screen = FixedAtlasScreen.from_checkpoint(1024)
        self.assertEqual(atlas.call_count, len(SAVED_RUN_WITNESSES[1024]['markers']))
        self.assertEqual(screen.witness()['source_exact_witness']['potential'],
                         SAVED_RUN_WITNESSES[1024]['potential'])
        with self.assertRaises(ValueError):
            FixedAtlasScreen.from_checkpoint(123)

    def test_exact_witness_constructor_rejects_raw_only_and_preserves_literals(self):
        witness = dict(construction=CONSTRUCTION, pair_swaps=True,
            orientation_aggregation='independent-fair-swaps', grid=64,
            mark=.9, alpha=.8, marker_log_odds=((9., 8.),), log2_potential_ratio=7.)
        with self.assertRaises(ValueError):
            FixedAtlasScreen.from_witness(witness)
        witness['source_exact_witness'] = dict(z='1/2', alpha='1/10', potential='2/3', markers=['1/2 1/8'])
        with patch('hd_shifted_potential_band_screen.alias_table', side_effect=self.tables) as atlas:
            screen = FixedAtlasScreen.from_witness(witness, outer_measure=TOEPLITZ_OUTER_MEASURE)
        self.assertEqual((screen.z, screen.alpha, screen.odds), (.5, .1, ((-1., -3.),)))
        self.assertEqual(atlas.call_args.args, (.5, .1, (-1., -3.), 64))
        self.assertEqual(screen.source_witness, witness['source_exact_witness'])
        self.assertIsNot(screen.source_witness, witness['source_exact_witness'])
        self.assertEqual(screen.witness()['outer_measure'], TOEPLITZ_OUTER_MEASURE)

    def test_witness_cli_roundtrip_is_directly_replayable_without_heavy_work(self):
        from hd_shifted_potential_band import replay_saved_band
        source = dict(z='1/2', alpha='1/10', potential='2/3', markers=['1/2 1/8'])
        witness = dict(construction=CONSTRUCTION, pair_swaps=True,
            orientation_aggregation='independent-fair-swaps', grid=64,
            outer_measure=TOEPLITZ_OUTER_MEASURE, source_exact_witness=source)
        module = Path(__file__).with_name('hd_shifted_potential_band_screen.py')
        with tempfile.TemporaryDirectory() as directory:
            input_path, output_path = Path(directory)/'witness.json', Path(directory)/'proposal.json'
            input_path.write_text(json.dumps(witness), encoding='utf-8')
            argv = [str(module), '--witness', str(input_path), '--first', '10', '--last', '12',
                    '--output', str(output_path)]
            with patch('sys.argv', argv), redirect_stdout(StringIO()), \
                    patch('hd_shifted_projected_screen.alias_table', side_effect=self.tables), \
                    patch('hd_shifted_run_potential.run_log_screen', return_value=-100.):
                runpy.run_path(str(module), run_name='__main__')
            proposal = json.loads(output_path.read_text(encoding='utf-8'))
        self.assertEqual(proposal['fixed_witness']['source_exact_witness'], source)
        self.assertEqual(proposal['fixed_witness']['outer_measure'], TOEPLITZ_OUTER_MEASURE)
        self.assertEqual([row['active_groups'] for row in proposal['results']], [10, 11, 12])
        with patch('hd_shifted_potential_band.band_receipt', return_value='recomputed') as replay:
            self.assertEqual(replay_saved_band(proposal), 'recomputed')
        self.assertEqual(replay.call_args.args,
                         (10, 12, F(1, 2), F(1, 10), ((F(1, 2), F(1, 8)),), (F(1), F(2, 3))))
        self.assertEqual(replay.call_args.kwargs['outer_measure'], TOEPLITZ_OUTER_MEASURE)

    def test_one_central_fit_retains_tables_and_exactifies_markers(self):
        def fit(weight, z, alpha, centers, grid):
            full, = centers
            return {}, [dict(odds=(float(full), 2.*full), objective_evaluations=3)]
        with patch('hd_shifted_swap_potential.atlas_bank_screen', side_effect=fit) as fitted, \
                patch('hd_shifted_potential_band_screen.alias_table', side_effect=self.tables) as atlas:
            screen = FixedAtlasScreen.fit_central_once(10, .5, .1, (0, 1), -.5, grid=64)
            screen.regional(9)
            screen.regional(11)
        self.assertEqual(fitted.call_count, 2)
        self.assertEqual(atlas.call_count, 2)
        self.assertEqual(screen.central_marker_fits, 2)
        self.assertEqual(screen.source_witness['markers'], ['1 1', '2 4'])
        self.assertEqual(screen.source_witness['central_objective_evaluations'], 6)

    def test_invalid_fixed_witnesses_and_weight_caps(self):
        for arguments in ((0, .1, ((0, 0),), 0.), (.5, 2., ((0, 0),), 0.),
                          (.5, .1, (), 0.), (.5, .1, ((0, float('inf')),), 0.),
                          (.5, .1, ((0, 0),), float('nan'))):
            with self.assertRaises(ValueError):
                FixedAtlasScreen(*arguments)
        with patch('hd_shifted_potential_band_screen.alias_table', side_effect=self.tables):
            screen = FixedAtlasScreen(.5, .1, ((0, 0),), 0., grid=64)
        for first, last, cap in ((1, 34, 33), (3, 2, 33), (0, 2, 33), (1, 2, 0)):
            with self.assertRaises(ValueError):
                screen.band(first, last, max_weights=cap)
        for weight in (0, 8193, 1.5):
            with self.assertRaises(ValueError):
                screen.evaluate(weight)

    def test_optional_band_fixture_keeps_failures_and_exact_source(self):
        path = Path(__file__).with_name('hd_shifted_potential_band_1024.json')
        if not path.exists():
            self.skipTest('Optional FLOAT artifact absent; no FFT or search runs during unit tests')
        with path.open(encoding='utf-8') as source:
            result = json.load(source)
        self.assertEqual([row['active_groups'] for row in result['results']], list(range(1008, 1041)))
        self.assertEqual(result['failed_weights'], [1039, 1040])
        self.assertFalse(result['all_weights_pass_threshold'])
        self.assertFalse(result['union_passes_threshold'])
        self.assertEqual(result['fft_tables_built'], 22)
        self.assertEqual(result['witness_searches'], 0)
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        values = [row['log2_upper_screen'] for row in result['results']]
        self.assertAlmostEqual(float(np.logaddexp2.reduce(values)), result['union_log2_upper_screen'])
        prefix = [row['log2_upper_screen'] for row in result['results'] if row['active_groups'] <= 1038]
        self.assertLess(float(np.logaddexp2.reduce(prefix)), -40.)
        self.assertEqual(result['fixed_witness']['source_exact_witness']['active_groups'], 1024)

    def test_optional_frontier_band_is_complete_with_a_fixed_exact_witness(self):
        path = Path(__file__).with_name('hd_shifted_potential_band_1024_frontier.json')
        if not path.exists():
            self.skipTest('Optional stronger FLOAT artifact absent; no fitting runs during unit tests')
        with path.open(encoding='utf-8') as source:
            result = json.load(source)
        self.assertEqual([row['active_groups'] for row in result['results']], list(range(1008, 1041)))
        self.assertEqual(result['failed_weights'], [])
        self.assertTrue(result['all_weights_pass_threshold'])
        self.assertTrue(result['union_passes_threshold'])
        self.assertLess(result['union_log2_upper_screen'], -100.)
        self.assertEqual(result['fft_tables_built'], 22)
        self.assertEqual(result['central_marker_fits'], 22)
        self.assertEqual(result['neighbor_witness_searches'], 0)
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        witness = result['fixed_witness']
        exact = witness['source_exact_witness']
        self.assertEqual(F(exact['z']), F(witness['mark']))
        self.assertEqual(F(exact['alpha']), F(witness['alpha']))
        self.assertEqual(F(exact['potential']), F(2**witness['log2_potential_ratio']))
        for row, odds in zip(exact['markers'], witness['marker_log_odds']):
            self.assertEqual(tuple(map(F, row.split())), tuple(F(2**value) for value in odds))


if __name__ == '__main__':
    unittest.main()

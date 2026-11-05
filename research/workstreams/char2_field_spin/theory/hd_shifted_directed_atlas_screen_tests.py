from copy import deepcopy
from contextlib import redirect_stderr
from fractions import Fraction as F
from io import StringIO
import json
from math import comb
import runpy
import unittest
from unittest.mock import patch
import warnings

import numpy as np
from flint import arb, ctx

from hd_shifted_coverage import RegionalAliases
from hd_shifted_directed_atlas_screen import (DirectedAtlasScreen, atlas_log_table,
                                             checked_window, TABLE_SOURCE, COUNT_METHOD)
from hd_shifted_good_pair_certificate import component_power
from hd_shifted_outer_density import DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE
from hd_shifted_potential_band import replay_saved_band
from hd_shifted_potential_band_screen import FixedAtlasScreen
from hd_shifted_run_count_screen import run_log_screen_float
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION
from hd_sumrank_paired_certificate import scalar_power


MODULE = 'hd_shifted_directed_atlas_screen.'


def witness():
    return dict(construction=SWAP_CONSTRUCTION, pair_swaps=True,
        orientation_aggregation='independent-fair-swaps', grid=64,
        # These convenience FLOAT values deliberately disagree with the source.
        mark=.9, alpha=.8, log2_potential_ratio=4,
        outer_measure=TOEPLITZ_OUTER_MEASURE,
        source_exact_witness=dict(z='2/4', alpha='2/20', potential='6/9',
            markers=['1/1 1', '2 4/1', '2/2 1'], note={'retained': [1, 2]}))


class DirectedAtlasScreenTests(unittest.TestCase):
    def test_cli_weight_budget_rejected_before_source_read_or_atlas_build(self):
        self.assertIsNone(checked_window(1150, 1450, max_weights=301))
        for args in ((0, 1), (2, 1), (True, 1), (1, 34)):
            with self.assertRaises(ValueError):
                checked_window(*args)
        with patch('sys.argv', ['hd_shifted_directed_atlas_screen', '--witness', 'unused.json',
                               '--first', '1150', '--last', '1450']), \
                patch('pathlib.Path.open') as source_read, \
                patch('hd_shifted_swap_certificate.swap_unmarked') as packet_build, \
                redirect_stderr(StringIO()), warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            with self.assertRaises(SystemExit) as stopped:
                runpy.run_module('hd_shifted_directed_atlas_screen', run_name='__main__')
            self.assertEqual(stopped.exception.code, 2)
            source_read.assert_not_called()
            packet_build.assert_not_called()

    def build_mock_screen(self, source=None):
        source = witness() if source is None else source
        observed = []
        def power(*pieces, dimensions):
            observed.append((pieces, dimensions, ctx.prec))
            return pieces[0], None, None
        def atlas(row, vector, rho):
            value = 1 if row == (F(1), F(1)) else 5
            return RegionalAliases(np.full((4, 64, 64), value, dtype=float), (arb(1),)*4)
        with patch(MODULE+'swap_unmarked', return_value=({}, {})) as packet, \
                patch(MODULE+'unsigned_components', side_effect=lambda s, p, row: (row,)), \
                patch(MODULE+'component_power', side_effect=power) as transform, \
                patch(MODULE+'RegionalAliases.from_power', side_effect=atlas), \
                patch('hd_shifted_potential_band_screen.alias_table') as floating_fft:
            with ctx.workprec(160):
                screen = DirectedAtlasScreen.from_witness(source, precision=256)
                self.assertEqual(ctx.prec, 160)
            packet.assert_called_once_with(F(1, 2), F(1, 10), 256, 60)
            self.assertEqual(transform.call_count, 2)
            floating_fft.assert_not_called()
        self.assertTrue(all(row[1:] == (2, 256) for row in observed))
        return screen

    def test_exact_source_and_one_transform_per_distinct_marker(self):
        source = witness()
        exact = deepcopy(source['source_exact_witness'])
        screen = self.build_mock_screen(source)
        self.assertEqual(screen.z, .5)
        self.assertEqual(screen.alpha, .1)
        self.assertEqual(screen.outer_measure, TOEPLITZ_OUTER_MEASURE)
        source['source_exact_witness']['note']['retained'].append(3)
        saved = screen.witness()
        self.assertEqual(saved['source_exact_witness'], exact)
        saved['source_exact_witness']['note']['retained'].append(4)
        self.assertEqual(screen.witness()['source_exact_witness'], exact)
        self.assertEqual(saved['table_source']['precision'], 256)
        self.assertEqual(saved['table_source']['packet_mantissa_bits'], 60)
        self.assertEqual(saved['table_source']['distinct_markers'], 2)
        self.assertTrue(saved['requires_fresh_directed_replay'])

    def test_histogram_normalization_precedes_minimum_and_orientation_survives(self):
        screen = self.build_mock_screen()
        matrix = screen.regional(10)[3]
        expected = np.log2(5.)-np.log2(float(comb(4096, 4)*comb(4092, 3)))-10
        np.testing.assert_allclose(matrix, np.full((2, 2), expected), rtol=0, atol=2.e-10)
        atlas = RegionalAliases(np.ones((4, 2, 2)), tuple(arb(2)**v for v in (1, 3, 5, 7)))
        table = atlas_log_table(atlas)
        np.testing.assert_array_equal(table[:, 1, 0].reshape(2, 2), [[1, 3], [5, 7]])

    def test_factor_log_avoids_float_overflow_and_underflow_without_floor(self):
        with ctx.workprec(256):
            atlas = RegionalAliases(np.ones((4, 2, 2)),
                (arb(2)**10000, arb(2)**-10000, arb(2)**20000, arb(2)**-20000))
            atlas.upper[0, 1, 0] = 0
            atlas.upper[1, 0, 1] = np.nextafter(0., 1.)
            table = atlas_log_table(atlas, precision=256)
        self.assertEqual(table[0, 0, 0], 10000)
        self.assertEqual(table[1, 0, 0], -10000)
        self.assertEqual(table[2, 0, 0], 20000)
        self.assertEqual(table[3, 0, 0], -20000)
        self.assertEqual(table[0, 1, 0], -np.inf)
        self.assertEqual(table[1, 0, 1], -11074)
        self.assertFalse(table.flags.writeable)

    def test_invalid_atlas_entries_and_uncertain_factors_rejected(self):
        for value in (-1., np.nan, np.inf, -np.inf):
            upper = np.ones((4, 2, 2))
            upper[1, 0, 0] = value
            with self.assertRaises(ValueError):
                atlas_log_table(RegionalAliases(upper, (arb(1),)*4))
        for factor in (arb(0), arb(-1), arb('nan'), arb('inf'), arb('0 +/- 1')):
            with self.assertRaises(ArithmeticError):
                atlas_log_table(RegionalAliases(np.ones((4, 2, 2)), (factor,)+(arb(1),)*3))

    def test_negative_infinity_allowed_but_nan_and_positive_infinity_rejected(self):
        table = np.zeros((4, 64, 64))
        table[1] = -np.inf
        screen = FixedAtlasScreen(.5, .1, ((0., 0.),), 0., grid=64,
                                  _tables={(0., 0.): table})
        self.assertEqual(screen.regional(2)[0][0, 1], -np.inf)
        self.assertTrue(np.isfinite(screen.evaluate(2)['log2_upper_screen']))
        table[:] = -np.inf
        self.assertEqual(screen.evaluate(2)['log2_upper_screen'], -np.inf)
        for value in (np.nan, np.inf):
            table[0, 0, 0] = value
            with self.assertRaises(ArithmeticError):
                FixedAtlasScreen(.5, .1, ((0., 0.),), 0., grid=64, _tables={(0., 0.): table})

    def test_invalid_identity_mode_and_source_rejected_before_build(self):
        variants = []
        for key, value in (('construction', 'different'), ('pair_swaps', False),
                           ('orientation_aggregation', 'fixed'), ('outer_measure', 'unproved')):
            source = witness()
            source[key] = value
            variants.append(source)
        source = witness()
        source['source_exact_witness']['z'] = .5
        variants.append(source)
        source = witness()
        source['source_exact_witness']['markers'] = ['1']
        variants.append(source)
        with patch(MODULE+'swap_unmarked') as build:
            for source in variants:
                with self.assertRaises(ValueError):
                    DirectedAtlasScreen.from_witness(source)
            # An override must not conceal an invalid saved mode.
            with self.assertRaises(ValueError):
                DirectedAtlasScreen.from_witness(variants[3], outer_measure=DEFAULT_OUTER_MEASURE)
            for kwargs in ({'precision': True}, {'precision': 64}, {'bits': 15}, {'grid': 63}):
                with self.assertRaises(ValueError):
                    DirectedAtlasScreen.from_witness(witness(), **kwargs)
            build.assert_not_called()

    def test_distinct_exact_markers_cannot_silently_collide(self):
        source = witness()
        source['source_exact_witness']['markers'] = ['2 1', str(F(2)+F(1, 2**300))+' 1']
        with patch(MODULE+'swap_unmarked') as build:
            with self.assertRaisesRegex(ValueError, 'collide'):
                DirectedAtlasScreen.from_witness(source)
            build.assert_not_called()

    def test_band_scope_json_and_fresh_replay_exact_literals(self):
        screen = self.build_mock_screen()
        proposal = screen.band(2, 4, max_weights=3)
        self.assertEqual([row['active_groups'] for row in proposal['results']], [2, 3, 4])
        self.assertTrue(proposal['every_weight_evaluated'])
        self.assertFalse(proposal['directed'])
        self.assertFalse(proposal['whole_code_certificate'])
        self.assertTrue(proposal['requires_fresh_directed_replay'])
        self.assertEqual(proposal['table_source']['kind'], TABLE_SOURCE)
        self.assertEqual(proposal['count_method'], COUNT_METHOD)
        self.assertEqual(proposal['fft_tables_built'], 2)
        self.assertEqual(proposal['witness_searches'], 0)
        saved = json.loads(json.dumps(proposal, allow_nan=False))
        with patch('hd_shifted_potential_band.band_receipt', return_value='fresh') as replay:
            self.assertEqual(replay_saved_band(saved, precision=384), 'fresh')
        self.assertEqual(replay.call_args.args[:4], (2, 4, F(1, 2), F(1, 10)))
        self.assertEqual(replay.call_args.args[-1], (F(1), F(2, 3)))
        self.assertEqual(replay.call_args.kwargs['outer_measure'], TOEPLITZ_OUTER_MEASURE)
        self.assertEqual(replay.call_args.kwargs['precision'], 384)

    def test_float_count_backend_matches_tiny_exact_counts_without_mutating_base(self):
        screen = self.build_mock_screen()
        reference = [FixedAtlasScreen.evaluate(screen, weight) for weight in range(1, 8)]
        with patch(MODULE+'run_log_screen_float', wraps=run_log_screen_float) as counted, \
                patch('hd_shifted_potential_band_screen.run_log_screen',
                      side_effect=AssertionError('Bigint count backend must not be called')):
            result = [screen.evaluate(weight) for weight in range(1, 8)]
        self.assertEqual(counted.call_count, 7)
        for actual, expected in zip(result, reference):
            self.assertAlmostEqual(actual['log2_upper_screen'], expected['log2_upper_screen'], delta=1.e-8)
            self.assertEqual(actual['count_method'], COUNT_METHOD)
        # The base remains usable after the override; no global backend patch persists.
        self.assertEqual(FixedAtlasScreen.evaluate(screen, 1), reference[0])

    def test_all_zero_operator_propagates_and_nonfinite_json_fails_closed(self):
        screen = self.build_mock_screen()
        screen.tables = tuple(np.full(table.shape, -np.inf) for table in screen.tables)
        result = screen.band(1, 3, max_weights=3)
        self.assertTrue(all(row['log2_upper_screen'] == -np.inf for row in result['results']))
        self.assertEqual(result['union_log2_upper_screen'], -np.inf)
        self.assertFalse(result['directed'])
        self.assertTrue(result['requires_fresh_directed_replay'])
        with self.assertRaises(ValueError):
            json.dumps(result, allow_nan=False)

    def test_default_mode_explicitly_retained(self):
        source = witness()
        del source['outer_measure']
        screen = self.build_mock_screen(source)
        self.assertEqual(screen.witness()['outer_measure'], DEFAULT_OUTER_MEASURE)
        self.assertEqual(screen.band(1, 1)['outer_measure'], DEFAULT_OUTER_MEASURE)

    def test_constant_matrix_actual_directed_transform_and_residue_queries(self):
        # This small constant polynomial uses real directed FFT machinery,
        # without constructing the expensive production packet coefficients.
        with ctx.workprec(256):
            matrix = tuple(arb(value)/16 for value in (1, 2, 3, 4))
            powered, vector, rho = component_power({}, {(0, 0): matrix}, {}, 64, dimensions=2)
            atlas = RegionalAliases.from_power(powered, vector, rho)
            table = atlas_log_table(atlas, precision=256)
            exact = scalar_power(matrix, 128)
            numerators = atlas.numerator(((0, 0), (0, 0)))
            for actual, expected in zip(numerators, exact):
                self.assertGreaterEqual(actual.upper(), expected.upper())
            for s, f in ((0, 0), (1, 2), (63, 63)):
                for entry, value in enumerate(atlas.numerator(((s, s), (f, f)))):
                    expected = float((value.log()/arb(2).log()).upper()) if not value.is_zero() else -np.inf
                    self.assertAlmostEqual(table[entry, s, f], expected, delta=5.e-12)
            for i in range(2):
                for j in range(2):
                    expected = rho**128*vector[i]/vector[j]
                    self.assertEqual(atlas.factors[2*i+j].lower(), expected.lower())
                    self.assertEqual(atlas.factors[2*i+j].upper(), expected.upper())


if __name__ == '__main__':
    unittest.main()

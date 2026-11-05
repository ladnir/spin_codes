from fractions import Fraction as F
import json
from math import comb
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, ctx

from hd_shifted_potential_band import (band_receipt, directed_band_regions,
    replay_saved_band, checked_inputs, singleton_numerators)
from hd_shifted_coverage import RegionalAliases
from hd_shifted_outer_density import (DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE,
                                     log_density_float)
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION
from hd_shifted_swap_potential import run_bank_receipt


class DirectedBandTests(unittest.TestCase):
    def test_atlas_built_once_per_marker_and_normalization_precedes_min(self):
        def atlas(powered, vector, rho):
            amount = 1 if powered == (F(1), F(1)) else 5
            return RegionalAliases(np.full((4, 128, 128), amount, dtype=float), (arb(1),)*4)
        def power(*pieces, dimensions):
            return pieces[0], None, None
        with patch('hd_shifted_potential_band.swap_unmarked', return_value=({}, {})), \
                patch('hd_shifted_potential_band.unsigned_components', side_effect=lambda s, p, row: (row,)), \
                patch('hd_shifted_potential_band.component_power', side_effect=power) as build, \
                patch('hd_shifted_potential_band.RegionalAliases.from_power', side_effect=atlas):
            with ctx.workprec(160):
                regions = directed_band_regions(10, 12, F(1, 2), F(1, 10),
                    ((F(1), F(1)), (F(2), F(4)), (F(1), F(1))), precision=256)
                self.assertEqual(ctx.prec, 160)
        self.assertEqual(build.call_count, 2)
        self.assertEqual(set(regions), {10, 11, 12})
        with ctx.workprec(256):
            mult = comb(4096, 3)*comb(4093, 4)
            expected = arb(5)/(1024*mult)
            for entry in regions[10][3]:
                self.assertGreaterEqual(entry, expected)
                self.assertLess(entry, expected*arb('1.00000000000000001'))

    def test_singleton_cache_matches_general_query_for_every_small_residue(self):
        with ctx.workprec(256):
            for grid in (2, 4, 8):
                upper = np.arange(4*grid*grid, dtype=float).reshape(4, grid, grid)/17.
                atlas = RegionalAliases(upper, tuple(arb(index+1)/7 for index in range(4)))
                residues = {(s, f) for s in range(grid) for f in range(grid)}
                cached = singleton_numerators(atlas, residues)
                self.assertEqual(len(cached), grid*grid)
                for s, f in residues:
                    expected = atlas.numerator(((s+3*grid, s+3*grid), (f+grid, f+grid)))
                    for actual, reference in zip(cached[s, f], expected):
                        self.assertEqual(actual.lower(), reference.lower())
                        self.assertEqual(actual.upper(), reference.upper())

    def test_equal_residues_keep_distinct_histogram_normalization(self):
        atlas = RegionalAliases(np.ones((4, 64, 64)), (arb(1),)*4)
        with patch('hd_shifted_potential_band.swap_unmarked', return_value=({}, {})), \
                patch('hd_shifted_potential_band.unsigned_components', return_value=()), \
                patch('hd_shifted_potential_band.component_power', return_value=(None, None, None)), \
                patch('hd_shifted_potential_band.RegionalAliases.from_power', return_value=atlas), \
                patch('hd_shifted_potential_band.singleton_numerators', wraps=singleton_numerators) as cached:
            regions = directed_band_regions(2, 66, F(1, 2), F(1, 10), ((F(2), F(3)),),
                                             grid=64, max_weights=65)
        self.assertEqual(cached.call_count, 1)
        # (s,f)=(2,0) and (66,0) use the same residue but different monomials/multinomials.
        with ctx.workprec(256):
            for weight in (2, 66):
                expected = arb(1)/(arb(comb(4096, weight))*arb(2)**weight)
                self.assertGreaterEqual(regions[weight][0][0], expected)
                self.assertLess(regions[weight][0][0], expected*arb('1.00000000000000001'))
            self.assertGreater(regions[2][0][0], regions[66][0][0])

    def test_one_weight_matches_existing_run_receipt(self):
        weight = 4
        regions = {f: tuple(map(arb, (1+f, 2, 3, 4+f))) for f in range(3)}
        z, alpha, markers, vector = F(1, 2), F(1, 10), ((F(1), F(1)),), (F(2), F(3))
        with patch('hd_shifted_potential_band.directed_band_regions', return_value={weight: regions}):
            band = band_receipt(weight, weight, z, alpha, markers, vector, precision=256)
        with patch('hd_shifted_swap_potential.directed_bank_regions', return_value=regions):
            single = run_bank_receipt(weight, z, alpha, markers, vector, precision=256)
        self.assertEqual(band['log2_upper_rational'], single['log2_upper_rational'])
        self.assertEqual(band['results'][0]['log2_upper_rational'], single['log2_upper_rational'])
        self.assertFalse(band['whole_code_certificate'])
        self.assertNotIn('outer_measure', band)
        self.assertNotIn('outer_measure', band['witness'])

    def test_sharper_outer_measure_changes_only_prefix_and_records_identity(self):
        args = (16, 16, F(1, 2), F(1, 10), ((1, 1),), (1, 1))
        with patch('hd_shifted_potential_band.directed_band_regions', return_value={16: {}}), \
                patch('hd_shifted_potential_band.run_sum_arb', return_value=arb(1)):
            legacy = band_receipt(*args)
            sharp = band_receipt(*args, outer_measure=TOEPLITZ_OUTER_MEASURE)
        difference = float(F(legacy['log2_upper_rational'])-F(sharp['log2_upper_rational']))
        expected = 1.6*(log_density_float()-log_density_float(TOEPLITZ_OUTER_MEASURE))
        self.assertAlmostEqual(difference, expected, delta=2.e-12)
        self.assertEqual(sharp['outer_measure'], TOEPLITZ_OUTER_MEASURE)
        self.assertEqual(sharp['witness']['outer_measure'], TOEPLITZ_OUTER_MEASURE)
        self.assertEqual(sharp['outer_measure_bound']['density_formula'], '256^160/255^32')
        self.assertEqual(sharp['construction'], legacy['construction'])
        with patch('hd_shifted_potential_band.directed_band_regions') as build:
            with self.assertRaises(ValueError):
                band_receipt(*args, outer_measure='caller-chosen-density')
            build.assert_not_called()

    def test_every_weight_replayed_once_and_summed(self):
        regions = {r: {f: (arb(1), arb(0), arb(0), arb(1))
                       for f in range(r//2+1)} for r in range(2, 5)}
        observed = []
        def support_sum(pairs, weight, matrices, vector):
            observed.append((weight, tuple(vector)))
            self.assertGreaterEqual(ctx.prec, 256)
            return arb(weight)
        with patch('hd_shifted_potential_band.directed_band_regions', return_value=regions), \
                patch('hd_shifted_potential_band.run_sum_arb', side_effect=support_sum):
            result = band_receipt(2, 4, F(1, 2), F(1, 10000), ((1, 1),), (1, 2), precision=256)
        self.assertEqual(observed, [(r, (F(1), F(2))) for r in range(2, 5)])
        self.assertEqual(result['weights_replayed'], 3)
        self.assertEqual(result['weight_range'], [2, 4])
        self.assertTrue(result['every_integer_weight_replayed'])
        self.assertGreater(F(result['log2_upper_rational']),
                           max(F(row['log2_upper_rational']) for row in result['results']))

    def test_exact_input_and_hard_weight_cap(self):
        for first, last, z, markers in ((1, 34, F(1, 2), ((1, 1),)),
                (0, 1, F(1, 2), ((1, 1),)), (1, 2, .5, ((1, 1),)),
                (1, 2, F(1, 2), ((1., 1),)), (1, 2, F(1, 2), ((0, 1),))):
            with self.assertRaises(ValueError):
                checked_inputs(first, last, z, F(1, 10), markers, 128, 256, 60)
        with self.assertRaises(ValueError):
            replay_saved_band(dict(construction='wrong'))

    def test_explicit_larger_budget_enumerates_every_weight(self):
        z, alpha, markers = F(1, 2), F(1, 10000), ((1, 1),)
        self.assertEqual(checked_inputs(1, 34, z, alpha, markers, 128, 256, 60, max_weights=34)[:2],
                         (z, alpha))
        for cap in (0, 8193, 34., True):
            with self.assertRaises(ValueError):
                checked_inputs(1, 34, z, alpha, markers, 128, 256, 60, max_weights=cap)
        regions = {r: {f: (arb(1), arb(0), arb(0), arb(1))
                       for f in range(r//2+1)} for r in range(1, 35)}
        seen = []
        def support_sum(pairs, weight, matrices, vector):
            seen.append(weight)
            return arb(1)
        with patch('hd_shifted_potential_band.directed_band_regions', return_value=regions) as normalize, \
                patch('hd_shifted_potential_band.run_sum_arb', side_effect=support_sum):
            result = band_receipt(1, 34, z, alpha, markers, (1, 1), max_weights=34)
        self.assertEqual(seen, list(range(1, 35)))
        self.assertEqual(result['weights_replayed'], 34)
        self.assertEqual(result['max_weights'], 34)
        self.assertEqual(normalize.call_args.kwargs['max_weights'], 34)

    def test_saved_replay_ignores_float_bounds_and_uses_exact_literals(self):
        proposal = dict(construction=SWAP_CONSTRUCTION, requested_weight_range=[10, 12],
            union_log2_upper_screen=-1.e300,
            fixed_witness=dict(construction=SWAP_CONSTRUCTION, pair_swaps=True,
                orientation_aggregation='independent-fair-swaps', grid=128,
                mark=.9, alpha=.8, source_exact_witness=dict(
                    z='1/2', alpha='1/10', potential='2/3', markers=['3/4 4/5'])))
        with patch('hd_shifted_potential_band.band_receipt', return_value='recomputed') as replay:
            self.assertEqual(replay_saved_band(proposal, precision=384), 'recomputed')
        self.assertEqual(replay.call_args.args,
                         (10, 12, F(1, 2), F(1, 10), ((F(3, 4), F(4, 5)),), (F(1), F(2, 3))))
        self.assertEqual(replay.call_args.kwargs['precision'], 384)
        self.assertEqual(replay.call_args.kwargs['max_weights'], 33)
        self.assertEqual(replay.call_args.kwargs['outer_measure'], DEFAULT_OUTER_MEASURE)
        proposal['fixed_witness']['outer_measure'] = TOEPLITZ_OUTER_MEASURE
        with patch('hd_shifted_potential_band.band_receipt') as replay:
            replay_saved_band(proposal)
            self.assertEqual(replay.call_args.kwargs['outer_measure'], TOEPLITZ_OUTER_MEASURE)
            replay_saved_band(proposal, outer_measure=DEFAULT_OUTER_MEASURE)
            self.assertEqual(replay.call_args.kwargs['outer_measure'], DEFAULT_OUTER_MEASURE)
        proposal['fixed_witness']['outer_measure'] = 'unproved'
        with self.assertRaises(ValueError):
            replay_saved_band(proposal)
        del proposal['fixed_witness']['outer_measure']
        # A FLOAT proposal cannot increase the caller's resource budget.
        proposal.update(requested_weight_range=[1, 34], max_weights=8192)
        with self.assertRaises(ValueError):
            replay_saved_band(proposal)

    def test_optional_precision_receipts_agree_for_every_weight(self):
        paths = [Path(__file__).with_name('hd_shifted_potential_band_1024_directed%d.json' % precision)
                 for precision in (256, 384)]
        if not all(path.exists() for path in paths):
            self.skipTest('Optional directed receipts absent; unit tests never launch a proof replay')
        receipts = [json.loads(path.read_text(encoding='utf-8')) for path in paths]
        for receipt in receipts:
            self.assertEqual([row['active_groups'] for row in receipt['results']], list(range(1008, 1041)))
            self.assertEqual(receipt['weights_replayed'], 33)
            self.assertEqual(receipt['atlas_builds'], 22)
            self.assertEqual(receipt['integer_profile_triples'], 8677104)
            self.assertTrue(receipt['directed'])
            self.assertFalse(receipt['whole_code_certificate'])
            self.assertEqual(receipt['construction'], SWAP_CONSTRUCTION)
            self.assertLess(F(receipt['log2_upper_rational']), -100)
        self.assertEqual(receipts[0]['log2_upper_rational'], receipts[1]['log2_upper_rational'])
        self.assertEqual([row['log2_upper_rational'] for row in receipts[0]['results']],
                         [row['log2_upper_rational'] for row in receipts[1]['results']])


if __name__ == '__main__':
    unittest.main()

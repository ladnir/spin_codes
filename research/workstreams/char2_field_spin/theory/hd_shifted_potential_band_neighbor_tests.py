from fractions import Fraction as F
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from hd_shifted_potential_band_neighbor import (checked_window, fixed_screen_from_artifact,
    fit_one_central_bank, artifact_outer_measure)
from hd_shifted_outer_density import DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE
from hd_shifted_swap_screen import CONSTRUCTION


class NeighborTests(unittest.TestCase):
    def test_source_bound_is_preserved_unless_explicitly_overridden(self):
        self.assertEqual(artifact_outer_measure({}), DEFAULT_OUTER_MEASURE)
        source = dict(fixed_witness=dict(outer_measure=TOEPLITZ_OUTER_MEASURE))
        self.assertEqual(artifact_outer_measure(source), TOEPLITZ_OUTER_MEASURE)
        self.assertEqual(artifact_outer_measure(source, DEFAULT_OUTER_MEASURE), DEFAULT_OUTER_MEASURE)
        with self.assertRaises(ValueError):
            artifact_outer_measure(source, 'unproved')
        with self.assertRaises(ValueError):
            artifact_outer_measure(dict(fixed_witness=dict(outer_measure='unproved')))

    def test_explicit_window_budget_preserves_default(self):
        self.assertEqual(checked_window(1, 33), (1, 33, 33))
        self.assertEqual(checked_window(1, 34, 34), (1, 34, 34))
        self.assertEqual(checked_window(567, 974, 408), (567, 974, 408))
        for args in ((1, 34), (1, 1, 0), (1, 1, 8193), (1, 1, True),
                     (0, 33, 34), (8192, 8193, 33), (2, 1, 33)):
            with self.assertRaises(ValueError):
                checked_window(*args)

    def test_exact_literals_not_float_bound_are_reused(self):
        source = dict(z='1/2', alpha='1/10', potential='1/4', markers=['1/2 1/8'])
        artifact = dict(construction=CONSTRUCTION, union_log2_upper_screen=-1.e300,
            fixed_witness=dict(construction=CONSTRUCTION, pair_swaps=True,
                orientation_aggregation='independent-fair-swaps', grid=128,
                mark=.9, alpha=.8, source_exact_witness=source))
        target = SimpleNamespace()
        with patch('hd_shifted_potential_band_neighbor.FixedAtlasScreen', return_value=target) as build:
            self.assertIs(fixed_screen_from_artifact(artifact), target)
        self.assertEqual(build.call_args.args, (F(1, 2), F(1, 10), ((-1., -3.),), -2.))
        self.assertEqual(target.source_witness, source)
        self.assertIsNot(target.source_witness, source)

    def test_construction_and_exact_source_required(self):
        for artifact in ({}, dict(construction='wrong'),
                dict(construction=CONSTRUCTION, fixed_witness=dict(construction=CONSTRUCTION,
                    pair_swaps=True, orientation_aggregation='independent-fair-swaps'))):
            with self.assertRaises(ValueError):
                fixed_screen_from_artifact(artifact)

    def test_only_one_central_fit_keeps_alpha_and_potential(self):
        source = dict(z='1/2', alpha='1/10', potential='1/4', markers=['1/2 1/8'], active_groups=1024)
        artifact = dict(construction=CONSTRUCTION,
            fixed_witness=dict(construction=CONSTRUCTION, pair_swaps=True,
                orientation_aggregation='independent-fair-swaps', grid=128, source_exact_witness=source))
        target = SimpleNamespace(source_witness={})
        with patch('hd_shifted_potential_band_neighbor.FixedAtlasScreen.fit_central_once',
                   return_value=target) as fitted:
            self.assertIs(fit_one_central_bank(artifact, 992), target)
        self.assertEqual(fitted.call_count, 1)
        weight, z, alpha, centers, ratio = fitted.call_args.args
        self.assertEqual((weight, alpha, ratio), (992, F(1, 10), -2.))
        self.assertEqual(z, .5**(992/1024))
        self.assertEqual((len(centers), min(centers), max(centers)), (22, 0, 496))

    def test_optional_neighbor_receipts_cover_complete_bands_at_both_precisions(self):
        cases = (('512_496_528', 496, 528, 2176240, F(-738148811533431, 10**12)),
                 ('992_975_1007', 975, 1007, 8127200, F(-143355167233523, 10**12)))
        paths = [[Path(__file__).with_name('hd_shifted_potential_band_%s_directed%d.json' % (name, p))
                  for p in (256, 384)] for name, *_ in cases]
        if not all(path.exists() for pair in paths for path in pair):
            self.skipTest('Optional directed artifacts absent; no proof replay runs during unit tests')
        for (_, first, last, profiles, ceiling), pair in zip(cases, paths):
            receipts = [json.loads(path.read_text(encoding='utf-8')) for path in pair]
            for result in receipts:
                self.assertEqual(result['weight_range'], [first, last])
                self.assertEqual([row['active_groups'] for row in result['results']], list(range(first, last+1)))
                self.assertEqual(result['weights_replayed'], 33)
                self.assertEqual(result['integer_profile_triples'], profiles)
                self.assertEqual(F(result['log2_upper_rational']), ceiling)
                self.assertTrue(result['directed'])
                self.assertFalse(result['whole_code_certificate'])
            self.assertEqual([row['log2_upper_rational'] for row in receipts[0]['results']],
                             [row['log2_upper_rational'] for row in receipts[1]['results']])

    def test_explicit_marker_step_is_bounded_and_does_not_change_default(self):
        source = dict(z='1/2', alpha='1/10', potential='1/4', markers=['1 1'], active_groups=3072)
        artifact = dict(construction=CONSTRUCTION, fixed_witness=dict(construction=CONSTRUCTION,
            pair_swaps=True, orientation_aggregation='independent-fair-swaps',
            grid=128, source_exact_witness=source))
        with patch('hd_shifted_potential_band_neighbor.FixedAtlasScreen.fit_central_once',
                   return_value=SimpleNamespace(source_witness={})) as fit:
            screen = fit_one_central_bank(artifact, 3584, marker_step=32)
            centers = fit.call_args.args[3]
            self.assertEqual((len(centers), min(centers), max(centers)), (62, 0, 1792))
            self.assertEqual(screen.source_witness['central_marker_step'], 32)
            fit.reset_mock()
            for step in (0, -1, True, 1.5, 1):
                with self.assertRaises(ValueError):
                    fit_one_central_bank(artifact, 3584, marker_step=step)
            fit.assert_not_called()

    def test_central_fit_preserves_nonbinary_exact_alpha_and_potential(self):
        source = dict(z='1/2', alpha='1/10', potential='1/3', markers=['1 1'], active_groups=3072)
        artifact = dict(construction=CONSTRUCTION, fixed_witness=dict(construction=CONSTRUCTION,
            pair_swaps=True, orientation_aggregation='independent-fair-swaps',
            grid=128, source_exact_witness=source))
        target = SimpleNamespace(source_witness=dict(alpha=str(F(.1)), potential=str(F(1./3))))
        with patch('hd_shifted_potential_band_neighbor.FixedAtlasScreen.fit_central_once',
                   return_value=target):
            screen = fit_one_central_bank(artifact, 3584, marker_step=32)
        self.assertEqual(screen.source_witness['alpha'], '1/10')
        self.assertEqual(screen.source_witness['potential'], '1/3')

    def test_optional_new_broad_receipts_match_at_both_precisions(self):
        cases = ((320, 280, 349, 1754585, F(-86242628866867, 10**12)),
                 (768, 700, 809, 15744905, F(-16459469287143, 200000000000)))
        paths = [[Path(__file__).with_name('hd_shifted_potential_band_%d_%d_%d_directed%d.json' %
                  (center, first, last, p)) for p in (256, 384)]
                 for center, first, last, *_ in cases]
        if not all(path.exists() for pair in paths for path in pair):
            self.skipTest('Optional broad directed artifacts absent; no replay occurs during tests')
        for (_, first, last, profiles, ceiling), pair in zip(cases, paths):
            receipts = [json.loads(path.read_text(encoding='utf-8')) for path in pair]
            for result in receipts:
                self.assertEqual(result['weight_range'], [first, last])
                self.assertEqual([row['active_groups'] for row in result['results']], list(range(first, last+1)))
                self.assertEqual(result['max_weights'], last-first+1)
                self.assertEqual(result['weights_replayed'], last-first+1)
                self.assertEqual(result['atlas_builds'], 22)
                self.assertEqual(result['integer_profile_triples'], profiles)
                self.assertEqual(F(result['log2_upper_rational']), ceiling)
                self.assertTrue(result['directed'])
                self.assertTrue(result['every_integer_weight_replayed'])
                self.assertFalse(result['whole_code_certificate'])
            self.assertEqual(receipts[0]['witness']['z'], receipts[1]['witness']['z'])
            self.assertEqual([row['log2_upper_rational'] for row in receipts[0]['results']],
                             [row['log2_upper_rational'] for row in receipts[1]['results']])


if __name__ == '__main__':
    unittest.main()

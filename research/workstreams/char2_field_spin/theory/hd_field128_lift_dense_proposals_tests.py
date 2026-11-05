"""Tiny geometry, artifact, and mock-atlas tests; no dense numerical gate."""

from fractions import Fraction as F
from contextlib import redirect_stdout
from io import StringIO
import json
from math import log2
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import hd_field128_lift_dense_proposals as dense
import hd_field128_lift_alias_cover as cover


class DenseProposalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = dense.read_source()

    def test_identity_center_preserves_all_original_markers(self):
        markers, records = dense.recentered_markers(self.source, 3936)
        self.assertEqual(markers[:67], self.source['markers'])
        self.assertEqual([r['source_marker_index'] for r in records[:67]], list(range(67)))
        self.assertTrue(all(r['kind'] == 'warped-original' for r in records[:67]))

    def test_dense_warp_retains_whole_source_grid_and_boundary_proxies(self):
        for center, exponent in dense.DEFAULT_CANDIDATES:
            markers, records = dense.recentered_markers(self.source, center)
            self.assertLessEqual(len(markers), 128)
            self.assertGreaterEqual(len(markers), 67)
            originals = records[:67]
            self.assertEqual(F(originals[0]['target_full']), center-dense.PAIRS)
            self.assertEqual(F(originals[-1]['target_full']), F(center, 2))
            self.assertTrue(originals[0]['uses_empty_proxy'])
            self.assertTrue(originals[-1]['uses_zero_marker'])
            self.assertGreater(markers[0][0], 1 << 32)
            self.assertEqual(markers[66][0], dense.ZERO_MARKER)
            empty_counts = {F(row['target_histogram'][0]) for row in records}
            self.assertTrue(set(map(F, (0, 1, 2, 4, 8, 16, 32))) <= empty_counts)
            self.assertTrue(all(min(pair) > 0 for pair in markers))

    def test_short_exact_proposal_accepted_by_frozen_directed_interface(self):
        artifact = dense.make_witness(self.source, 4428, 18)
        witness = cover.witness_from_artifact(artifact)
        self.assertEqual(witness['z'], dense.rounded(self.source['z']**18))
        self.assertEqual(witness['alpha'], F(1, 160))
        self.assertEqual(len(witness['markers']), artifact['marker_count'])
        self.assertEqual(artifact['source']['sha256'], dense.OLD_SHA256)
        self.assertFalse(artifact['directed'])
        self.assertFalse(artifact['whole_code_certificate'])
        self.assertLess(len(str(witness['z'].numerator)), 19)

    def test_complete_mock_profile_diagnostic_and_floor_reporting(self):
        grid = 64
        logs = np.full((4, grid, grid), 2.)
        tables = (SimpleNamespace(logs=logs, reliable=np.ones_like(logs, dtype=bool),
                                  positive_numerator=np.full(4, 3.)),
                  SimpleNamespace(logs=logs-10, reliable=np.zeros_like(logs, dtype=bool),
                                  positive_numerator=np.full(4, 4.)))
        bank = SimpleNamespace(odds=((0., 0.), (0., 0.)), tables=tables, grid=grid,
                               floor_policy='positive-fallback', alpha=.00625, density=16384.,
                               forbidden=734003, z=.125, constant=13762.55625)
        budget = dense.tail.Budget(10.)
        report = dense.diagnose(bank, 8, -11., budget=budget)
        self.assertEqual(report['complete_histogram_count'], 5)
        self.assertEqual(report['selected_below_floor_entries'], 0)
        self.assertAlmostEqual(sum(report['decomposition'].values()), report['log2_largest_summand_screen'], places=8)
        for histogram in report['dominant_cyclic_run_summand']:
            self.assertEqual(histogram['single']+2*histogram['full'], 8)
            self.assertTrue(all(entry['marker_index'] == 0 for entry in histogram['entries']))
        bank.floor_policy = 'legacy-cushion'
        report = dense.diagnose(bank, 8, -11., budget=budget)
        self.assertEqual(report['selected_below_floor_entries'], 20)
        self.assertTrue(all(entry['marker_index'] == 1 for row in report['dominant_cyclic_run_summand']
                            for entry in row['entries']))

    def test_timeout_keeps_requested_weights_pending_and_no_coverage_claim(self):
        witness = dense.make_witness(self.source, 4428, 18)
        with patch.object(dense.alias, 'AliasBankScreen', side_effect=TimeoutError('mock deadline')):
            result = dense.screen_witness(witness, (4428, 4490), max_seconds=1.)
        self.assertEqual(result['pending_weights'], [4428, 4490])
        self.assertEqual(result['results'], [])
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['directed'])

    def test_exact_ratio_ranges_output_pins_and_input_guards(self):
        for log_ratio in (-2048.5, -11.25, 2048.5):
            value = dense.exact_ratio(log_ratio)
            self.assertGreater(value, 0)
            self.assertAlmostEqual(dense.log_rational(value), log_ratio, places=10)
        for center in (0, 8192):
            with self.assertRaises(ValueError):
                dense.recentered_markers(self.source, center)
        with self.assertRaises(ValueError):
            dense.make_witness(self.source, 4428, 0)
        with self.assertRaises(ValueError):
            dense.recentered_markers(self.source, 4428, empty_proxy=0)
        with TemporaryDirectory() as directory:
            target = Path(directory)/'new.json'
            pin = dense.write_new(target, {'test': True})
            self.assertEqual(dense.sha256(target.read_bytes()).hexdigest(), pin['sha256'])
            with self.assertRaises(FileExistsError):
                dense.write_new(target, {})
        self.assertIn(Path(dense.__file__).name, dense.source_pins())

    def test_cli_deadline_after_check_keeps_candidate_pending_without_build(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)/'plan'
            argv = ['dense', '--mode', 'float', '--candidate', '4428', '18',
                    '--output-dir', str(output), '--max-seconds', '1']
            with patch('sys.argv', argv), patch.object(dense, 'monotonic', return_value=1e30), \
                    patch.object(dense, 'screen_witness') as run, redirect_stdout(StringIO()):
                dense.main()
            run.assert_not_called()
            result = json.loads((output/'summary.json').read_text())
            self.assertEqual(result['failure']['error_type'], 'TimeoutError')
            self.assertEqual(result['completed_candidates'], [])
            self.assertEqual(result['pending_candidate_indices'], [0])


if __name__ == '__main__':
    unittest.main()

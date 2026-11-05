from copy import deepcopy
from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction as F
from hashlib import sha256
from io import StringIO
import json
from math import log2
from pathlib import Path
import runpy
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from hd_grouped_packet_screen import LIMIT
from hd_shifted_density_proposal import reweight_proposal
from hd_shifted_outer_density import (DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE,
                                     log_density_float)
from hd_shifted_swap_screen import CONSTRUCTION


class DensityProposalTests(unittest.TestCase):
    def fixture(self):
        return dict(construction=CONSTRUCTION, directed=False, whole_code_certificate=False,
            method='fixed-atlas-common-potential-run-convolution',
            every_weight_evaluated=True, requested_weight_range=[2, 4], weights_evaluated=3,
            max_weights=3, threshold=-80., fft_tables_built=2, central_marker_fits=2,
            results=[dict(active_groups=r, log2_support_sum_screen=-100.,
                          log2_outer_prefix=1.e100, log2_upper_screen=1.e100) for r in range(2, 5)],
            fixed_witness=dict(construction=CONSTRUCTION, pair_swaps=True,
                orientation_aggregation='independent-fair-swaps', mark=.9, alpha=.8,
                source_exact_witness=dict(z='1/2', alpha='1/10000', potential='2/3',
                                           markers=['1/2 1/8'])))

    def test_exact_prefix_only_deepcopy_metadata_and_complete_union(self):
        source = self.fixture()
        original = deepcopy(source)
        result = reweight_proposal(source, TOEPLITZ_OUTER_MEASURE)
        self.assertEqual(source, original)
        self.assertEqual(result['fixed_witness']['source_exact_witness'], source['fixed_witness']['source_exact_witness'])
        values = []
        for row in result['results']:
            prefix = .0001*row['active_groups']*log_density_float(TOEPLITZ_OUTER_MEASURE)-.0001*(log2(255)+LIMIT*log2(.5))
            self.assertEqual(row['log2_outer_prefix'], prefix)
            self.assertEqual(row['log2_support_sum_screen'], -100.)
            self.assertEqual(row['log2_upper_screen'], prefix-100.)
            self.assertFalse(row['directed'])
            values.append(prefix-100.)
        self.assertEqual(result['union_log2_upper_screen'], float(np.logaddexp2.reduce(values)))
        self.assertEqual(result['failed_weights'], [2, 3, 4])
        self.assertFalse(result['all_weights_pass_threshold'])
        self.assertFalse(result['union_passes_threshold'])
        self.assertEqual(result['outer_measure_bound']['density_formula'], '256^160/255^32')
        self.assertEqual(result['fixed_witness']['outer_measure'], TOEPLITZ_OUTER_MEASURE)
        self.assertEqual((result['new_fft_tables_built'], result['new_witness_searches'], result['new_central_marker_fits']), (0, 0, 0))
        self.assertEqual(result['fft_tables_built'], 2)
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])

    def test_legacy_restore_omits_bound_metadata_and_threshold_can_change(self):
        sharp = reweight_proposal(self.fixture(), TOEPLITZ_OUTER_MEASURE)
        legacy = reweight_proposal(sharp, DEFAULT_OUTER_MEASURE, threshold=0.)
        for record in (legacy, legacy['fixed_witness']):
            self.assertNotIn('outer_measure', record)
            self.assertNotIn('outer_measure_bound', record)
        self.assertEqual(legacy['density_reweighted_from'], TOEPLITZ_OUTER_MEASURE)
        self.assertEqual(legacy['failed_weights'], [])
        self.assertTrue(legacy['all_weights_pass_threshold'])
        self.assertTrue(legacy['union_passes_threshold'])

    def test_malformed_construction_exact_witness_or_weight_rows_rejected(self):
        mutators = (lambda d: d.update(construction='other'),
            lambda d: d.update(directed=True),
            lambda d: d.update(every_weight_evaluated=False),
            lambda d: d.update(weights_evaluated=2),
            lambda d: d['results'].pop(),
            lambda d: d['results'][1].update(active_groups=2),
            lambda d: d['results'].reverse(),
            lambda d: d['results'][0].update(log2_support_sum_screen=float('nan')),
            lambda d: d['fixed_witness']['source_exact_witness'].update(z=.5),
            lambda d: d['fixed_witness']['source_exact_witness'].update(alpha='0'),
            lambda d: d['fixed_witness']['source_exact_witness'].update(markers=['0 1']),
            lambda d: d['fixed_witness'].update(outer_measure='unproved'))
        for mutate in mutators:
            source = self.fixture()
            mutate(source)
            with self.assertRaises(ValueError):
                reweight_proposal(source, TOEPLITZ_OUTER_MEASURE)
        with self.assertRaises(ValueError):
            reweight_proposal(self.fixture(), 'unproved')

    def test_optional_real_band_uses_exact_literals_without_fft_or_fit(self):
        path = Path(__file__).with_name('hd_shifted_potential_band_1280_zminus004_1150_1450.json')
        if not path.exists():
            self.skipTest('Optional FLOAT artifact absent; no atlas or fit runs during tests')
        source = json.loads(path.read_text(encoding='utf-8'))
        result = reweight_proposal(source, TOEPLITZ_OUTER_MEASURE)
        self.assertEqual([row['active_groups'] for row in result['results']], list(range(1150, 1451)))
        self.assertEqual(result['weights_evaluated'], 301)
        alpha = float(F(source['fixed_witness']['source_exact_witness']['alpha']))
        delta = log_density_float()-log_density_float(TOEPLITZ_OUTER_MEASURE)
        for old, new in zip(source['results'], result['results']):
            self.assertEqual(old['log2_support_sum_screen'], new['log2_support_sum_screen'])
            self.assertAlmostEqual(old['log2_upper_screen']-new['log2_upper_screen'],
                                   alpha*old['active_groups']*delta, delta=1.e-9)

    def test_cli_exclusive_output_and_source_hash_without_numerical_rebuild(self):
        module = Path(__file__).with_name('hd_shifted_density_proposal.py')
        raw = json.dumps(self.fixture()).encode('utf-8')
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory)/'source.json', Path(directory)/'proposal.json'
            source.write_bytes(raw)
            argv = [str(module), '--source', str(source), '--outer-measure', TOEPLITZ_OUTER_MEASURE,
                    '--output', str(output)]
            with patch('sys.argv', argv), redirect_stdout(StringIO()), \
                    patch('hd_shifted_projected_screen.alias_table', side_effect=AssertionError('No FFT allowed')):
                runpy.run_path(str(module), run_name='__main__')
            saved = output.read_bytes()
            result = json.loads(saved)
            self.assertEqual(result['source_screen_sha256'], sha256(raw).hexdigest())
            self.assertEqual(result['new_fft_tables_built'], 0)
            self.assertFalse(result['directed'])
            with patch('sys.argv', argv), redirect_stderr(StringIO()), self.assertRaises(SystemExit):
                runpy.run_path(str(module), run_name='__main__')
            self.assertEqual(output.read_bytes(), saved)


if __name__ == '__main__':
    unittest.main()

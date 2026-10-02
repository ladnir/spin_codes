"""Scope/authentication regression checks without replaying the large census."""
import copy
from fractions import Fraction as Q
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import dense_model as dm


class DenseScopeTests(unittest.TestCase):
    def test_fingerprint_survives_integer_spectrum_keys(self):
        spectrum = {'weights': {0: 1, 32: 40, 128: 1}}
        self.assertEqual(dm.fingerprint(spectrum), dm.fingerprint(json.loads(json.dumps(spectrum))))

    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.source = Path(self.temp.name)/'mixture.json'
        self.premises = {'canonical_cdf': [1, 2], 'test_identity': 'fresh'}
        self.caps, self.comparison = [3, 4], [5, 6]
        self.saved = dict(K=1 << 20, N=1 << 21, block_width=8,
            minimum_groups=33, maximum_groups=2048,
            comparison='direct-expected-shell-majorant',
            outer_premises=self.premises,
            expected_cdf_sha256=dm.cell_search.dense.fingerprint(self.caps),
            comparison_caps_sha256=dm.cell_search.dense.fingerprint(self.comparison),
            mixture=[{'mass': '1', 'activity': '1/2'}])
        self.write()
        self.auth = self.enter(patch.object(dm, 'authenticated_bch_cdf',
            return_value=(self.caps, self.premises)))
        self.enter(patch.object(dm, 'transport_shells', return_value=self.comparison))
        self.exact = self.enter(patch.object(dm.cell_search.dense, 'exact_mixture',
            return_value=([(Q(1), Q(1, 2))], {'checked_every_shell': True})))
        self.maps = self.enter(patch.object(dm.kernel_maps, 'prepare',
            return_value=({'test': 'fresh-data'}, {'test': 'fresh-maps'})))
        self.enter(patch.object(dm, 'Model',
            return_value=SimpleNamespace(root=(Q(1, 100), Q(1)), tilt=Q(3, 16))))

    def enter(self, context):
        value = context.start()
        self.addCleanup(context.stop)
        return value

    def write(self):
        self.source.write_text(json.dumps({'scope': self.saved, 'forged_upper': [1, -9999]}))

    def test_fresh_counts_and_exact_mixture_before_maps(self):
        _, scope, provenance = dm.fresh_model(self.source)
        self.auth.assert_called_once_with()
        self.exact.assert_called_once_with(self.comparison, self.saved['mixture'])
        self.maps.assert_called_once_with(feedback='bch16', seed=0, refresh='uniform', updates=8)
        self.assertEqual(scope['sampling'], {'kind': 'uniform_gl', 'bits': 16})
        self.assertEqual(scope['maps'], {'test': 'fresh-maps'})
        self.assertEqual(scope['threshold'], 209715)
        self.assertNotIn('forged_upper', scope)
        self.assertEqual(provenance['path'], str(self.source.resolve()))

    def test_altered_premises_rejected_before_inner(self):
        self.saved['outer_premises'] = {'canonical_cdf': [1, 3]}
        self.write()
        with self.assertRaisesRegex(ValueError, 'fresh outer count'):
            dm.fresh_model(self.source)
        self.maps.assert_not_called()

    def test_invalid_shell_mixture_rejected_before_inner(self):
        self.exact.side_effect = ValueError('shell majorant failed')
        with self.assertRaisesRegex(ValueError, 'shell majorant'):
            dm.fresh_model(self.source)
        self.maps.assert_not_called()

    def test_source_change_during_preparation_rejected(self):
        def changed(**kwargs):
            self.source.write_text('changed source')
            return {}, {'test': 'fresh-maps'}
        self.maps.side_effect = changed
        with self.assertRaisesRegex(ValueError, 'source changed'):
            dm.fresh_model(self.source)

    def test_transvection_scope_stays_distinct(self):
        _, scope, _ = dm.fresh_model(self.source, refresh='transvections', updates=16)
        self.assertEqual(scope['sampling'], {'kind': 'transvections', 'bits': 16, 'updates': 16})

    def test_wrong_geometry_rejected_before_authentication(self):
        self.saved['N'] = 1 << 20
        self.write()
        with self.assertRaisesRegex(ValueError, 'matching K20'):
            dm.fresh_model(self.source)
        self.auth.assert_not_called()


if __name__ == '__main__':
    unittest.main()

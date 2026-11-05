"""Contracts for literal-only retargeting, without numerical coefficient work."""

from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import hd_field128_gamma_retarget_positive as mode


class RetargetTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.path = Path(self.temporary.name)/'source.json'

    def tearDown(self):
        self.temporary.cleanup()

    def source(self, mutate=None):
        data = dict(active_groups=112, state_dimension=30, target_distance='2/5',
            construction=mode.adaptive.provider.gamma.construction(30),
            source_sha256={}, log2_upper_screen=-99999,
            witness=dict(z='1/4', alpha='1/256', markers=[['1/8', '1/16']],
                         grid=128, state_dimension=30, old_numerical_cap=-88888))
        if mutate:
            mutate(data)
        raw = json.dumps(data).encode()
        self.path.write_bytes(raw)
        return sha256(raw).hexdigest()

    def test_exact_positive_literals_only(self):
        digest = self.source()
        result = mode.retarget(self.path, digest, 112, 80)
        self.assertEqual(result['active_groups'], 80)
        self.assertEqual(result['source_weight'], 112)
        self.assertEqual(result['witness'], dict(z='1/4', alpha='1/256',
            markers=[['1/8', '1/16']], grid=128, state_dimension=30))
        self.assertNotIn('log2_upper_screen', result)
        for key in ('directed', 'selected_weight_certificate',
                    'whole_code_certificate', 'old_numerical_bounds_reused'):
            self.assertIs(result[key], False)
        self.assertEqual(result['source']['sha256'], digest)
        self.assertIn(Path(mode.__file__).name, result['source_sha256'])

    def test_bad_hash_and_original_weight(self):
        digest = self.source()
        for supplied_hash, source_weight in [('0'*64, 112), (digest, 113)]:
            with self.assertRaises(ValueError):
                mode.retarget(self.path, supplied_hash, source_weight, 80)
        for weight in (True, 1, 4096, 80.0):
            with self.assertRaises(ValueError):
                mode.retarget(self.path, digest, 112, weight)

    def test_construction_and_field_mismatch(self):
        mutations = [lambda d: d.update(state_dimension=31),
            lambda d: d.update(target_distance='9/20'),
            lambda d: d.update(construction='another-code'),
            lambda d: d['witness'].update(field_order=256),
            lambda d: d['witness'].update(inner_width=32)]
        for mutate in mutations:
            digest = self.source(mutate)
            with self.assertRaises(ValueError):
                mode.retarget(self.path, digest, 112, 80)

    def test_invalid_markers_and_scalar_parameters(self):
        mutations = [lambda d: d['witness'].update(markers=[]),
            lambda d: d['witness'].update(markers=[['0', '1']]),
            lambda d: d['witness'].update(markers=[['1']]),
            lambda d: d['witness'].update(z='1'),
            lambda d: d['witness'].update(alpha='2'),
            lambda d: d['witness'].update(grid=256)]
        for mutate in mutations:
            digest = self.source(mutate)
            with self.assertRaises(ValueError):
                mode.retarget(self.path, digest, 112, 80)

    def test_conflicting_helper_pin(self):
        digest = self.source(lambda d: d['source_sha256'].update(
            {'hd_field128_gamma_midpoint_adaptive.py': '0'*64}))
        with self.assertRaisesRegex(ValueError, 'Source closure conflict'):
            mode.retarget(self.path, digest, 112, 80)


if __name__ == '__main__':
    unittest.main()

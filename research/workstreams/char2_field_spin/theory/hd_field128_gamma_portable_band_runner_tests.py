"""Pure validation/provenance tests; never construct a numerical AliasBank."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath, PureWindowsPath
import tempfile
import unittest
from unittest.mock import patch

import hd_field128_gamma_portable_band_runner as runner


def seed():
    return dict(schema='gamma-portable-exact-witness-v1',
        field_order=1 << 128, construction=runner.cert.local.construction(30),
        active_groups=68, state_dimension=30, target_distance='2/5',
        witness=dict(z='1/4', alpha='1/8', markers=[['1/2', '3/4'], ['2', '3']],
                     grid=128, state_dimension=30),
        source=dict(path='C:/old/unavailable/original.json', sha256='1'*64),
        origin_role='Exact witness only; no bounds', directed=False, whole_code_certificate=False)


class ValidationTests(unittest.TestCase):
    def test_exact_seed_and_external_ancestry_not_loaded(self):
        value = seed()
        self.assertEqual(runner.validate_seed(value),
                         (F(1, 4), F(1, 8), ((F(1, 2), F(3, 4)), (F(2), F(3)))))
        with tempfile.TemporaryDirectory() as folder, patch.object(runner, 'HERE', Path(folder)):
            raw = (json.dumps(value)+'\n').encode()
            Path(folder, 'seed.json').write_bytes(raw)
            loaded, parameters, pin = runner.load_seed('seed.json', sha256(raw).hexdigest())
            self.assertEqual(loaded, value)
            self.assertEqual(parameters, runner.validate_seed(value))
            self.assertEqual(pin, dict(path='seed.json', sha256=sha256(raw).hexdigest()))
            with self.assertRaises(ValueError):
                runner.load_seed('seed.json', '0'*64)

    def test_wrong_construction_or_scope_rejected(self):
        changes = dict(schema='other', field_order=256, construction='other',
                       state_dimension=31, target_distance='1/3', active_groups=True,
                       directed=True, whole_code_certificate=True)
        for key, value in changes.items():
            with self.subTest(key=key):
                data = seed()
                data[key] = value
                with self.assertRaises(ValueError):
                    runner.validate_seed(data)
        for key, value in (('grid', 64), ('state_dimension', 31), ('z', '1'),
                           ('z', '0'), ('alpha', '2'), ('alpha', '0')):
            data = seed()
            data['witness'][key] = value
            with self.subTest(witness_key=key, value=value), self.assertRaises(ValueError):
                runner.validate_seed(data)

    def test_literal_and_marker_guards(self):
        for value in (0.25, True, '0.25', '2/8', 'nan', '1/0', None):
            data = seed()
            data['witness']['z'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                runner.validate_seed(data)
        for rows in ([], [['0', '1']], [['-1', '1']], [['1']], [[1, '2']],
                     [['1', '2'], ['1', '2']], [['1', str(k+1)] for k in range(129)]):
            data = seed()
            data['witness']['markers'] = rows
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                runner.validate_seed(data)

    def test_portable_names_and_assignment(self):
        for value in ('/tmp/output', 'C:/output', 'C:\\output', '../output',
                      'a/b', 'a\\b', '.', '..', 'output.', ''):
            with self.subTest(value=value), self.assertRaises(ValueError):
                runner.relative_name(value)
        self.assertEqual(runner.relative_name('seed-v1.json', suffix='.json'), 'seed-v1.json')
        for args in ((65, 100, 512), (68, 68, 2048)):
            runner.validate_band(*args)
        for args in ((64, 68, 512), (68, 101, 512), (70, 68, 512),
                     (True, 68, 512), (65, 68, 1024), (65, 68, 512.0)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                runner.validate_band(*args)

    def test_relative_metadata_matches_typed_audit_contract(self):
        pins = {runner.NAME: '2'*64, 'helper.py': '3'*64}
        provenance = dict(path='seed.json', sha256='4'*64)
        saved = runner.driver_record(pins, provenance, 65, 72, 512, 68)
        self.assertEqual(saved['witness_source'], provenance)
        self.assertEqual(saved['requested_weight_range'], [65, 72])
        self.assertEqual((saved['cell_failure_bits'], saved['aggregate_failure_bits']), (72, 54))
        self.assertEqual(saved['numerical_provider'], 'hd_field128_gamma_midpoint_weight_bands.py')
        self.assertEqual(saved['source_sha256'][runner.NAME], saved['driver']['sha256'])
        self.assertEqual((saved['precision'], saved['coefficient_bits'], saved['grid']), (256, 80, 128))
        self.assertEqual((saved['max_build_seconds'], saved['max_prepare_seconds'],
                          saved['max_cover_seconds']), (600, 180, 60))
        for pin in (saved['driver'], saved['witness_source']):
            self.assertFalse(PurePosixPath(pin['path']).is_absolute())
            self.assertFalse(PureWindowsPath(pin['path']).is_absolute())
            self.assertEqual(len(PurePosixPath(pin['path']).parts), 1)
        self.assertFalse(saved['directed'])
        self.assertFalse(saved['whole_code_certificate'])
        self.assertTrue(saved['fresh_bank_required'])
        self.assertFalse(saved['fresh_bank_rebuilt'])

    def test_runtime_metadata_is_explicit_and_not_a_source_link(self):
        with patch.dict(runner.os.environ, {'PYTHONPATH': '/execution/example',
                                          'OPENBLAS_NUM_THREADS': '1'}):
            metadata = runner.runtime_metadata()
        self.assertEqual(metadata['PYTHONPATH'], '/execution/example')
        self.assertEqual(metadata['OPENBLAS_NUM_THREADS'], '1')
        for key in ('python_version', 'platform', 'python_flint_version',
                    'flint_version', 'numpy_version', 'scipy_version'):
            self.assertIsInstance(metadata[key], str)
            self.assertTrue(metadata[key])
        self.assertNotIn('path', metadata)
        self.assertNotIn('sha256', metadata)
        self.assertIn('not source authentication', metadata['scope'])

    def test_source_seed_mutation_and_per_write_guard(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(runner, 'HERE', Path(folder)):
            a = Path(folder, 'runner.py')
            b = Path(folder, 'seed.json')
            a.write_bytes(b'source')
            b.write_bytes(b'{}')
            pins = {'runner.py': sha256(a.read_bytes()).hexdigest()}
            pin = dict(path=b.name, sha256=sha256(b.read_bytes()).hexdigest())
            runner.check_pins(pins, pin)
            b.write_bytes(b'{"changed":true}')
            with self.assertRaises(RuntimeError):
                runner.check_pins(pins, pin)
            b.write_bytes(b'{}')
            a.write_bytes(b'changed source')
            with self.assertRaises(RuntimeError):
                runner.check_pins(pins, pin)
            destination = Path(folder, 'artifact.json')
            usage = type('Usage', (), {'free': runner.RESERVE-1})()
            with patch.object(runner.shutil, 'disk_usage', return_value=usage):
                with self.assertRaises(RuntimeError):
                    runner.write_new(destination, {'validation_only': True})
            self.assertFalse(destination.exists())
            usage.free = runner.RESERVE
            with patch.object(runner.shutil, 'disk_usage', return_value=usage):
                self.assertIsNone(runner.write_new(destination, {'validation_only': True}))
                with self.assertRaises(FileExistsError):
                    runner.write_new(destination, {'validation_only': False})
            self.assertEqual(json.loads(destination.read_text()), {'validation_only': True})


if __name__ == '__main__':
    unittest.main()

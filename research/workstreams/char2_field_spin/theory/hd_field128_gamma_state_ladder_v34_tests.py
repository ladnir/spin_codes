"""Memory-only serialization regression; no coefficients or output artifacts."""
import inspect
import io
import json
from pathlib import Path
import unittest
from unittest.mock import Mock, patch

import numpy as np
import hd_field128_gamma_state_ladder_v33 as old
import hd_field128_gamma_state_ladder_v34 as pilot


class Checks(unittest.TestCase):
    def test_numpy_diagnostics(self):
        value = {'largest_profile': {'full_even': np.int64(1581), 'full_odd': np.int64(1614)},
                 'diagnostics': [np.uint64(4198401), np.float64(-3.5), np.bool_(True)],
                 'tuple': (np.int32(30), np.float32(1.25))}
        ready = pilot.json_ready(value)
        self.assertEqual(ready['largest_profile'], {'full_even': 1581, 'full_odd': 1614})
        self.assertIs(type(ready['largest_profile']['full_even']), int)
        self.assertIs(type(ready['diagnostics'][1]), float)
        self.assertIs(type(ready['diagnostics'][2]), bool)
        self.assertEqual(json.loads(json.dumps(ready, allow_nan=False)), ready)
        self.assertIsInstance(value['largest_profile']['full_even'], np.int64)

    def test_writer_reproduces_real_failure_without_writes(self):
        value = {'comparisons': [{'largest_profile': {'full_even': np.int64(1735)}}]}
        with self.assertRaises(TypeError):
            json.dumps(value)
        writer = object.__new__(pilot.Writer)
        writer.directory, writer.written, writer.pins = pilot.HERE/'not-created', 0, {}
        stream = io.BytesIO()
        context = Mock()
        context.__enter__ = Mock(return_value=stream)
        context.__exit__ = Mock(return_value=False)
        with patch.object(Path, 'open', return_value=context) as opening, patch.object(
                pilot.shutil, 'disk_usage', return_value=Mock(free=pilot.RESERVE+pilot.OUTPUT_LIMIT)):
            receipt = writer.write('result.json', value)
        opening.assert_called_once_with('xb')
        raw = stream.getvalue()
        self.assertEqual(json.loads(raw), {'comparisons': [{'largest_profile': {'full_even': 1735}}]})
        self.assertEqual(writer.written, len(raw))
        self.assertEqual(receipt['bytes'], len(raw))

    def test_reject_arrays_and_nonfinite(self):
        for value in (np.array([1, 2]), np.float64(float('nan')), np.float64(float('inf'))):
            with self.assertRaises((TypeError, ValueError)):
                json.dumps(pilot.json_ready(value), allow_nan=False)

    def test_frozen_arithmetic_and_scope_unchanged(self):
        # The entire computation body differs only in the new output schema tag.
        self.assertEqual(inspect.getsource(pilot._run).replace(
            'gamma-state-ladder-conditional-FLOAT-v34', 'gamma-state-ladder-conditional-FLOAT-v33'),
            inspect.getsource(old._run))
        for state in (30, 31, 32):
            self.assertEqual(pilot.state_metadata(state), old.state_metadata(state))
        self.assertEqual(pilot.RESERVE, old.RESERVE)
        self.assertEqual(pilot.OUTPUT_LIMIT, old.OUTPUT_LIMIT)

    def test_both_source_identities_in_transitive_closure(self):
        pins = pilot.source_pins()
        self.assertEqual(pins[Path(old.__file__).name], pilot.FROZEN_V33_SHA256)
        self.assertIn(Path(pilot.__file__).name, pins)
        for key, digest in old.source_pins().items():
            self.assertEqual(pins[key], digest)


if __name__ == '__main__':
    unittest.main()

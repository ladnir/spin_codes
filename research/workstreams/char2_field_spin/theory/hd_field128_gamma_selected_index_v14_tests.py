"""Tiny fixture checks for selected-event admission; no coefficients."""

from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import hd_field128_gamma_selected_index_v14 as selected
from hd_field128_gamma_near_face_index_update_tests import audited_fixture


class SelectedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = TemporaryDirectory()
        summary, unused, unused_pin = audited_fixture(Path(cls.temporary.name))
        context = json.loads(Path(summary['context']['path']).read_bytes())
        row = json.loads(Path(summary['weight_receipts'][0]['path']).read_bytes())
        cls.record = dict(context, **row)
        cls.record['source_pins_unchanged_at_exit'] = True
        cls.pins = cls.record['source_sha256']

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_positive_complete_selected_admission(self):
        bound = selected.complete_selected(self.record, 1, self.pins)
        self.assertLessEqual(bound, -54)
        self.assertIn(Path(selected.__file__).name, selected.source_pins())
        self.assertIn('../HD_FIELD128_MIDPOINT_COUNTS.md', selected.source_pins())

    def test_scope_row_prefix_and_bound_fail_closed(self):
        for change in (lambda x: x.update(complete_selected_weight=False),
                       lambda x: x.update(state_dimension=31),
                       lambda x: x['rows'].pop(),
                       lambda x: x['rows'][0].update(even_full=9),
                       lambda x: x.update(log2_upper_rational='-10000000')):
            row = deepcopy(self.record)
            change(row)
            with self.assertRaises(ValueError):
                selected.complete_selected(row, 1, self.pins)
        with patch.object(selected.index, 'scalar_ceiling', return_value=selected.F(-53)):
            with self.assertRaisesRegex(ValueError, '54-bit'):
                selected.complete_selected(self.record, 1, self.pins)

    def test_only_requested_near_audit_set_allowed(self):
        for audits in ({}, {128: 'x'}, {128: 'x', 384: 'y', 768: 'z'}):
            with self.assertRaisesRegex(ValueError, 'three requested'):
                selected.extend_v13('not-read.json', audits, 'not-read.json')


if __name__ == '__main__':
    unittest.main()

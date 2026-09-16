"""Retained-grid shape, evidence scope, and generated-figure regression tests."""
import copy
import json
import unittest
from unittest.mock import patch

import build_imt_parameter_figures as figures


class IMTParameterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base = figures.evidence.BASE
        cls.grid = json.loads((base / figures.GRID_PIN[0]).read_text())
        cls.replay = json.loads((base / figures.REPLAY_PIN[0]).read_text())

    def validate_changed(self, change):
        grid, replay = copy.deepcopy(self.grid), copy.deepcopy(self.replay)
        change(grid, replay)
        with self.assertRaises(ValueError):
            figures.validate_grid(grid, replay)

    def test_coverage(self):
        self.assertEqual(len(figures.validate_grid(self.grid, self.replay)), 130)

    def test_missing_row(self):
        self.validate_changed(lambda g, r: g['rows'].pop())

    def test_duplicate_row(self):
        self.validate_changed(lambda g, r: g['rows'].append(g['rows'][0]))

    def test_no_full_claim(self):
        self.validate_changed(lambda g, r: g.update(full_distance_proved=True))

    def test_wrong_replay(self):
        self.validate_changed(lambda g, r: r.update(producer_sha256='0'*64))

    def test_wrong_inner(self):
        self.validate_changed(lambda g, r: g['maps']['t64_s20'].update(transvection_rounds=2))

    def test_wrong_cutoff(self):
        self.validate_changed(lambda g, r: g['rows'][0].update(bad_weight=0))

    def test_changed_nested_map(self):
        self.validate_changed(lambda g, r: g['maps']['t64_s20']['feedback_columns'].__setitem__(0, 1))

    def test_authenticated_generated_figures(self):
        report = figures.check()
        self.assertEqual(report['matched_full_anchors'], 5)
        self.assertFalse(report['full_grid_certified'])

    def test_tampered_receipt(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'wrong.json').write_text('{}')
            with patch.object(figures.evidence, 'BASE', root):
                with self.assertRaisesRegex(ValueError, 'Receipt changed'):
                    figures.evidence.authenticate(('wrong.json', '0'*64), {})


if __name__ == '__main__':
    unittest.main()

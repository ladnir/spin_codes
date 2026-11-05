"""Count-rescaling arithmetic and scope guards for midpoint cell covers."""

from fractions import Fraction as F
from pathlib import Path
import unittest
from unittest.mock import patch

import hd_field128_gamma_midpoint_cell_cover as mode
from hd_field128_gamma_cell_pinned_batch_tests import evaluator as frozen_fixture


def receipt(count=16, ceiling='-5'):
    return dict(complete_cell=True, directed=True, status='complete',
                count=dict(support_count_upper=count, method='fixture'),
                log2_upper_rational=ceiling, log2_upper_enclosure=ceiling,
                zero=False, failure_bits=7, meets_cell_budget=False)


def evaluator():
    original, bank = frozen_fixture()
    regional = mode.frozen.regional
    rows = dict(enumerate(original.table.levels[0]))
    with patch.object(regional.restored, 'RegionalCache') as cache:
        cache.return_value.regional.return_value = rows
        result = mode.Evaluator(bank, 4)
    bank.check_pins.reset_mock()
    return result, bank


class RefineTests(unittest.TestCase):
    def test_exact_power_of_two_rescaling(self):
        original = receipt()
        row = mode.refine_receipt(original, 1, 256)
        self.assertLessEqual(F(row['log2_upper_rational']), -9+F(1, 10**12))
        self.assertGreaterEqual(F(row['log2_upper_rational']), -9)
        self.assertTrue(row['meets_cell_budget'])
        self.assertEqual(original, receipt())
        self.assertEqual(row['count']['support_count_upper'], 1)

    def test_looser_zero_and_guards(self):
        row = mode.refine_receipt(receipt(), 128, 256)
        self.assertEqual(row['log2_upper_rational'], '-5')
        self.assertFalse(row['meets_cell_budget'])
        row = mode.refine_receipt(receipt(), 0, 256)
        self.assertTrue(row['zero'])
        self.assertTrue(row['meets_cell_budget'])
        self.assertIsNone(row['log2_upper_rational'])
        for value in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                mode.refine_receipt(receipt(), value, 256)
        with self.assertRaises(ValueError):
            mode.refine_receipt(dict(receipt(), complete_cell=False), 1, 256)
        with self.assertRaises(ValueError):
            mode.refine_receipt(receipt(count=0), 0, 256)
        self.assertIn('../HD_FIELD128_MIDPOINT_COUNTS.md', mode.source_pins())

    def test_public_cover_uses_new_caps_with_amortized_checks(self):
        value, bank = evaluator()
        result = value.run_cover(max_nodes=10, max_seconds=30., cell_seconds=5.)
        self.assertTrue(result['complete_selected_weight'])
        self.assertTrue(result['batch_source_verified_at_exit'])
        self.assertEqual(bank.check_pins.call_count, 2)
        for node in result['nodes']:
            if node['status'] == 'accepted':
                row = node['receipt']
                self.assertIn('midpoint_rectangle_count_upper', row)
                self.assertLessEqual(row['count']['support_count_upper'], row['unrefined_count']['support_count_upper'])
        self.assertIn('hd_field128_gamma_midpoint_boxes.py', result['batch_source_sha256'])

    def test_new_source_mutation_rejects_preparation_and_cover(self):
        read = Path.read_bytes
        changed = [False]
        def altered(path):
            raw = read(path)
            return raw+b'\n# mock mutation\n' if changed[0] and path.name == Path(mode.boxes.__file__).name else raw
        value, bank = evaluator()
        def checkpoint(row):
            changed[0] = True
        with patch.object(Path, 'read_bytes', altered):
            with self.assertRaisesRegex(RuntimeError, 'pinned'):
                value.run_cover(max_nodes=10, max_seconds=30., checkpoint=checkpoint)
        self.assertFalse(value._batch_active)
        changed[0] = False
        original, bank = frozen_fixture()
        original_boxes = mode.boxes.MidpointBoxes
        def mutate_after_preparation(*args, **kwargs):
            result = original_boxes(*args, **kwargs)
            changed[0] = True
            return result
        with patch.object(Path, 'read_bytes', altered), patch.object(mode.boxes, 'MidpointBoxes', mutate_after_preparation):
            with patch.object(mode.frozen.regional.restored, 'RegionalCache') as cache:
                cache.return_value.regional.return_value = dict(enumerate(original.table.levels[0]))
                with self.assertRaisesRegex(RuntimeError, 'pinned'):
                    mode.Evaluator(bank, 4)


if __name__ == '__main__':
    unittest.main()

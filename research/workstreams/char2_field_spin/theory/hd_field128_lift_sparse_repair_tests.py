"""Tiny q256 repair/provenance checks; no extension-field source-weight gate."""

from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

import hd_field128_lift_sparse_bands as bands
import hd_field128_lift_sparse_repair as repair


class SparseRepairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt = bands.run_schedule(first=1, last=3, q=256, max_seconds=20,
                                         witnesses={c: (F(7, 8), F(1)) for c in (1, 2, 4)})
        cls.directory = tempfile.TemporaryDirectory(prefix='field128_sparse_repair_')
        cls.path = Path(cls.directory.name)/'input.json'
        cls.path.write_text(json.dumps(cls.receipt), encoding='utf-8')

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()

    def test_validated_input_keeps_every_weight_once(self):
        rows = repair.validate_schedule(self.receipt, first=1, last=3, q=256)
        self.assertEqual([row['active_groups'] for row in rows], [1, 2, 3])
        self.assertEqual(rows[2]['origin'], dict(kind='input-band-row', band_index=2, row_index=0))
        self.assertEqual(self.receipt['completed_bands'][2]['weight_receipts'][0]['active_groups'], 3)
        self.assertNotIn('origin', self.receipt['completed_bands'][2]['weight_receipts'][0])

    def test_rejects_incompatible_or_incomplete_receipts(self):
        changes = [
            lambda r: r.update(field_order=65536),
            lambda r: r.update(directed=False),
            lambda r: r.update(whole_code_certificate=True),
            lambda r: r.update(target_distance='2/5'),
            lambda r: r.update(pending_weight_ranges=[[3, 3]]),
            lambda r: r.update(source_pins_unchanged_at_exit=False),
            lambda r: r['source_sha256'].update({'hd_field128_lift_sparse_bands.py': '0'*64}),
            lambda r: r['source_sha256'].pop('hd_field128_lift_certificate.py'),
            lambda r: r['completed_bands'][1].update(weight_range=[1, 2]),
            lambda r: r['completed_bands'][1].update(coefficient_maximum=1),
            lambda r: r['completed_bands'][1].update(pair_swaps=False),
            lambda r: r['completed_bands'][1].update(coefficient_method='FLOAT'),
            lambda r: r['completed_bands'][1]['weight_receipts'].clear(),
            lambda r: r['completed_bands'][1]['weight_receipts'][0].update(profiles=0),
            lambda r: r['completed_bands'][1]['weight_receipts'][0].update(covers_all_supports_at_weight=False),
            lambda r: r['completed_bands'][1]['weight_receipts'][0]['witness'].update(alpha='1/2'),
            lambda r: r['completed_bands'][1]['weight_receipts'][0].update(log2_upper_rational='-999999999'),
            lambda r: r['completed_bands'][1].update(log2_upper_rational='-999999999'),
            lambda r: r.update(completed_union_log2_upper_rational='-999999999'),
        ]
        for index, change in enumerate(changes):
            with self.subTest(index=index):
                receipt = deepcopy(self.receipt)
                change(receipt)
                with self.assertRaises(ValueError):
                    repair.validate_schedule(receipt, first=1, last=3, q=256)

    def test_own_weight_replays_and_best_bound_selection(self):
        snapshots = []
        result = repair.repair_schedule(self.path, first=1, last=3, q=256, repairs=(1, 3),
                                        max_seconds=20, checkpoint=snapshots.append)
        self.assertEqual(result['pending_repairs'], [])
        self.assertTrue(result['all_requested_repairs_completed'])
        self.assertEqual(result['selected_weight_count'], 3)
        self.assertEqual(result['duplicate_weight_count'], 0)
        self.assertEqual([row['active_groups'] for row in result['selected_weight_receipts']], [1, 2, 3])
        original = repair.validate_schedule(self.receipt, first=1, last=3, q=256)
        candidates = {row['active_groups']: F(row['log2_upper_rational']) for row in original}
        for row in result['completed_repairs']:
            r = row['active_groups']
            z, alpha = bands.default_witness(r, 256)
            self.assertEqual((row['witness']['z'], row['witness']['alpha']), (str(z), str(alpha)))
            self.assertEqual(row['witness']['coefficient_maximum'], r)
            candidates[r] = min(candidates[r], F(row['log2_upper_rational']))
        self.assertEqual([F(row['log2_upper_rational']) for row in result['selected_weight_receipts']],
                         [candidates[r] for r in (1, 2, 3)])
        self.assertEqual(F(result['log2_upper_rational']), repair.union_ceiling(result['selected_weight_receipts'], 256))
        self.assertEqual(result['input_provenance']['sha256'], sha256(self.path.read_bytes()).hexdigest())
        self.assertEqual(result['completed_repair_provenance'], repair.repair_provenance(result['completed_repairs']))
        self.assertEqual(len(result['completed_repair_provenance']), 2)
        self.assertTrue(result['source_pins_unchanged_at_exit'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(snapshots[0]['pending_repairs'], [1, 3])
        self.assertEqual(snapshots[1]['pending_repairs'], [3])

    def test_timeout_preserves_original_coverage(self):
        result = repair.repair_schedule(self.path, first=1, last=3, q=256, repairs=(1, 3), max_seconds=1e-12)
        self.assertTrue(result['complete_requested_range'])
        self.assertFalse(result['all_requested_repairs_completed'])
        self.assertEqual(result['pending_repairs'], [1, 3])
        self.assertEqual(result['completed_repairs'], [])
        self.assertEqual(result['failure']['error_type'], 'TimeoutError')
        self.assertEqual(result['selected_weight_count'], 3)
        self.assertTrue(all(row['origin']['kind'] == 'input-band-row' for row in result['selected_weight_receipts']))

    def test_strict_json_and_source_closure(self):
        path = Path(self.directory.name)/'invalid.json'
        for contents in ('{"directed":true,"directed":false}', '{"x":NaN}'):
            path.write_text(contents, encoding='utf-8')
            with self.assertRaises(ValueError):
                repair.load_receipt(path)
        pins = repair.source_pins()
        for name in ('hd_field128_lift_sparse_repair.py', 'hd_field128_lift_sparse_bands.py',
                     'hd_field128_lift_certificate.py', 'hd_common_refresh_principal.py'):
            self.assertIn(name, pins)
        for weights in ((1, 1), (0,), (4,)):
            with self.assertRaises(ValueError):
                repair.repair_schedule(self.path, first=1, last=3, q=256, repairs=weights)


if __name__ == '__main__':
    unittest.main()

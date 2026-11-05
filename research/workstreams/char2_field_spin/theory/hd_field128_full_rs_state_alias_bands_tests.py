"""Tiny q256 band compatibility and mock timeouts; no full-field gate."""

from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import hd_field128_full_rs_state_alias_bands as bands


class StateAliasBandTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = bands.provider.AliasBank(F(7, 8), F(1, 4), ((F(1, 512), F(1, 512)),),
                                            q=256, state=24, grid=64, max_seconds=60.)

    def test_one_cache_replays_frozen_cover_for_every_tiny_weight(self):
        saved = []
        with patch.object(bands.cover, 'RegionalCache', wraps=bands.cover.RegionalCache) as cache:
            result = bands.band_receipt(self.bank, 1, 3, (F(1, 256), F(1, 128)),
                                        max_seconds=30., checkpoint=saved.append)
        cache.assert_called_once_with(self.bank)
        frozen = bands.cover.band_receipt(self.bank, 1, 3, (F(1, 256), F(1, 128)), max_seconds=30.)
        self.assertEqual([row['active_groups'] for row in saved], [1, 2, 3])
        for row, control in zip(saved, frozen['weight_receipts']):
            for key in ('active_groups', 'regional_histograms', 'log2_upper_rational', 'state_potential_ratio'):
                self.assertEqual(row[key], control[key])
            self.assertEqual(row['construction'], bands.provider.local.construction(24))
            self.assertEqual(row['state_dimension'], 24)
            self.assertEqual(row['target_distance'], '2/5')
            self.assertIn('hd_field128_full_rs_state_alias_bands.py', row['source_sha256'])
            self.assertTrue(row['directed'])
            self.assertFalse(row['whole_code_certificate'])
        self.assertEqual(result['log2_upper_rational'], frozen['log2_upper_rational'])
        self.assertEqual(result['pending_weight_ranges'], [])
        self.assertTrue(result['every_integer_weight_evaluated'])
        self.assertTrue(result['covers_all_supports_in_range'])

    def test_later_timeout_retains_full_immutable_first_row(self):
        original = bands.cover.RegionalCache.regional
        def timeout(cache, weight, budget):
            if weight == 2:
                raise TimeoutError('mock later row deadline')
            return original(cache, weight, budget)
        saved = []
        with patch.object(bands.cover.RegionalCache, 'regional', new=timeout):
            result = bands.band_receipt(self.bank, 1, 3, (F(1, 256),), checkpoint=saved.append)
        self.assertEqual([row['active_groups'] for row in saved], [1])
        self.assertEqual(saved[0]['state_potential_ratio'], '1/256')
        self.assertEqual(result['pending_weight_ranges'], [[2, 3]])
        self.assertEqual(result['completed_weight_range'], [1, 1])
        self.assertFalse(result['covers_all_supports_in_range'])
        self.assertIsNone(result['log2_upper_rational'])
        self.assertIsNotNone(result['completed_union_log2_upper_rational'])
        self.assertEqual(result['failure']['active_groups'], 2)

    def test_expired_budget_leaves_every_requested_weight_pending(self):
        result = bands.band_receipt(self.bank, 2, 4, (F(1, 256),), max_seconds=1e-12)
        self.assertEqual(result['weights_evaluated'], 0)
        self.assertEqual(result['pending_weight_ranges'], [[2, 4]])
        self.assertFalse(result['complete_requested_range'])

    def test_guards_and_exclusive_receipts(self):
        for first, last, maximum in ((0, 1, 1), (2, 1, 2), (1, 4, 3), (1, 8193, 8192)):
            with self.assertRaises(ValueError):
                bands.band_receipt(self.bank, first, last, (F(1),), max_weights=maximum)
        for ratios in ((), (F(0),), (.5,)):
            with self.assertRaises(ValueError):
                bands.band_receipt(self.bank, 1, 1, ratios)
        with self.assertRaises(ValueError):
            bands.band_receipt(object(), 1, 1, (F(1),))
        with TemporaryDirectory() as directory:
            path = Path(directory)/'weight_0001.json'
            receipt = bands.write_new(path, {'scope': 'mock only', 'directed': True})
            self.assertEqual(receipt['sha256'], sha256(path.read_bytes()).hexdigest())
            with self.assertRaises(FileExistsError):
                bands.write_new(path, {})


if __name__ == '__main__':
    unittest.main()

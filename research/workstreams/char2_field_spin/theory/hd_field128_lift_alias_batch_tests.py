"""Tiny batch/provenance tests; all bank calculations below are mocked."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import hd_field128_lift_alias_batch as batch
import hd_field128_lift_alias_cover as cover


MANIFEST = Path(__file__).with_name('hd_shifted_closed_cover_manifest.json')


def fake_receipt(bank,first,last,ratios,**kwargs):
    rows = [dict(active_groups=r,log2_upper_rational='-100') for r in range(first,last+1)]
    return dict(weight_receipts=rows,complete_requested_range=True,
                completed_union_log2_upper_rational=cover.union_ceiling(rows,256),failure=None)


class AliasBatchTests(unittest.TestCase):
    def test_all33_old_geometry_jobs_and_rounded_exact_witnesses(self):
        jobs,pin = batch.load_jobs(MANIFEST,257,4067)
        self.assertEqual(len(jobs),33)
        self.assertEqual([r for job in jobs for r in range(job['weight_range'][0],job['weight_range'][1]+1)],
                         list(range(257,4068)))
        for job in jobs:
            artifact,_ = cover.read_json(job['witness']['source']['path'])
            exact = cover.witness_from_artifact(artifact,old=True)
            rounded = cover.witness_from_artifact(job['witness'])
            self.assertEqual(rounded['markers'],exact['markers'])
            self.assertEqual(rounded['ratio'],exact['ratio'])
            self.assertEqual(rounded['z'],F(str(float(exact['z']))))
            self.assertEqual(rounded['alpha'],F(str(float(exact['alpha']))))
            self.assertFalse(job['witness']['directed'])
            self.assertIn('new arbitrary exact witnesses',job['witness']['transform']['equality_claim'])

    def test_max_bands_saves_completed_band_immediately(self):
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)/'batch'
            with (patch.object(cover.interval,'AliasBank',return_value=object()) as build,
                  patch.object(cover,'band_receipt',side_effect=fake_receipt)):
                result = batch.run_batch(MANIFEST,directory,first=257,last=280,max_bands=1,max_seconds=20)
            self.assertEqual(build.call_count,1)
            self.assertEqual(result['covered_weight_count'],23)
            self.assertEqual(result['pending_weight_ranges'],[[280,280]])
            self.assertFalse(result['complete_requested_range'])
            self.assertTrue((directory/'band_000_257_279_receipt.json').is_file())
            self.assertTrue((directory/'band_000_257_279_witness.json').is_file())
            self.assertTrue((directory/'snapshot_001.json').is_file())
            self.assertTrue((directory/'final_summary.json').is_file())
            pin = result['completed_bands'][0]['receipt']
            self.assertEqual(pin['sha256'],sha256(Path(pin['path']).read_bytes()).hexdigest())
            with self.assertRaises(FileExistsError):
                batch.write_new(directory/'final_summary.json',{})

    def test_full_tiny_mock_batch_sums_each_weight_once(self):
        with TemporaryDirectory() as temporary:
            with (patch.object(cover.interval,'AliasBank',return_value=object()),
                  patch.object(cover,'band_receipt',side_effect=fake_receipt)):
                result = batch.run_batch(MANIFEST,Path(temporary)/'batch',first=257,last=280,max_bands=2,max_seconds=20)
            self.assertTrue(result['complete_requested_range'])
            self.assertEqual(result['covered_weight_count'],24)
            self.assertEqual(result['duplicate_weight_count'],0)
            self.assertTrue(result['requested_range_meets_target'])
            self.assertEqual(result['log2_upper_rational'],cover.union_ceiling(result['weight_receipts'],256))

    def test_partial_band_receipt_survives_and_remainder_is_pending(self):
        def partial(bank,first,last,ratios,**kwargs):
            result = fake_receipt(bank,first,first+1,ratios)
            result.update(complete_requested_range=False,failure=dict(error_type='TimeoutError'))
            return result
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)/'batch'
            with (patch.object(cover.interval,'AliasBank',return_value=object()),
                  patch.object(cover,'band_receipt',side_effect=partial)):
                result = batch.run_batch(MANIFEST,directory,first=257,last=280,max_bands=2,max_seconds=20)
            self.assertEqual(result['covered_weight_count'],2)
            self.assertEqual(result['pending_weight_ranges'],[[259,280]])
            self.assertEqual(result['failure']['reason'],'partial-band')
            self.assertIsNone(result['log2_upper_rational'])
            saved = json.loads((directory/'band_000_257_279_receipt.json').read_text())
            self.assertEqual(len(saved['weight_receipts']),2)

    def test_deadline_before_first_bank_and_invalid_parameters(self):
        with TemporaryDirectory() as temporary:
            with patch.object(cover.interval,'AliasBank',return_value=object()) as build:
                result = batch.run_batch(MANIFEST,Path(temporary)/'batch',first=257,last=279,max_seconds=1e-12)
            self.assertEqual(build.call_count,0)
            self.assertEqual(result['pending_weight_ranges'],[[257,279]])
            self.assertEqual(result['failure']['error_type'],'TimeoutError')
        with self.assertRaises(ValueError):
            batch.load_jobs(MANIFEST,256,279)
        with self.assertRaises(ValueError):
            batch.run_batch(MANIFEST,'unused',factors=(.5,))
        self.assertIn('hd_field128_lift_alias_cover.py',batch.source_pins())


if __name__=='__main__':
    unittest.main()

"""Orchestration guards; these tests do not certify numerical bounds."""
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

import hd_field128_gamma_root_v22_runner_20261010 as runner


class RootRunnerTests(unittest.TestCase):
    def test_scope(self):
        for lo, hi in ((65, 100), (140, 240), (273, 383)):
            self.assertTrue(runner.assigned(lo, hi))
        for lo, hi in ((64, 80), (100, 140), (240, 273), (383, 384),
                       (80, 79), (True, 80)):
            self.assertFalse(runner.assigned(lo, hi))

    def fixture(self, folder):
        path = Path(folder) / 'source.json'
        row = dict(active_groups=80, state_dimension=30, target_distance='2/5',
                   construction=runner.cert.local.construction(30), source_sha256={},
                   witness=dict(z='1/4', alpha='1/128', marker_log_odds_exact=[['-4', '-5']]))
        raw = json.dumps(row).encode()
        path.write_bytes(raw)
        return path, sha256(raw).hexdigest()

    def test_exact_literal_conversion_and_noncertificate(self):
        with tempfile.TemporaryDirectory() as folder:
            path, digest = self.fixture(folder)
            result = runner.positive_proposal(path, digest, 80)
            self.assertEqual(result['witness']['markers'], [['1/16', '1/32']])
            self.assertFalse(result['directed'])
            self.assertFalse(result['whole_code_certificate'])
            self.assertEqual(result['source']['sha256'], digest)

    def test_source_hash_and_original_weight(self):
        with tempfile.TemporaryDirectory() as folder:
            path, digest = self.fixture(folder)
            for bad_hash, weight in (('0'*64, 80), (digest, 81), (digest, 384)):
                with self.assertRaises(ValueError):
                    runner.positive_proposal(path, bad_hash, weight)

    def test_bank_type_guard(self):
        with self.assertRaises(TypeError):
            runner.run_band(object(), 65, 80, 'unused', {})


if __name__ == '__main__':
    unittest.main()

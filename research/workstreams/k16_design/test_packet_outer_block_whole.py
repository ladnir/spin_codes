"""Light receipt validation checks; numerical data are optional local fixtures."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from flint import ctx
import packet_outer_block_transfer as transfer
import packet_outer_block_whole as whole

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / 'tmp/rs16-s22-k262144-whole-fresh-p256.json'
SPARSE = {m: ROOT / f'tmp/rs-s22-outer-x{m}-sparse-p256.json' for m in (2, 4)}


class WholeArithmeticTests(unittest.TestCase):
    def test_exact_dyadic_sum(self):
        value = transfer._sum_dyadics([(7, -21), (17, -4096), (13, 3)])
        expected = sum(transfer._fraction(v) for v in [(7, -21), (17, -4096), (13, 3)])
        self.assertEqual(transfer._fraction(value), expected)

    def test_reject_insufficient_precision_before_reading(self):
        with self.assertRaises(ValueError):
            whole.assemble('missing-source', 'missing-sparse', multiplier=2, precision=192)


@unittest.skipUnless(SOURCE.exists() and all(p.exists() for p in SPARSE.values()),
                     'retained local numerical receipts not present')
class LocalReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source, cls.endpoints, cls.choices = transfer.authenticate(SOURCE)

    def test_authenticate_all_tail_endpoints_and_source_pins(self):
        self.assertEqual(set(self.endpoints), set(range(3, 2049)))
        self.assertEqual(len(self.source['source_sha256']), 103)
        self.assertEqual(self.choices[89]['tilt'], '17/200')

    def test_both_full_assemblies_close_and_restore_precision(self):
        previous = ctx.prec
        for m in (2, 4):
            result = whole.assemble(SOURCE, SPARSE[m], multiplier=m)
            self.assertEqual(set(result['occupancy_uppers']), {str(q) for q in range(1, 2049)})
            self.assertTrue(result['whole_code_certificate'])
            self.assertTrue(result['derived_certificate'])
            self.assertFalse(result['fresh_replay'])
            self.assertEqual(result['tail_source']['source']['sha256'], transfer.BASE_SHA256)
            self.assertEqual(ctx.prec, previous)

    def test_wrong_outer_rejected(self):
        with self.assertRaises(ValueError):
            whole.authenticate_sparse(SPARSE[2], multiplier=4, source_whole=self.source)

    def test_sparse_scope_mutations_rejected(self):
        original = json.loads(SPARSE[2].read_text())
        mutations = [
            ('cutoff', original['cutoff']+1), ('fresh_computation', False),
            ('occupancy_covered', [1]), ('beta', '1'), ('state_bits', 20),
            ('source_sha256', {}), ('map_record', {}),
            ('occupancy_choices', {'1': '0', '2': '.0016'}),
            ('occupancy_uppers', {'1': [0, 0], '2': [1, -90]}),
        ]
        with tempfile.TemporaryDirectory(prefix='spin-outer-transfer-test-') as folder:
            target = Path(folder) / 'mutated.json'
            for key, value in mutations:
                mutated = copy.deepcopy(original)
                mutated[key] = value
                target.write_text(json.dumps(mutated))
                with self.subTest(field=key), self.assertRaises(ValueError):
                    whole.authenticate_sparse(target, multiplier=2, source_whole=self.source)

    def test_source_hash_anchor_rejects_reformatting(self):
        with tempfile.TemporaryDirectory(prefix='spin-outer-transfer-test-') as folder:
            target = Path(folder) / 'source.json'
            target.write_text(json.dumps(self.source))
            with self.assertRaises(ValueError):
                transfer.authenticate(target)


if __name__ == '__main__':
    unittest.main()

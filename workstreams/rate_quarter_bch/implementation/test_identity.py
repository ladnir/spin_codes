"""Check implementation provenance and equality with the certified constituents."""
import hashlib
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
import smaller_outer


class Identity(unittest.TestCase):
    def setUp(self):
        self.manifest = json.loads((HERE/'generated/MANIFEST.json').read_text())

    def test_sources(self):
        for name, expected in self.manifest['source_sha256'].items():
            self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(), expected, name)
        for name, expected in self.manifest['generated_sha256'].items():
            self.assertEqual(hashlib.sha256((HERE/'generated'/name).read_bytes()).hexdigest(), expected, name)

    def test_outer_span(self):
        self.assertEqual(self.manifest['outer'], smaller_outer.construction())
        rows = [int(r,16) for r in self.manifest['generator_rows_hex']]
        self.assertEqual(smaller_outer.r.rank(rows), 32)
        raw = [int(r,16) for r in self.manifest['outer']['generator_rows_hex']]
        self.assertEqual(smaller_outer.r.rank(rows+raw), 32)
        self.assertEqual(smaller_outer.spectrum()[32], 588)
        self.assertTrue(self.manifest['symbolic_circuit_verified'])

    def test_inner_selection(self):
        bare = json.loads((ROOT/'workstreams/bare_bch_rm2sub/generated/MANIFEST.json').read_text())
        self.assertEqual(self.manifest['inner'], bare['t128_s19'])
        path = ROOT/'workstreams/bch_rm2sub_bridge/generated/larger_state_inputs_v1/t128_s19_selection.json'
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), self.manifest['inner']['selection_sha256'])
        columns = [int(c,16) for c in json.loads(path.read_text())['selected']['B_columns_hex']]
        self.assertEqual(columns, self.manifest['inner']['columns'])


if __name__ == '__main__':
    unittest.main()

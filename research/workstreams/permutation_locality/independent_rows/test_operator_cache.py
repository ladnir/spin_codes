"""Exact endpoint round-trip and corruption rejection for local memoization."""
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from flint import arb, arb_mat, ctx
from operator_cache import cache_path, load, save, parameters, source_digest


class CacheTests(unittest.TestCase):
    def test_options_and_schema_invalidate_old_keys(self):
        args = SimpleNamespace(groups=40,precision=192,full_feedback=6,
                               window_histogram=8,joint_cancellation=True)
        original = parameters(args,'.052','.95','source')
        self.assertEqual(original['schema'],2)
        for name,value in (('column_density',True),('feedback_density',6),('weight_tilt','9/8')):
            changed = SimpleNamespace(**vars(args),**{name:value})
            self.assertNotEqual(cache_path('memo',original),
                                cache_path('memo',parameters(changed,'.052','.95','source')))

    def test_new_sibling_sources_are_hashed(self):
        original = source_digest()
        read = Path.read_bytes
        for filename in ('universal_density.py','column_density.py','feedback_density.py'):
            def changed(path):
                return read(path)+(b'\n# changed' if path.name == filename else b'')
            with patch.object(Path,'read_bytes',changed):
                self.assertNotEqual(source_digest(),original)

    def test_round_trip_and_rejection(self):
        ctx.prec = 192
        key = dict(groups=1, precision=192, test=True)
        matrix = arb_mat([[arb(i+j).exp() if i != j else arb(0)
                           for j in range(11)] for i in range(11)])
        with tempfile.TemporaryDirectory() as directory:
            self.assertIsNone(load(directory, key))
            save(directory, key, [matrix, matrix])
            result = load(directory, key)
            for actual in result:
                for i in range(11):
                    for j in range(11):
                        self.assertTrue(actual[i,j] == matrix[i,j].upper())
                        self.assertTrue(actual[i,j].is_exact())
            path = cache_path(directory, key)
            record = json.loads(path.read_bytes())
            record['body']['matrices'][0][0][0] = [1, 0]
            path.write_text(json.dumps(record))
            with self.assertRaises(ValueError):
                load(directory, key)


if __name__ == '__main__':
    unittest.main()

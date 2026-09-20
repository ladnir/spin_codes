import hashlib
import json
from pathlib import Path
import unittest

import general_occupancies as g


class NoConstantTests(unittest.TestCase):
    def test_complete_map_audit(self):
        path=Path(__file__).resolve().parent/'NO_CONSTANT_MAP.json'
        record=json.loads(path.read_text())
        for name,digest in record['source_sha256'].items():
            self.assertEqual(hashlib.sha256((g.tv.ROOT/name).read_bytes()).hexdigest(),digest,name)
        columns=record['columns'];rows=g.tv.fixed.maps.generators(columns,19)
        exact=g.tv.fixed.maps.spectrum(rows)
        self.assertEqual(exact,{int(w):int(n) for w,n in record['spectrum'].items()})
        numpy_weights=g.wm.all_weights(columns,19)
        self.assertEqual(int(sum(numpy_weights==0)),1)
        self.assertEqual(int(max(numpy_weights)),80)
        self.assertEqual(min(w for w in exact if w),48)
        self.assertEqual(len(set(columns)),128)
        self.assertNotIn(0,columns)
        self.assertFalse(any((a&b).bit_count()&1 for a in rows for b in rows))
        kernel=g.tv.fixed.maps.dual_spectrum(exact,128,19)
        self.assertEqual(kernel,{int(w):int(n) for w,n in record['kernel'].items()})
        self.assertEqual(min(w for w in kernel if w),5)


if __name__=='__main__':unittest.main()

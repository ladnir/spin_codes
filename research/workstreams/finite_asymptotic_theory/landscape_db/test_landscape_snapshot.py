import tempfile
import unittest
from pathlib import Path

import landscape_snapshot as snapshot


class SnapshotTest(unittest.TestCase):
    def test_reproducible_roundtrip_and_mismatch_rejection(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory=Path(temporary);source=directory/'landscape_export.csv'
            source.write_bytes(b'a,b\r\n1,2\r\n'*1000)
            archive=directory/(source.name+'.gz')
            record=snapshot.compress(source,archive)
            original=archive.read_bytes()
            self.assertEqual(snapshot.compress(source,archive),record)
            self.assertEqual(archive.read_bytes(),original)
            source.unlink();snapshot.restore_one(directory,record)
            self.assertEqual(snapshot.sha(source),record['plain_sha256'])
            source.write_text('keep this file')
            with self.assertRaisesRegex(ValueError,'already exists'):
                snapshot.restore_one(directory,record)
            self.assertEqual(source.read_text(),'keep this file')
            archive.write_bytes(original+b'changed')
            with self.assertRaisesRegex(ValueError,'compressed snapshot hash'):
                snapshot.restore_one(directory,record)


if __name__=='__main__':unittest.main()

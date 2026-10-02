import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import proof_archive as archive


class ProofArchiveTests(unittest.TestCase):
    def make_archive(self, root, name='research/example.py', payload=b'original', recorded=b'original'):
        path = root / 'proof.zip'
        manifest = dict(schema='spin-proof-archive-1', files={name: dict(
            bytes=len(recorded), sha256=hashlib.sha256(recorded).hexdigest())})
        with zipfile.ZipFile(path, 'w') as output:
            output.writestr(name, payload)
            output.writestr('MANIFEST.json', json.dumps(manifest))
        return path

    def test_intact_archive(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.make_archive(Path(directory))
            self.assertEqual(len(archive.verify(path)['files']), 1)

    def test_changed_member(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.make_archive(Path(directory), payload=b'modified')
            with self.assertRaises(ValueError):
                archive.verify(path)

    def test_unsafe_member(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.make_archive(Path(directory), name='../outside.py')
            with self.assertRaises(ValueError):
                archive.verify(path)

    def test_no_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.make_archive(Path(directory))
            original = path.read_bytes()
            with patch.object(archive, 'selected_files', side_effect=AssertionError('must not collect inputs')):
                with self.assertRaises(FileExistsError):
                    archive.create(path)
            self.assertEqual(path.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()

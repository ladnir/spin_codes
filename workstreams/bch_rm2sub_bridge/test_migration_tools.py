"""Small filesystem checks; no numerical sweep or benchmark."""
import json
from pathlib import Path
import tempfile
import unittest

import migrate_legacy_workspace as migration


class MigrationToolsTest(unittest.TestCase):
    def test_copy_preserves_bytes_and_rejects_changed_targets(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, target = root / 'source.py', root / 'nested/target.py'
            source.write_bytes(b'line1\r\nline2\n')
            expected = migration.digest(source)
            migration.copy_checked(source, target, expected)
            self.assertEqual(target.read_bytes(), source.read_bytes())
            migration.copy_checked(source, target, expected)
            target.write_bytes(b'local changes')
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                migration.copy_checked(source, target, expected)
            self.assertEqual(target.read_bytes(), b'local changes')
            with self.assertRaisesRegex(ValueError, 'Source changed'):
                migration.copy_checked(source, root / 'absent', '0' * 64)
            self.assertFalse((root / 'absent').exists())

    def test_manifest_restore_verify_and_path_escape(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            source, target = parent / 'source', parent / 'target'
            source.mkdir()
            (source / 'input.json').write_bytes(b'{"value": 1}\r\n')
            manifest = target / migration.MANIFEST
            manifest.parent.mkdir(parents=True)
            manifest.write_text(json.dumps({'files': [dict(
                path='input.json', sha256=migration.digest(source / 'input.json'))]}))
            migration.restore_or_verify(target, source)
            migration.restore_or_verify(target)
            with self.assertRaisesRegex(ValueError, 'leaves workspace'):
                migration.contained(target.resolve(), '../escape')
            (target / 'input.json').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'changed migration file'):
                migration.restore_or_verify(target)


if __name__ == '__main__':
    unittest.main()

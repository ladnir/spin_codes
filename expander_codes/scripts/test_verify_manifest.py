from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_manifest


class ManifestTest(unittest.TestCase):
    def test_checked_in_manifest(self) -> None:
        data = verify_manifest.load_manifest(verify_manifest.DEFAULT_MANIFEST)
        self.assertEqual(verify_manifest.check_manifest(data, strict_versions=True), [])

    def test_rejects_changed_artifact(self) -> None:
        data = verify_manifest.load_manifest(verify_manifest.DEFAULT_MANIFEST)
        changed = json.loads(json.dumps(data))
        relative = next(iter(changed["artifacts"]))
        changed["artifacts"][relative] = "0" * 64
        errors = verify_manifest.check_manifest(changed, strict_versions=False)
        self.assertTrue(any("artifact digest mismatch" in error for error in errors))

    def test_rejects_path_escape(self) -> None:
        with self.assertRaises(ValueError):
            verify_manifest.resolve_artifact("../outside.json")

    def test_rejects_unpinned_command_input(self) -> None:
        data = verify_manifest.load_manifest(verify_manifest.DEFAULT_MANIFEST)
        changed = json.loads(json.dumps(data))
        relative = next(iter(changed["artifacts"]))
        del changed["artifacts"][relative]
        errors = verify_manifest.check_manifest(changed, strict_versions=False)
        self.assertTrue(any("unpinned artifact" in error for error in errors))


if __name__ == "__main__":
    unittest.main()

"""Fast artifact regression tests; no benchmarks or interval replay."""
import hashlib
from pathlib import Path
import tempfile
import unittest
import shutil
import subprocess
import sys
import re
from unittest.mock import patch

import reproduce


class ArtifactChecks(unittest.TestCase):
    def test_contained_path(self):
        self.assertEqual(reproduce.checked_path('paper/main.tex'), reproduce.ROOT / 'paper/main.tex')

    def test_escape_rejected(self):
        with self.assertRaises(ValueError):
            reproduce.checked_path('../outside')

    def test_manifest_authentication_rejects_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'manifest.json'
            path.write_bytes(b'{}')
            digest = hashlib.sha256(b'{}').hexdigest()
            self.assertEqual(reproduce.authenticated_json(path, digest), {})
            path.write_bytes(b'{"changed": true}')
            with self.assertRaises(ValueError):
                reproduce.authenticated_json(path, digest)

    def test_asymptotic_pins(self):
        reproduce.check_asymptotic()

    def test_archive_requires_complete_evidence(self):
        with patch.object(reproduce, 'dependency_inventory', return_value={'missing': 1, 'mismatched': 0}):
            with self.assertRaises(ValueError):
                reproduce.pack_evidence(Path('unused.zip'))

    def test_archive_never_overwrites(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'existing.zip'
            path.write_bytes(b'preserve me')
            with patch.object(reproduce, 'dependency_inventory', return_value={'missing': 0, 'mismatched': 0, 'files': []}):
                with self.assertRaises(FileExistsError):
                    reproduce.pack_evidence(path)
            self.assertEqual(path.read_bytes(), b'preserve me')

    def test_document_links(self):
        docs = [reproduce.ROOT / 'README.md'] + list((reproduce.ROOT / 'artifact').glob('*.md'))
        for doc in docs:
            for target in re.findall(r'\]\(([^)]+)\)', doc.read_text(encoding='utf-8')):
                if '://' in target or target.startswith('#'):
                    continue
                self.assertTrue((doc.parent / target.split('#')[0]).exists(), f'{doc}: {target}')

    def test_historical_compact_tree_does_not_certify_imt(self):
        # Historical ledgers must not silently stand in for current IMT evidence.
        # Explicit inputs only: no Git metadata or parent-worktree fallback.
        root = reproduce.ROOT
        names = {Path('artifact/reproduce.py'), Path('paper/build_parameter_figures.py'),
                 Path('paper/check_finite_integration.py'),
                 Path('paper/imt_results.py'), Path('paper/build_imt_comparison.py'),
                 Path('paper/build_imt_parameter_figures.py'),
                 Path('paper/build_imt_length_figure.py'),
                 Path('workstreams/transposed_comparison/report.py'),
                 Path('workstreams/transposed_comparison/run.py'),
                 reproduce.BUNDLE / 'SINGLE_SAMPLED_BA_RM2SUB_CERTIFICATE_MANIFEST.json',
                 Path('workstreams/finite_asymptotic_theory/landscape_db/CURRENT_RESULTS_TABLES.md'),
                 Path('workstreams/bare_bch_rm2sub/PERFORMANCE.json'),
                 Path('workstreams/bare_bch_rm2sub/generated/MANIFEST.json')}
        names.update(p.relative_to(root) for p in (root / 'paper').glob('*.tex'))
        names.update(p.relative_to(root) for p in (root / 'paper/figures').glob('*.tex'))
        names.update(name for name, _ in reproduce.asymptotic_entries())
        bridge = Path('workstreams/bch_rm2sub_bridge/generated')
        for m in (16, 18, 20, 22, 24):
            suffix = 'full_split_coverage' if m < 20 else 'ladder_full'
            names.add(bridge / f't128_s19_m{m}_{suffix}_v1.json')
        for suffix in ('selection', 'a_spectrum', 'b_kernel_spectrum'):
            names.add(bridge / f'larger_state_inputs_v1/t128_s19_{suffix}.json')
        with tempfile.TemporaryDirectory(prefix='spin-artifact-check-') as directory:
            destination = Path(directory)
            for name in sorted(names):
                target = destination / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(root / name, target)
            result = subprocess.run([sys.executable, '-B', str(destination / 'artifact/reproduce.py'), 'quick'],
                                    cwd=destination, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Missing IMT evidence', result.stdout + result.stderr)
            self.assertNotIn('QUICK CHECK PASSED', result.stdout)


if __name__ == '__main__':
    unittest.main()

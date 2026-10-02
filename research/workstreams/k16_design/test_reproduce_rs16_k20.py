import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import reproduce_rs16_k20 as recipe


class K20ReproductionTests(unittest.TestCase):
    def test_exact_disjoint_schedule_and_no_bounds(self):
        choices = recipe.prefix_choices()
        self.assertEqual(set(choices), {str(q) for q in range(1, 443)})
        self.assertEqual(choices['1'], {'tilt': '1/2500'})
        self.assertEqual(choices['128'], {'tilt': '23/400'})
        self.assertEqual(choices['442'], {'tilt': '4/25'})
        self.assertTrue(all(set(value) == {'tilt'} for value in choices.values()))
        record = recipe.proposal_record({'test': 'metadata'})
        self.assertTrue(record['proposal_only'])
        self.assertFalse(record['whole_code_certificate'])
        self.assertFalse(record['numerical_endpoints_present'])

    def test_schedule_matches_completed_receipt_when_available(self):
        base = Path(__file__).resolve().parents[3]/'tmp'/'rs-k20-outer256-s20-whole-fresh-p256.json'
        if not base.exists():
            self.skipTest('local completed receipt is intentionally not a repository dependency')
        exact = json.loads(Path(str(base)+'.exact.json').read_text())
        dense = json.loads(Path(str(base)+'.dense.json').read_text())
        self.assertEqual({q: value['tilt'] for q, value in recipe.prefix_choices().items()},
            exact['occupancy_choices'])
        self.assertEqual(list(recipe.DENSE_TILTS), dense['tilts'])
        self.assertEqual(list(recipe.MARKERS), dense['marker_probabilities'])

    def test_dry_run_without_site_packages_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/'new-directory'/'whole.json'
            script = Path(recipe.__file__).resolve()
            result = subprocess.run([sys.executable, '-B', '-S', str(script),
                '--output', str(output), '--dry-run'], check=True, capture_output=True, text=True)
            plan = json.loads(result.stdout)
            self.assertEqual(plan['occupancy_covered'], [1, 4096])
            self.assertEqual(plan['exact_prefix'], [1, 442])
            self.assertEqual(plan['dense_suffix'], [443, 4096])
            self.assertFalse(plan['whole_code_certificate'])
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_each_existing_output_blocks_before_numerical_work(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/'whole.json'
            whole, proposal, components = recipe.output_paths(output)
            for path in (whole, proposal, *components):
                with self.subTest(path=path):
                    path.write_text('preserve me')
                    with self.assertRaises(FileExistsError):
                        recipe.reproduce(output)
                    self.assertEqual(path.read_text(), 'preserve me')
                    path.unlink()
            with self.assertRaises(ValueError):
                recipe.checked_plan(Path(tmp)/'wrong.txt')


if __name__ == '__main__':
    unittest.main()

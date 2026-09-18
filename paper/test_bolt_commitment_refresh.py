"""Check that the x86 refresh changes commitment, not opening accounting."""
import json
import unittest
from pathlib import Path

from build_application_tables import flock_projection, render

ROOT = Path(__file__).resolve().parent


class BoltCommitmentRefreshTests(unittest.TestCase):
    def setUp(self):
        if not (ROOT / "data" / "bolt_x86_commitment.json").is_file():
            self.skipTest("Requires local author benchmark input, outside the core-code artifact")
        def read(name):
            return json.loads((ROOT / "data" / name).read_text(encoding="utf-8"))
        self.data = read("application_results.json")
        self.opening = read("bolt_opening_projection.json")
        self.commit = read("bolt_x86_commitment.json")
        self.ligerito = read("ligerito_standalone.json")

    def test_headline_keeps_opening_and_communication(self):
        total = (self.commit["cases"]["512"]["commit_ms"]
                 + self.opening["case"]["one_opening_amortized_limit_projection_ms"])
        self.assertAlmostEqual(total, 2614.3967585)
        tables = render(self.data, self.opening, self.ligerito, self.commit)
        self.assertIn("$2614^{*}$ & --- & $11.1^{*}$", tables[ROOT / "tables/pcs_standalone.tex"])
        self.assertIn("Calibrated projection & 1273", tables[ROOT / "tables/pcs_opening.tex"])
        self.assertIn("Calibrated projection & 1729", tables[ROOT / "tables/pcs_opening.tex"])

    def test_flock_change_is_exactly_commitment_delta(self):
        for case in self.opening["flock_cases"]:
            new = self.commit["cases"][str(case["input_mib"])]
            old_total = flock_projection(self.data, case, case["commit_measured_ms"])
            new_total = flock_projection(self.data, case, new["commit_ms"])
            self.assertAlmostEqual(old_total - new_total,
                                   case["commit_measured_ms"] - new["commit_ms"])
            self.assertEqual(round(new_total), {32: 228, 128: 983}[case["input_mib"]])

    def test_receipt_scope(self):
        self.assertEqual(self.commit["threads"], 1)
        self.assertEqual(self.commit["hash"], "sha256")
        self.assertEqual(set(self.commit["cases"]), {"32", "128", "512"})
        for case in self.commit["cases"].values():
            self.assertEqual(len(case["receipt_sha256"]), 64)
            self.assertGreater(case["commit_ms"], case["expander_ms"])


if __name__ == "__main__":
    unittest.main()

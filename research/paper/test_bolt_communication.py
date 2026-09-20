"""Small exhaustive checks and regression tests for the explicit Bolt scenario."""
import itertools
import json
import unittest
from pathlib import Path

from bolt_communication import (communication_estimate, expected_siblings,
                                outer_payload, queries, sibling_count)


class BoltCommunicationTests(unittest.TestCase):
    def test_expected_siblings_exhaustively(self):
        for n in (1, 2, 4, 8):
            for q in range(1, n + 1):
                counts = [sibling_count(xs, n.bit_length() - 1)
                          for xs in itertools.combinations(range(n), q)]
                self.assertAlmostEqual(expected_siblings(n, q),
                                       sum(counts) / len(counts), places=10)

    def test_odd_tree_join(self):
        for xs in itertools.combinations(range(8), 3):
            for ys in ((0,), (1,), (0, 1)):
                self.assertEqual(sibling_count(list(xs) + [8+y for y in ys], 4),
                                 sibling_count(xs, 3) + sibling_count(ys, 1) + 2)

    def test_query_counts(self):
        self.assertEqual(queries(.013, 3), 15962)
        self.assertEqual(queries(.5, 2), 241)

    def test_outer_regression(self):
        self.assertAlmostEqual(outer_payload()["expected_outer_bytes"] / 2**20,
                               10.512014092602715, places=9)

    def test_full_scenario_retains_inner_and_messages(self):
        result = communication_estimate()
        self.assertAlmostEqual(result["estimated_bytes"],
                               sum(result["components_bytes"].values()))
        self.assertEqual(f'{result["estimated_mib"]:.1f}', "11.1")
        self.assertGreater(result["estimated_bytes"] - result["expected_outer_bytes"],
                           600000)
        self.assertIn("not a bound", result["scope"])

    def test_table_marks_estimates_and_omits_unmeasured_verification(self):
        from build_application_tables import render
        root = Path(__file__).resolve().parent
        if not (root / "data" / "bolt_x86_commitment.json").is_file():
            self.skipTest("Requires local author benchmark input, outside the core-code artifact")
        data, bolt, ligerito = [json.loads((root / "data" / name).read_text())
                               for name in ("application_results.json",
                                            "bolt_opening_projection.json",
                                            "ligerito_standalone.json")]
        bolt_commit = json.loads((root / "data" / "bolt_x86_commitment.json").read_text())
        table = render(data, bolt, ligerito, bolt_commit)[root / "tables" / "pcs_standalone.tex"]
        self.assertIn(r"Bolt-max (amortized)$^{*}$ & $2614^{*}$ & --- & $11.1^{*}$", table)
        self.assertIn("Bolt communication model LF SHA-256", table)


if __name__ == "__main__":
    unittest.main()

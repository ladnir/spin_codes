"""Retained-data checks only; never starts a benchmark."""
import copy
import json
import unittest
import report


class ComparisonReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(report.RESULTS.read_text())

    def test_complete(self):
        report.validate(self.data)

    def test_generated_table(self):
        self.assertEqual(report.table(self.data),
                         (report.ROOT / "paper/figures/transposed_comparison.tex").read_text())

    def test_partial_run_rejected(self):
        data = copy.deepcopy(self.data)
        data["status"] = "incomplete"
        with self.assertRaises(ValueError):
            report.validate(data)

    def test_raa_reference_scope(self):
        content = report.table(self.data)
        self.assertIn(r"(0.19,{\approx}13)^{*}", content)
        self.assertNotIn(r"(\delta,\lambda)_{20}", content)
        prose = (report.ROOT / "paper/implementation.tex").read_text()
        self.assertIn(r"RAA reference bound at $K=2^{22}$", prose)
        self.assertIn(r"$41.5$ bits at this size", prose)

    def test_changed_sample_rejected(self):
        data = copy.deepcopy(self.data)
        data["runs"][0]["samples_ms"] = [1.0] * 101
        with self.assertRaises(ValueError):
            report.validate(data)

    def test_wrong_rate_rejected(self):
        data = copy.deepcopy(self.data)
        data["runs"][0]["n"] += 1
        with self.assertRaises(ValueError):
            report.validate(data)

    def test_input_reset_rejected(self):
        data = copy.deepcopy(self.data)
        data["runs"][0]["input_policy"] = "reset"
        with self.assertRaises(ValueError):
            report.validate(data)

    def test_wrong_ec_rejected(self):
        data = copy.deepcopy(self.data)
        data["exconv_parameters"]["random_taps"] = 25
        with self.assertRaises(ValueError):
            report.validate(data)


if __name__ == "__main__":
    unittest.main()

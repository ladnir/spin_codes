"""Fast exact checks; no benchmark or large enumeration."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
import reconstruct as r


class QuarterBchTests(unittest.TestCase):
    def test_construction(self):
        c = r.construction()
        self.assertEqual((c["length"], c["dimension"], c["distance_lower_bound"]),
                         (256, 64, 60))
        self.assertEqual(c["nonzero_coset_orbit"], 255)

    def test_retained_generator(self):
        receipt = json.loads((r.HERE / "CONSTRUCTION_AUDIT.json").read_text())
        self.assertEqual(receipt["construction"], r.construction())
        self.assertEqual(receipt["status"], "construction_only")
        self.assertNotIn("spectrum", receipt)

    def test_published_small_table(self):
        a = r.read_spectrum(r.HERE / "sources/EBCH256_63.wd")
        info = r.audit_spectrum(a, 63, 64)
        self.assertEqual(info["minimum_distance"], 64)
        self.assertEqual(a[64], 43180)

    def test_toy_cosets(self):
        # Repetition [4,1,4] inside even-weight [4,3,2], quotient orbit size 3.
        small = [1, 0, 0, 0, 1]
        parent = [1, 0, 6, 0, 1]
        selected = r.reconstruct(small, parent, cosets=3)
        self.assertEqual(selected, [1, 0, 2, 0, 1])
        r.audit_spectrum(selected, 2, 2)
        self.assertEqual(r.macwilliams(parent, 3), small)

    def test_published_parent_and_selected(self):
        small = r.read_spectrum(r.HERE / "sources/EBCH256_63.wd")
        parent = r.read_spectrum(r.HERE / "sources/EBCH256_71.wd")
        selected = r.reconstruct(small, parent)
        self.assertEqual(r.audit_spectrum(parent, 71, 60)["minimum_distance"], 62)
        self.assertEqual(r.audit_spectrum(selected, 64, 60)["minimum_distance"], 62)
        self.assertEqual(parent[62], 130560)
        self.assertEqual(selected[62], 512)
        self.assertEqual(selected[64], 67372)
        self.assertEqual(selected[128], 919000829635375910)
        self.assertEqual(sum(selected), 1 << 64)
        self.assertEqual(selected, r.read_spectrum(r.HERE / "BCH256_64.wd"))

    def test_retained_spectrum_receipt(self):
        receipt = json.loads((r.HERE / "SPECTRUM_AUDIT.json").read_text())
        small_path = r.HERE / "sources/EBCH256_63.wd"
        parent_path = r.HERE / "sources/EBCH256_71.wd"
        small, parent = r.read_spectrum(small_path), r.read_spectrum(parent_path)
        selected = r.reconstruct(small, parent)
        self.assertEqual(receipt["construction"], r.construction())
        self.assertEqual(receipt["status"], "reconstructed_from_supplied_parent_table")
        self.assertEqual(receipt["spectrum"], {str(w): str(v) for w, v in enumerate(selected) if v})
        self.assertEqual(receipt["source63_audit"], r.audit_spectrum(small, 63, 64))
        self.assertEqual(receipt["parent_audit"], r.audit_spectrum(parent, 71, 60))
        self.assertEqual(receipt["selected_audit"], r.audit_spectrum(selected, 64, 60))
        self.assertEqual(receipt["source63_sha256"], hashlib.sha256(small_path.read_bytes()).hexdigest())
        self.assertEqual(receipt["parent_sha256"], hashlib.sha256(parent_path.read_bytes()).hexdigest())
        provenance = json.loads((r.HERE / "sources/PAPER_TABLE_PROVENANCE.json").read_text())
        self.assertEqual(provenance["source"], receipt["parent_source"])
        self.assertEqual(provenance["coefficients_sha256"], receipt["parent_sha256"])
        self.assertEqual((provenance["table"], provenance["pdf_page"]), (7, 6))
        self.assertFalse(provenance["pdf_redistributed"])

    def test_coset_divisibility(self):
        with self.assertRaisesRegex(ValueError, "divisible"):
            r.reconstruct([1, 0, 1], [1, 1, 1])

    def test_negative_coset(self):
        with self.assertRaisesRegex(ValueError, "negative"):
            r.reconstruct([1, 255, 1], [1, 0, 1])

    def test_bad_mass(self):
        with self.assertRaisesRegex(ValueError, "mass"):
            r.audit_spectrum([1, 0, 3, 0, 1], 2, 2)

    def test_bad_macwilliams(self):
        with self.assertRaisesRegex(ValueError, "MacWilliams"):
            r.macwilliams([1, 0, 1, 0, 2], 2)

    def test_duplicate_weight(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.wd"
            path.write_text("0 1\n0 1\n")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                r.read_spectrum(path)


if __name__ == "__main__":
    unittest.main()

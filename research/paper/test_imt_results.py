"""Fast rejection tests for the paper's selected-evidence adapter."""
from fractions import Fraction
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import imt_results as e


class EvidenceTests(unittest.TestCase):
    def test_exact_dyadic_union(self):
        parts = [dict(mantissa='1', exponent=-43), dict(mantissa='3', exponent=-44)]
        total = dict(numerator='5', denominator=str(2**44))
        self.assertAlmostEqual(e.exact_union(parts, total, 40), 44-e.math.log2(5))

    def test_bad_union_and_target(self):
        for total, target in ((dict(numerator='6', denominator=str(2**44)), 40),
                              (dict(numerator='5', denominator=str(2**44)), 42)):
            with self.assertRaises(ValueError):
                e.exact_union([dict(mantissa='5', exponent=-44)], total, target)

    def test_negative_component_rejected(self):
        with self.assertRaises(ValueError):
            e.exact_union([dict(mantissa='-1', exponent=-50)],
                          dict(mantissa='-1', exponent=-50), 40)

    def test_outward_rounded_union(self):
        parts = [dict(mantissa='3', exponent=-50)]
        ceiling = dict(mantissa='4', exponent=-50)
        self.assertEqual(e.exact_union(parts, ceiling, 40, rounded=True), 48)
        with self.assertRaises(ValueError):
            e.exact_union(parts, ceiling, 40)
        with self.assertRaises(ValueError):
            e.exact_union(parts, dict(mantissa='2', exponent=-50), 40, rounded=True)

    def test_large_dyadic_span(self):
        parts = [dict(mantissa='1', exponent=-100000), dict(mantissa='1', exponent=-50)]
        total = dict(mantissa=str((1 << 9950)+1), exponent=-10000)
        # Avoid decimal-string conversion of huge test fixtures.
        with self.assertRaises(ValueError):
            e.exact_union(parts, total, 40)
        self.assertEqual(e.fraction(parts[1]), Fraction(1, 2**50))

    def test_missing_and_changed_receipts(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            with patch.object(e, 'BASE', base):
                with self.assertRaisesRegex(ValueError, 'Missing IMT evidence'):
                    e.authenticate(('absent.json', '0'*64), {})
                # No generated data fixture is needed for the missing-file check.
                with patch.object(Path, 'is_file', return_value=True), \
                     patch.object(Path, 'read_bytes', return_value=b'{}'):
                    with self.assertRaisesRegex(ValueError, 'Receipt changed'):
                        e.authenticate(('changed.json', '0'*64), {})


if __name__ == '__main__':
    unittest.main()

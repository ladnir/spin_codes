"""Independent verifier checks, including rejection of material receipt edits."""
from copy import deepcopy
from fractions import Fraction
import unittest

import verify_fractional as verify


class ArithmeticTests(unittest.TestCase):
    def test_exact_dyadic_and_rounding(self):
        self.assertEqual(verify.dyadic([3, -4]), Fraction(3, 16))
        self.assertEqual(verify.dyadic([3, 4]), 48)
        self.assertEqual(verify.upward_round(Fraction(13, 16), 3), Fraction(7, 8))
        self.assertEqual(verify.upward_round(Fraction(7, 16), 3), Fraction(7, 16))
        with self.assertRaises(ValueError):
            verify.dyadic([True, -1])

    def test_large_union_without_decimal_conversion(self):
        numerator = (1 << 20000) + 17
        record = dict(numerator_hex=hex(numerator), exponent=-20100)
        self.assertEqual(verify.exact_union(record), Fraction(numerator, 1 << 20100))
        with self.assertRaises(ValueError):
            verify.exact_union(dict(numerator_hex='0x01', exponent=-1))


@unittest.skipUnless(verify.DEFAULT_RECEIPT.exists(), 'local ignored final receipt not present')
class ReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = verify.load_json(verify.DEFAULT_RECEIPT)

    def test_current_complete_receipt(self):
        result = verify.audit_record(self.original)
        self.assertTrue(result['target_met'])
        self.assertEqual(result['authenticated_occupancies'], 2048)
        self.assertFalse(result['new_fractional_endpoints_independently_recomputed'])

    def test_reject_changed_endpoint(self):
        modified = deepcopy(self.original)
        modified['selected_occupancy_bounds']['89'][0] += 1
        with self.assertRaisesRegex(ValueError, 'selection'):
            verify.audit_record(modified)

    def test_reject_changed_exact_union(self):
        modified = deepcopy(self.original)
        encoded = modified['exact_union']['numerator_hex']
        modified['exact_union']['numerator_hex'] = hex(int(encoded, 16) + 1)
        with self.assertRaisesRegex(ValueError, 'exact union'):
            verify.audit_record(modified)

    def test_reject_missing_new_occupancy(self):
        modified = deepcopy(self.original)
        del modified['new_occupancy_bounds']['62']
        with self.assertRaisesRegex(ValueError, 'new occupancy coverage'):
            verify.audit_record(modified)

    def test_reject_missing_proof_source_pin(self):
        modified = deepcopy(self.original)
        del modified['source_sha256'][str(verify.HERE / 'SLACK_REVIEW.md')]
        with self.assertRaisesRegex(ValueError, 'required receipt or source pin'):
            verify.audit_record(modified)

    def test_reject_changed_probability_scope(self):
        modified = deepcopy(self.original)
        modified['group_union_outside_power'] = False
        with self.assertRaisesRegex(ValueError, 'scope'):
            verify.audit_record(modified)


if __name__ == '__main__':
    unittest.main()

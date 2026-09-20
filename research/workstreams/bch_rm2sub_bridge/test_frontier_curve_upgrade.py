"""Curve upgrades must preserve old points and never overwrite a certificate."""
import unittest
from frontier_curve_upgrade import upgrade


class CurveUpgradeTest(unittest.TestCase):
    def test_upgrade_preserves_uncertified_evidence(self):
        rows=[dict(message_exponent=26,full_certificate_margin_bits=44.5),
              dict(message_exponent=28,full_certificate_margin_bits=None,modeled_q1_bits=62.5)]
        result=upgrade(rows,28,42.5,62.6)
        self.assertEqual(result[0],rows[0])
        self.assertEqual(result[1]['prior_uncertified_point'],rows[1])
        self.assertIsNone(rows[1]['full_certificate_margin_bits'])
        self.assertIsNone(result[1]['modeled_full_margin_bits'])
        with self.assertRaises(AssertionError):upgrade(result,28,42.6,62.6)
        with self.assertRaises(AssertionError):upgrade(rows,30,40.5,60.6)


if __name__=='__main__':unittest.main()

import unittest
from fractions import Fraction as F
import certify_wider24 as c


class CertificateTests(unittest.TestCase):
    def test_dyadic_roundtrip(self):
        for x in (F(0),F(1),F(7,32),F(1593,2**219),F(1,2**200)):
            self.assertEqual(c.read_dyadic(c.dyadic_record(x)),x)
        with self.assertRaises(ValueError):c.dyadic_record(F(1,3))

    def test_exact_union_not_display_rounding(self):
        below={q:F(1,2**48) for q in range(1,257)}
        total,passed=c.exact_union(below)
        self.assertEqual(total,F(1,2**40));self.assertTrue(passed)
        below[7]+=F(1,2**200)
        self.assertFalse(c.exact_union(below)[1])
        del below[2]
        with self.assertRaises(ValueError):c.exact_union(below)

    def test_caps_only_weaken(self):
        for x in (F(0),F(1,2**400),F(1,2**100),F(3,2)):
            self.assertGreaterEqual(c.round_probability(x),min(F(1),x))


if __name__=='__main__':unittest.main()

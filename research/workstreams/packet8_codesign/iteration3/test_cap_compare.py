from fractions import Fraction
import unittest

import cap_compare as cc


class CapCompareTests(unittest.TestCase):
    def test_cartesian_order(self):
        got=cc.candidates([Fraction(1,2)],[Fraction(0),Fraction(1)],[1,3])
        self.assertEqual(len(got),4)
        self.assertEqual(got[0],(Fraction(1,2),(Fraction(0),Fraction(0))))
        self.assertEqual(got[-1],(Fraction(1,2),(Fraction(1),Fraction(1))))

    def test_invalid_choices(self):
        for caps in ([],[1,1],[0],[9]):
            with self.assertRaises(ValueError):
                cc.candidates([Fraction(1,2)],[Fraction(1)],caps)


if __name__=='__main__':
    unittest.main()

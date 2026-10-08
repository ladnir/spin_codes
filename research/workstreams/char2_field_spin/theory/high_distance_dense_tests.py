from fractions import Fraction as F
import math
import unittest

from gf8_exact_spectrum_dense import excess_bracket
from high_distance_dense import log_membership_row, ball_bits_per_symbol, margin, search


class DenseTests(unittest.TestCase):
    def test_exact_rational_row_agreement(self):
        for bits in (8,16,32,128):
            for alpha in (F(1,2),F(4,5),F(1)):
                for tilt in (F(1),F(3,2),F(4)):
                    exact=excess_bracket(alpha,alpha,tilt,16,4,1<<bits)
                    expected=math.log(exact.numerator)-math.log(exact.denominator)
                    actual=log_membership_row(float(alpha),math.log(float(tilt)),16,4,bits)
                    self.assertAlmostEqual(actual,expected,places=10)

    def test_ball_and_threshold_direction(self):
        self.assertAlmostEqual(ball_bits_per_symbol(F(1,2),16),
                               1+.5*math.log2(65535),places=12)
        self.assertGreater(margin(.95,0,distance=F(2,5)),
                           margin(.95,0,distance=F(41,100)))

    def test_completion_and_endpoint(self):
        result=search(1.,16,F(2,5))
        self.assertEqual(result['output_n'],2**21)
        self.assertEqual(result['input_n'],2**22)
        self.assertEqual(result['rate'],'1/2')
        self.assertFalse(result['whole_code'])
        self.assertTrue(math.isfinite(result['margin_bits']))
        with self.assertRaises(ValueError): log_membership_row(0,0,64,16,16)


if __name__=='__main__': unittest.main()

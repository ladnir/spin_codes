import unittest
from fractions import Fraction as F
from close_larger_state_gap import dyadic_power, as_fraction


class GapDyadicTests(unittest.TestCase):
    def test_upper_rounding_and_cap(self):
        for exponent in range(-150,151):
            for coefficient in (F(1),F(3,2),F(7,4),F(17,16)):
                value = coefficient*F(2)**exponent
                upper = as_fraction(dyadic_power(value))
                self.assertGreaterEqual(upper,value)
                self.assertGreaterEqual(upper,F(1,1<<80))
                if value >= F(1,1<<80):self.assertLess(upper,2*value)
                else:self.assertEqual(upper,F(1,1<<80))


if __name__ == '__main__':unittest.main()

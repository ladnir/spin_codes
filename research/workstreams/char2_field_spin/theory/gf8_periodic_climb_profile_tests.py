from fractions import Fraction as F
import unittest
from gf8_periodic_climb_profile import metadata


class PeriodicProfileTests(unittest.TestCase):
    def test_rates(self):
        for d,pattern,rate in ((84,(6,6,6,7),F(21,25)),(74,(5,6),F(37,44)),(84,(7,),F(3,4))):
            g=metadata(d,pattern)
            self.assertEqual(F(g['nominal_rate']),rate)
            self.assertEqual(F(g['padded_k'],g['output_n']),rate)
            self.assertEqual(g['rows']%(8*len(pattern)),0)
            self.assertLess(10*g['output_threshold'],g['output_n'])
            self.assertEqual(g['kept_outputs_per_phase'],[list(range(x)) for x in pattern])

    def test_invalid(self):
        for pattern in ((),(0,5),(9,),tuple([6]*3)):
            with self.assertRaises(ValueError):metadata(84,pattern)


if __name__=='__main__':unittest.main()

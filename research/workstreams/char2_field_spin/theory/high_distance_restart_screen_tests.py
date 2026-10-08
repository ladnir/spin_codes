from fractions import Fraction as F
import unittest
from gf8_punctured_screen import profile,sparse_bound,middle_bound
from high_distance_restart_screen import adjusted_bound,gv_distance
from higher_distance_screen import output_threshold


class RestartTests(unittest.TestCase):
    def test_reference_agreement(self):
        self.assertEqual(adjusted_bound((5.,2.),8,distance=F(1,10)),
                         sparse_bound(8,(5.,2.),4,4,4,64,32))
        point=(-2.,-.3,2.,0.)
        self.assertEqual(adjusted_bound(point,100000,False,F(1,10)),
                         middle_bound(point,100000,4,4,4,64,32))

    def test_distance_monotonicity_and_threshold(self):
        lo=adjusted_bound((5.,2.),8,distance=F(1,5))
        hi=adjusted_bound((5.,2.),8,distance=F(3,10))
        self.assertGreater(hi,lo)
        self.assertEqual(output_threshold(100,F(3,10)),29)
        self.assertEqual(output_threshold(101,F(3,10)),30)

    def test_gv_benchmark(self):
        self.assertTrue(.37<gv_distance(8)<.39)
        self.assertTrue(.49<gv_distance(128)<.5)
        self.assertTrue(all(gv_distance(a)<gv_distance(b) for a,b in zip((8,16,32,64),(16,32,64,128))))


if __name__=='__main__':
    unittest.main()

import unittest
from fractions import Fraction
from itertools import product

from flint import ctx

from hd_shifted_box_normalizer import scalar_envelope,point_scalar,histogram_affines
from hd_shifted_general_count import ProfileBox,joint_is_feasible,matching_histograms


class BoxNormalizerTests(unittest.TestCase):
    markers = ((Fraction(1,3),Fraction(2,7),Fraction(3,5)),
               (Fraction(2,5),Fraction(1,4),Fraction(1,2)))
    options = dict(repetitions=2,message_symbols_per_group=2,weight_limit=5,q=16)

    def setUp(self):
        # Endpoint extraction itself uses ctx.prec. Keep the 512-bit direct
        # reference from being rounded to the default 53 bits in assertions.
        self.precision = ctx.workprec(512)
        self.precision.__enter__()

    def tearDown(self):
        self.precision.__exit__(None,None,None)

    def test_affines_match_oriented_histograms(self):
        for profile in product(range(5),repeat=4):
            if joint_is_feasible(4,*profile):
                actual = tuple(tuple(offset+sum(c*x for c,x in zip(coeff,profile))
                                     for offset,coeff in row) for row in histogram_affines(4))
                self.assertEqual(actual,matching_histograms(4,*profile))

    def test_singleton_box_recovers_scalar(self):
        for profile in ((1,0,0,0),(2,2,1,1),(4,4,4,4)):
            box = ProfileBox(*((n,n) for n in profile))
            upper = scalar_envelope(4,box,Fraction(3,5),Fraction(1,3),self.markers,**self.options)
            direct = point_scalar(4,profile,Fraction(3,5),Fraction(1,3),self.markers,precision=512,**self.options)
            self.assertLess(abs(float(upper)-float(direct)),1e-12)
            self.assertTrue(upper>=direct.upper())

    def test_boxes_cover_every_feasible_integer_point(self):
        boxes = (ProfileBox(*((0,4),)*4), ProfileBox((1,3),(0,2),(0,2),(0,2)),
                 ProfileBox((2,4),(2,4),(1,3),(1,4)))
        for box in boxes:
            upper = scalar_envelope(4,box,Fraction(3,5),Fraction(1,3),self.markers,**self.options)
            for profile in product(*(range(lo,hi+1) for lo,hi in box.ranges)):
                if joint_is_feasible(4,*profile):
                    direct = point_scalar(4,profile,Fraction(3,5),Fraction(1,3),self.markers,precision=512,**self.options)
                    self.assertTrue(upper>=direct.upper(),(box,profile,upper,direct))

    def test_bad_parameters_are_rejected(self):
        for box in (ProfileBox((4,3),(0,0),(0,0),(0,0)),
                    ProfileBox((0,5),(0,0),(0,0),(0,0))):
            with self.assertRaises(ValueError):
                scalar_envelope(4,box,Fraction(3,5),Fraction(1,3),self.markers)


if __name__=='__main__':
    unittest.main()

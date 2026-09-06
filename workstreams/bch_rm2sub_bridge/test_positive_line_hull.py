import itertools
import unittest
from fractions import Fraction as F
import positive_line_hull as hull


class HullTests(unittest.TestCase):
    def test_all_small_sets_and_ratios(self):
        choices=((0.,0.),(0.,1.),(1.,0.),(1.,1.),(.5,.75),(.75,.5))
        for pairs in itertools.combinations_with_replacement(choices,4):
            keep=hull.indices([p[0] for p in pairs],[p[1] for p in pairs])
            exact=[tuple(F.from_float(x) for x in p) for p in pairs]
            for x in (F(0),F(1,100),F(1,3),F(1),F(3),F(100)):
                self.assertEqual(max(a+b*x for a,b in exact),max(exact[i][0]+exact[i][1]*x for i in keep))


if __name__=='__main__':unittest.main()

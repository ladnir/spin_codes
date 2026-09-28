from fractions import Fraction as Q
import unittest

from flint import arb,ctx

from cover import outward_occupancy_masses


class OccupancyPolynomialTests(unittest.TestCase):
    def test_grouped_polynomial_dominates_exact_recurrence(self):
        ctx.prec=192
        for numbers in ([],[1,8,1,3],[7]*288,[1]*180+[8]*140):
            exact=[Q(1)]
            for number in numbers:
                p=Q(number,10);following=[Q(0)]*(len(exact)+1)
                for i,mass in enumerate(exact):
                    following[i]+=mass*(1-p);following[i+1]+=mass*p
                exact=following
            upper=outward_occupancy_masses(numbers,10)
            self.assertEqual(len(upper),len(exact))
            for x,y in zip(upper,exact):
                self.assertTrue(x>=arb(y.numerator)/y.denominator)
            self.assertTrue(sum(upper,arb(0))<1+arb(2)**-120)


if __name__=='__main__':unittest.main()

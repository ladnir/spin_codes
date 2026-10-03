import unittest
from fractions import Fraction
import numpy as np
import outward_positive as op


class PositiveTests(unittest.TestCase):
    def test_runtime(self): op.check_runtime()

    def test_exact_dot_enclosures_and_support(self):
        values=[0.,1.,np.nextafter(1.,0.),np.finfo(float).tiny,op.ETA,2.**-600,2.**-100]
        rng=np.random.default_rng(917)
        for _ in range(25):
            a=rng.choice(values,size=(3,10));b=rng.choice(values,size=(10,4))
            bound=op.matmul_upper(a,b)
            for r in range(3):
                for c in range(4):
                    exact=sum(Fraction(float(a[r,k]))*Fraction(float(b[k,c])) for k in range(10))
                    self.assertGreaterEqual(Fraction(float(bound[r,c])),exact)
                    if not exact:self.assertEqual(bound[r,c],0.)

    def test_root_integer_witnesses(self):
        for alpha in (Fraction(2,5),Fraction(7,20),Fraction(7,10),Fraction(1,2)):
            for value in (Fraction(1,3),Fraction(1,2),Fraction(17,19),Fraction(32)):
                upper=op.root_power_upper(value,alpha.numerator,alpha.denominator)
                self.assertGreaterEqual(upper**alpha.denominator,value**alpha.numerator)

    def test_fractional_enclosure_extreme_exponents(self):
        x=np.array([0.,op.ETA,2.**-1000,2.**-600,.001,.3,.5,.999,1.,4.])
        for alpha in (Fraction(2,5),Fraction(7,20),Fraction(7,10),Fraction(1,2),Fraction(1)):
            upper=op.fractional_power_upper(x,alpha)
            for value,bound in zip(x,upper):
                self.assertGreaterEqual(Fraction(float(bound))**alpha.denominator,
                    Fraction(float(value))**alpha.numerator)
            self.assertEqual(upper[0],0.)


if __name__=='__main__':unittest.main()

import math
import unittest
from fractions import Fraction as F
from collections import Counter
from christoffel_caps import kraw


class BivariateTests(unittest.TestCase):
    def test_even_code_orthogonality_and_point_caps(self):
        n=8;code=[x for x in range(1<<n) if x.bit_count()%2==0]
        # The even code has OA strength 7, enough for degree-2 squares.
        degree=2;indices=[(a,b) for a in range(3) for b in range(3-a)]
        for c in (2,4,6):
            counts=Counter(((x&((1<<c)-1)).bit_count(),(x>>c).bit_count()) for x in code)
            for a,b in indices:
                for d,e in indices:
                    moment=sum(count*kraw(c,a,i)*kraw(n-c,b,j)*kraw(c,d,i)*kraw(n-c,e,j)
                        for (i,j),count in counts.items())
                    expected=len(code)*math.comb(c,a)*math.comb(n-c,b) if (a,b)==(d,e) else 0
                    self.assertEqual(moment,expected)
            for (i,j),count in counts.items():
                kernel=sum((F(kraw(c,a,i)**2*kraw(n-c,b,j)**2,math.comb(c,a)*math.comb(n-c,b))
                    for a,b in indices),F(0))
                self.assertLessEqual(count,F(len(code))/kernel)


if __name__=='__main__':unittest.main()

from fractions import Fraction as Q
from itertools import product
from math import comb
import unittest

from row_mixture import envelope,pair_components


class RowMixtureTests(unittest.TestCase):
    def test_row_and_nonzero_pair_measure_domination(self):
        n=4;caps=[1,2,3,2,1]
        rows=envelope(caps,4,Q(1,4));pairs=pair_components(rows)
        self.assertEqual(len(pairs),15)
        # Before random lane swaps, test every ordered pair of concrete words.
        for x,y in product(range(1<<n),repeat=2):
            actual=Q(caps[x.bit_count()],comb(n,x.bit_count()))*Q(caps[y.bit_count()],comb(n,y.bit_count()))
            if x==y==0:actual-=1
            bound=Q(0)
            for i,(c,p) in enumerate(rows):
                for j,(d,r) in enumerate(rows):
                    if i==j==0:continue
                    bound+=c*d*p**x.bit_count()*(1-p)**(n-x.bit_count())*r**y.bit_count()*(1-r)**(n-y.bit_count())
            self.assertLessEqual(actual,bound)
        # The pair-component coefficients have the correct location factor.
        self.assertEqual(sum(c for _,c,*_ in pairs),sum(c for c,_ in rows)**2)
        for _,_,p0,p1,p2,_ in pairs:self.assertEqual(p0+p1+p2,1)


if __name__=='__main__':unittest.main()

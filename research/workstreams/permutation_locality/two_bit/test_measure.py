from fractions import Fraction as Q
from itertools import product
from math import comb
import unittest

from flint import arb,ctx

import model
from measure import DensityFolds
from profiles import pair_profiles


class MeasureTests(unittest.TestCase):
    def test_density_dominates_whole_support_sum(self):
        ctx.prec=192;n=4;p=Q(2,5)
        actual=[Q(0),Q(3),Q(8),Q(4),Q(1)]
        shells=[a+Q(1,3) for a in actual]
        cdf=[sum(actual[:u+1])+Q(1,7) for u in range(n+1)]
        fold=DensityFolds(cdf,shells,n=n)
        for lo in range(n+1):
            for hi in range(lo,n+1):
                # A non-monotone, position-sensitive inner cost, not just weight.
                cost=lambda x:Q((x*7)%13+1,17)
                exact=sum((actual[x.bit_count()]/comb(n,x.bit_count())*cost(x)
                           for x in range(1<<n) if lo<=x.bit_count()<=hi),Q(0))
                reference=sum((p**x.bit_count()*(1-p)**(n-x.bit_count())*cost(x) for x in range(1<<n)),Q(0))
                bound=fold.outward(lo,hi,arb(2)/5)*(arb(reference.numerator)/reference.denominator)
                self.assertTrue(bound>=arb(exact.numerator)/exact.denominator)

    def test_pair_profile_measure_pointwise(self):
        n=4;p,r=Q(2,5),Q(3,7);spectrum=[1,2,3,2,1]
        profiles=pair_profiles(spectrum)
        probabilities=[1-p,p*(1-r)/2,p*(1-r)/2,p*r]
        # Shuffled input measure, averaged over independent lane swaps.
        for packets in product(range(4),repeat=n):
            u=sum(x!=0 for x in packets);b=sum(x==3 for x in packets);a=u-b
            numerator=profiles.get((u,b),Q(0))
            actual=numerator/(comb(n,u)*comb(u,b)*2**a)
            reference=Q(1)
            for x in packets:reference*=probabilities[x]
            ratio=numerator/(comb(n,u)*comb(u,b)*(1-p)**(n-u)*(p*(1-r))**a*(p*r)**b)
            self.assertEqual(actual,ratio*reference)


if __name__=='__main__':unittest.main()

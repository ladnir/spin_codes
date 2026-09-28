from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import comb
import unittest

import model
from profiles import pair_profiles,weighted_shells
from support import weighted_union_shells


class ProfileTests(unittest.TestCase):
    def test_direct_enumeration(self):
        for n in range(1,5):
            for spectrum in ([1]+[int(w%2==0) for w in range(1,n+1)],
                             [comb(n,w) for w in range(n+1)],list(range(1,n+2))):
                expected=defaultdict(Q)
                for x,y in product(range(1<<n),repeat=2):
                    mass=Q(spectrum[x.bit_count()],comb(n,x.bit_count()))*Q(spectrum[y.bit_count()],comb(n,y.bit_count()))
                    if mass:expected[(x|y).bit_count(),(x&y).bit_count()]+=mass
                actual=pair_profiles(spectrum)
                self.assertEqual(actual,expected)
                for weight in (Q(1),Q(10,9),Q(2)):
                    self.assertEqual(weighted_shells(actual,n,weight),weighted_union_shells(spectrum,2,weight))

    def test_componentwise_caps(self):
        for spectrum in product(range(3),repeat=4):
            caps=[a+int(i%2==0) for i,a in enumerate(spectrum)]
            actual,bounds=pair_profiles(spectrum),pair_profiles(caps)
            for key,count in actual.items():self.assertGreaterEqual(bounds[key],count)

    def test_invalid_inputs(self):
        for spectrum in ([],[1,-1],[1,Q(1,2)]):
            with self.assertRaises(ValueError):pair_profiles(spectrum)


if __name__=='__main__':unittest.main()

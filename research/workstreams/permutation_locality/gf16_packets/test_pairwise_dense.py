import unittest
from fractions import Fraction as Q
from math import comb
from itertools import product
from pairwise_dense import pair_majorant, combine_pair_mixtures,probe_cell
from shared_mixture import verify


class PairwiseDenseTests(unittest.TestCase):
    def test_exact_probe_scope(self):
        self.assertEqual(probe_cell('.032'),(Q(4,125),Q(4,125)))
        self.assertEqual(probe_cell('.032','1/65536'),(Q(4,125)-Q(1,65536),Q(4,125)+Q(1,65536)))
        for value,radius in ((0,'1/65536'),(1,'1/65536'),('.032','-1/65536')):
            with self.assertRaises(ValueError):probe_cell(value,radius)

    def test_central_plus_residual_exact_domination(self):
        caps=[0,3,7,9,5]
        centers=[Q(1,4),Q(1,2),Q(3,4),Q(1)]
        for exponent in range(7):
            mixture=pair_majorant(caps,centers,zero_bits=None,central_bits=exponent)
            verify(caps,mixture)
        with self.assertRaises(ValueError):
            pair_majorant(caps,centers,zero_bits=10,central_bits=0)

    def test_union_of_pair_mixtures_matches_support_measure(self):
        mixture=[(Q(3),Q(1,4)),(Q(7,2),Q(2,3))]
        result=combine_pair_mixtures(mixture)
        n=3
        def mass(support):
            u=support.bit_count()
            return int(support==0)+sum(c*p**u*(1-p)**(n-u) for c,p in mixture)
        expected=[Q(0)]*(n+1)
        for a,b in product(range(1<<n),repeat=2):
            expected[(a|b).bit_count()] += mass(a)*mass(b)
        expected[0]-=1  # Only the genuine zero-pair product is inactive.
        actual=[comb(n,u)*sum(c*p**u*(1-p)**(n-u) for c,p in result) for u in range(n+1)]
        self.assertEqual(actual,expected)
        self.assertGreater(actual[0],0)

    def test_positive_centers_cover_full_support_without_endpoint(self):
        caps=[0,3,7,9,5]
        mixture=pair_majorant(caps,[Q(1,4),Q(1,2),Q(3,4)],zero_bits=None)
        verify(caps,mixture)
        self.assertTrue(all(p < 1 for _,p in mixture))
        self.assertTrue(all(p < 1 for _,p in combine_pair_mixtures(mixture)))


if __name__ == '__main__':
    unittest.main()

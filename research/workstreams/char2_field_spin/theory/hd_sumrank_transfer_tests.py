from fractions import Fraction as F
from math import comb
import unittest

from hd_randomized_outer_joint_input import split_enumerator
from hd_randomized_outer_proof import prescribed_mds_count
from hd_sumrank_transfer import profiles, minimum_live_weight


class SumRankMomentTests(unittest.TestCase):
    def test_unweighted_normalization(self):
        q = 65536
        mass = [0]*65
        for profile,j,multiplicity,k,t,v in profiles(F(1)):
            mass[j] += multiplicity
            self.assertEqual(k,F(1,q**min(4,j)))
            self.assertEqual(t,1)
            self.assertEqual(v,1)
        self.assertEqual(mass,[comb(64,j) for j in range(65)])

    def test_input_mark_against_mds_kernel(self):
        q, mark = 65536,F(3,2)
        for profile,j,multiplicity,k,t,v in profiles(F(1), input_mark=mark):
            dimension=max(0,j-4)
            expected=sum(comb(j,w)*prescribed_mds_count(w,j,dimension,q)*mark**w
                         for w in range(j+1))/q**j
            self.assertEqual(k,expected)
            self.assertEqual(t,((1+(q-1)*mark)/q)**j)
            self.assertEqual(v,t)

    def test_free_mark_against_independent_split_enumerator(self):
        q,z,input_mark=65536,F(1,4),F(3,2)
        local=[]
        for a in range(9):
            local.append(sum(count*input_mark**i*z**w
                             for i,row in enumerate(split_enumerator(a))
                             for w,count in enumerate(row))/q**a)
        for index,(profile,j,multiplicity,k,t,v) in enumerate(profiles(z,input_mark=input_mark)):
            expected=F(1)
            for a in profile:
                expected*=local[a]
            self.assertEqual(t,expected)

    def test_live_support_floor(self):
        self.assertEqual(minimum_live_weight((0,)*8),61)
        self.assertEqual(minimum_live_weight((1,)*8),53)
        self.assertEqual(minimum_live_weight((4,)*8),29)
        self.assertEqual(minimum_live_weight((8,)*8),0)
        self.assertEqual(minimum_live_weight((0,)*8,8),53)
        for profile,j,multiplicity,k,t,v,live in profiles(F(1,4),with_live=True):
            self.assertGreaterEqual(live,v)
            self.assertLessEqual(v,F(1,4)**minimum_live_weight(profile))
            self.assertLessEqual(live,t)

    def test_quiet_live_mask_cap(self):
        for state,tiles in ((4,1),(4,2),(8,2)):
            for profile,j,multiplicity,k,t,v,live,gamma in profiles(F(1,4),state=state,tiles=tiles,with_live=True,with_gamma=True):
                self.assertGreaterEqual(gamma,0)
                self.assertLessEqual(gamma,k)
                self.assertLessEqual(gamma,live)
                if not j:
                    self.assertEqual(gamma,live)
                    self.assertGreaterEqual(gamma,v)
        for profile,j,multiplicity,k,t,v,live,gamma in profiles(F(1),tiles=2,with_live=True,with_gamma=True):
            self.assertEqual(gamma,k)

    def test_natural_width_normalization(self):
        for tiles,state in ((1,4),(2,4),(4,8)):
            mass=[0]*(8*tiles+1)
            for profile,j,multiplicity,k,t,v in profiles(F(1),state=state,tiles=tiles):
                mass[j]+=multiplicity
                self.assertEqual(k,F(1,65536**min(state,sum(min(4,a) for a in profile))))
                self.assertEqual(t,1)
                self.assertEqual(v,1)
            self.assertEqual(mass,[comb(8*tiles,j) for j in range(8*tiles+1)])


if __name__=='__main__':
    unittest.main()

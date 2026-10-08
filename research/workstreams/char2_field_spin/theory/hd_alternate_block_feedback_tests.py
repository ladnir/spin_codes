"""Regression checks for pairing block4 emission with block8 feedback."""
from functools import lru_cache
from itertools import permutations
from math import comb
from fractions import Fraction as F
import random
import unittest
import numpy as np
from hd_block_emission_screen import paired_feedback_caps
from hd_bucket_feedback import maps
from hd_blockwise_shared_transfer import mm


class AlternateBlockFeedbackTests(unittest.TestCase):
    def test_exact_packet_entries_and_group_counts(self):
        from hd_rank_two_block_packets import matrix,grouped_packets
        from hd_block_emission_screen import exact_moments
        z=F(3,4);q=65536;d=q**16-1
        for profile in ((0,0,0,0),(0,0,1,2),(2,3,4,8),(4,4,4,4),(16,16,16,16)):
            k,t,v,gamma=exact_moments(profile,z,16,64,4,q,True);j=sum(profile)
            expected=(k,t-k,v/d if j else F(0),(v if j else (1-F(1,q**2))*v)+gamma/q**2)
            self.assertEqual(matrix(profile,z,4),expected)
            for feedback in (4,8):self.assertTrue(all(x>=0 for x in matrix(profile,z,feedback)))
        for feedback in (4,8):
            for cutoff in (4,8):
                groups=grouped_packets(z,feedback,maximum=8,cutoff=cutoff)
                for j in range(9):self.assertEqual(sum(n for k,n,_ in groups if k==j),comb(64,j))

    def test_larger_feedback_annihilates_block_emission(self):
        _,c4=maps(local_rank=4);a=tuple(zip(*c4))
        for global_nodes in (False,True):
            _,c8=maps(local_rank=8,global_nodes=global_nodes)
            self.assertFalse(any(x for row in mm(c8,a,65536) for x in row))

    def test_worst_pairing_caps_all_labels(self):
        q=65536;z=F(3,4);rng=random.Random(4816)
        @lru_cache(None)
        def quiet(j):return sum(F(comb(j,c),q**(c+min(8,j-c)))*z**(j-c)*(1-z)**c for c in range(j+1))
        for _ in range(100):
            row=tuple(sorted(rng.randrange(17) for _ in range(4)))
            cap,rho,lower=paired_feedback_caps(np.array([row]),float(z),4,16)
            for a,b,c,d in set(permutations(row)):
                self.assertGreaterEqual(cap[0],float(quiet(a+b)*quiet(c+d))*(1-1e-12))
                self.assertLessEqual(rho[0],min(8,a+b)+min(8,c+d))
                self.assertLessEqual(lower[0],float(quiet(a+b)*quiet(c+d))*(1+1e-12))

    def test_exact_matrix_dominates_every_pairing_entrywise(self):
        from hd_rank_two_block_packets import matrix,local_data
        q=65536;d=q**16-1;z=F(1,8);rng=random.Random(51639)
        _,k8,flat,powers=local_data(z)
        # In the first case the pairings have different quiet moments.
        # Using max(K) in T-K would underestimate the birth transition.
        rows=[(0,0,8,8)]+[tuple(sorted(rng.randrange(17) for _ in range(4))) for _ in range(100)]
        self.assertNotEqual(k8[0]*k8[16],k8[8]**2)
        for row in rows:
            upper=matrix(row,z,8);j=sum(row);t=powers[j]
            for a,b,c,e in set(permutations(row)):
                quiet=k8[a+b]*k8[c+e];rho=min(8,a+b)+min(8,c+e)
                v=t*(q**16*flat[16-a]*flat[16-b]*flat[16-c]*flat[16-e]-1)/d
                gamma=min(quiet,z**max(0,13-max(row))/q**rho,d*v)
                exact=(quiet,t-quiet,v/d if j else F(0),
                       (v if j else (1-F(1,q**2))*v)+gamma/q**2)
                self.assertTrue(all(x>=y for x,y in zip(upper,exact)))


if __name__=='__main__':unittest.main()

"""Exact toy/local tests for the full-RS candidate; no distance gate."""

from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest

import numpy as np

import hd_field128_full_rs_screen as full
import hd_field128_lift_screen as old
from hd_field128_lift_certificate import good_geometry
from hd_blockwise_maps import maps,matmul,matvec,power


def toy_moments(a,positions,z=F(1,2),q=16):
    c = tuple(zip(*a))
    quiet,returns = defaultdict(int),defaultdict(int)
    for values in product(range(q),repeat=len(positions)):
        x = [0]*len(a)
        for index,value in zip(positions,values):
            x[index] = value
        syndrome = matvec(c,x,q)
        if not any(syndrome):
            quiet[sum(bool(v) for v in x)] += 1
        else:
            emission = matvec(a,syndrome,q)
            returns[sum(bool(u^v) for u,v in zip(x,emission))] += 1
    moment = lambda counts:sum((count*z**weight for weight,count in counts.items()),F(0))/q**len(positions)
    return moment(quiet),moment(returns)


class FullRSTests(unittest.TestCase):
    def test_fixed_full_map_ca_and_additive_power_sums(self):
        a,c = maps(width=64,state=16,base=128,q=256)
        self.assertFalse(any(x for row in matmul(c,a,256) for x in row))
        for exponent in range(63):
            value = 0
            for point in range(64):
                value ^= power(point,exponent,256)
            self.assertEqual(value,0)

    def test_toy_one_block_versus_four_blocks_and_actual_return(self):
        a = tuple(tuple(power(point,exponent,16) for exponent in range(4)) for point in range(16))
        four = tuple(tuple(int(point//4==column) for column in range(4)) for point in range(16))
        self.assertFalse(any(x for row in matmul(tuple(zip(*a)),a,16) for x in row))
        self.assertFalse(any(x for row in matmul(tuple(zip(*four)),four,16) for x in row))
        for positions in ((),(0,),(0,1),(0,4),(1,7)):
            quiet,actual = toy_moments(a,positions)
            k,birth,ret,live = full.exact_moments(len(positions),F(1,2),16,16,4,True)
            self.assertEqual(quiet,k)
            self.assertGreaterEqual(ret,actual)
            self.assertGreaterEqual(birth,actual)
            self.assertGreaterEqual((16**4-1)*live,actual)
        one_k,one_j = toy_moments(a,(0,1))
        four_k,four_j = toy_moments(four,(0,1))
        self.assertLess(one_k,four_k)
        self.assertLess(one_j,four_j)
        self.assertEqual(one_k,F(1,16**2))
        self.assertEqual(four_k,F(1+15*F(1,2)**2,16**2))

    def test_positive_expansions_equal_exact_subtractions(self):
        q,n,h,z = 16,16,4,F(7,8)
        t = z+(1-z)/q
        for j in range(n+1):
            k,birth,ret,live = full.exact_moments(j,z,q,n,h,True)
            self.assertEqual(birth,t**j-k)
            inactive = n-j
            flat = sum((F(comb(inactive,e),q**min(h,e))*z**(inactive-e)*(1-z)**e
                        for e in range(inactive+1)),F(0))
            self.assertEqual(live,t**j*(q**h*flat-1)/(q**h-1))
            self.assertLessEqual(ret,full.exact_moments(j,z,q,n,h,False)[2])

    def test_log_operator_matches_exact_toy_and_selected_full_entries(self):
        for q,n,h,positions in ((16,16,4,range(17)),(256,64,16,(0,1,16,33,64))):
            for joint in (False,True):
                logs = full.local_logs(.5,q,n,h,joint).reshape((n+1,4))
                for j in positions:
                    exact = full.exact_packet(j,F(1,2),q,n,h,joint)
                    np.testing.assert_allclose(logs[j],[old.log_fraction(x) for x in exact],rtol=0,atol=2e-11)

    def test_q128_tiny_mark_has_no_positive_term_clipped_to_zero(self):
        z,q = F(1,1<<200),1<<128
        for joint in (False,True):
            logs = full.local_logs(float(z),q,joint=joint)[1].ravel()
            exact = full.exact_packet(1,z,q,joint=joint)
            np.testing.assert_allclose(logs,[old.log_fraction(x) for x in exact],rtol=0,atol=5e-10)
            self.assertTrue(np.all(np.isfinite(logs)))
            self.assertLess(logs[2],-2048)

    def test_collapsed_geometry_preserves_every_fair_swap_count(self):
        counts = defaultdict(F)
        for s,f,count,occ in old.sparse_geometry(6):
            counts[s,f] += F(count,2**s)
        self.assertEqual(dict(counts),{(s,f):F(n) for s,f,n in full.geometry(6)})
        counts = defaultdict(int)
        for s,f,count,occ in good_geometry(128):
            counts[s,f] += count
        self.assertEqual(dict(counts),{(s,f):n for s,f,n in full.geometry(128,good=True)})
        self.assertEqual(len(full.geometry(128,good=True)),65)

    def test_simultaneous_basis_change_does_not_change_return_map(self):
        a = tuple(tuple(power(point,exponent,16) for exponent in range(4)) for point in range(16))
        c = tuple(zip(*a))
        h = tuple(tuple(int((row^1)==column) for column in range(4)) for row in range(4))
        self.assertEqual(matmul(matmul(a,h,16),matmul(h,c,16),16),matmul(a,c,16))
        self.assertNotEqual(matmul(matmul(a,h,16),c,16),matmul(a,c,16))
        self.assertFalse(any(x for row in matmul(c,matmul(a,h,16),16) for x in row))


if __name__=='__main__':
    unittest.main()

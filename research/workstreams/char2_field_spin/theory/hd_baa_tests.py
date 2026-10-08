from fractions import Fraction as F
from itertools import combinations,product
import math
import unittest

import numpy as np
from scipy.special import logsumexp

from hd_baa_screen import (accumulator_iowe,accumulator_probability,band_gate,
                          distant_repeat_moment,double_prefix_transfer,
                          group_ceiling,local_transfer)
from hd_baa_moment import row_log_transfer


class BaaTests(unittest.TestCase):
    def test_iowe_exhaustive_small_fields(self):
        for q,n in ((2,7),(3,5),(4,5)):
            actual=[[0]*(n+1) for _ in range(n+1)]
            for word in product(range(q),repeat=n):
                state=0
                output_weight=0
                for x in word:
                    state=state^x if q in (2,4) else (state+x)%q
                    output_weight += state != 0
                actual[sum(x!=0 for x in word)][output_weight]+=1
            for w in range(n+1):
                for h in range(n+1):
                    self.assertEqual(accumulator_iowe(n,w,h,q),actual[w][h])

    def test_large_q_normalization(self):
        for q in (256,65536):
            for w in (0,1,2,3,8,12):
                self.assertEqual(sum(accumulator_probability(12,w,h,q)
                                     for h in range(13)),F(1))

    def test_distant_repeat_boundary_gluing(self):
        for q in (2,4,256):
            y=.37
            for k in range(1,5):
                direct=0.
                transitions=local_transfer(y,q)
                for mask in product(range(2),repeat=k):
                    v=np.array([1.,0.])
                    for active in mask+mask:
                        v=v@transitions[active]
                    direct+=(q-1)**sum(mask)*float(v.sum())
                self.assertAlmostEqual(distant_repeat_moment(k,y,q),direct,places=5)

    def test_double_prefix_is_stochastic_and_coupled(self):
        for q in (2,4,256,65536):
            for matrix in double_prefix_transfer(1.,q):
                np.testing.assert_allclose(matrix.sum(axis=1),np.ones(4))
            zero,live=double_prefix_transfer(.4,q)
            self.assertEqual(zero[0,0],1)
            self.assertEqual(live[0,3],.4)
            # A live first-prefix state injects the second: not independent scans.
            self.assertEqual(zero[2,3],.4)

    def test_row_transfer_counts_all_mds_messages(self):
        for b,k,q in ((4,2,4),(8,4,256)):
            row=row_log_transfer(b,k,q)
            np.testing.assert_allclose(logsumexp(row,axis=1),k*math.log(q)/2,
                                       rtol=1e-12,atol=1e-12)

    def test_row_transfer_against_explicit_gf4_mds(self):
        def mul(a,b):
            out=0
            while b:
                if b&1:out^=a
                b>>=1
                a<<=1
                if a&4:a^=7
            return out
        transitions=local_transfer(1.,4)
        expected=np.zeros((5,5))
        for r in range(5):
            incoming=tuple(int(i<r) for i in range(4))
            for a,b in product(range(4),repeat=2):
                support=tuple(int((a^mul(b,x))!=0) for x in range(4))
                for outgoing in product(range(2),repeat=4):
                    mass=math.prod(transitions[support[i]][incoming[i],outgoing[i]]
                                   for i in range(4))
                    expected[r,sum(outgoing)]+=mass
        np.testing.assert_allclose(np.exp(row_log_transfer(4,2,4))*4,expected,
                                   rtol=1e-12,atol=1e-12)

    def test_finite_group_ceilings(self):
        expected={16:832265,32:832273,64:832289,128:832321,256:832385}
        for b,d in expected.items():
            result=group_ceiling(b)
            self.assertEqual(result['distance_upper_bound'],d)
            self.assertLess(d,838861)
            g,m,D=result['check_groups'],result['selected_groups'],4
            self.assertGreater(2*math.comb(m-1,D-1),math.comb(g-1,D-1))
            self.assertLessEqual(2*math.comb(m-2,D-1),math.comb(g-1,D-1))
        self.assertEqual(group_ceiling(16,profile=(2,1))['distance_upper_bound'],912849)
        self.assertEqual(group_ceiling(16,profile=(3,1))['distance_upper_bound'],949729)

    def test_band_failures_exceed_forty_bit_target(self):
        for b,a in ((16,5),(32,6),(64,6)):
            self.assertGreater(band_gate(b,a)['failure_probability_lower'],2**-40)
        self.assertLess(band_gate(128,6)['failure_probability_lower'],2**-40)

    def test_disjoint_hypergeometric_second_moment_bound(self):
        # Exact small analogue of negatively correlated container-hit events.
        population=range(8)
        containers=({0,1},{2,3},{4,5})
        counts=[]
        for chosen in combinations(population,4):
            selected=set(chosen)
            counts.append(sum(c<=selected for c in containers))
        total=len(counts)
        p=F(math.comb(6,2),math.comb(8,4))
        claimed=3*p/(1+2*p)
        actual=F(sum(v>0 for v in counts),total)
        self.assertGreaterEqual(actual,claimed)


if __name__=='__main__':unittest.main()

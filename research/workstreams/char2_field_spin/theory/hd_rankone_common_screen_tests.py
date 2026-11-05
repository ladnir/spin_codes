from fractions import Fraction as F
from itertools import product
import unittest

import numpy as np

from hd_rankone_common_screen import (packet_matrix, log_product, log_power,
                                      pair_coefficients, exact_sparse_regions, arbitrary_path_lower_log)


def mul(a,b):
    value=0
    while b:
        if b&1:value^=a
        a=(a<<1)^(7 if a&2 else 0)
        b>>=1
    return value


def dot(a,b):
    value=0
    for x,y in zip(a,b):value^=mul(x,y)
    return value


class RankOneCommonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.states=tuple(product(range(4),repeat=2))
        cls.index={x:i for i,x in enumerate(cls.states)}
        cls.kernel=[[F(0) for _ in cls.states] for _ in cls.states]
        mixers=[(u,v) for u in cls.states[1:] for v in cls.states if dot(u,v)!=1]
        for i,x in enumerate(cls.states):
            for u,v in mixers:
                scalar=dot(v,x)
                y=tuple(a^mul(b,scalar) for a,b in zip(x,u))
                cls.kernel[i][cls.index[y]]+=F(1,len(mixers))

    def test_exact_invertible_refresh_law_including_scalar_lines(self):
        for i,x in enumerate(self.states[1:],1):
            line={tuple(mul(c,a) for a in x) for c in range(1,4)}
            for j,y in enumerate(self.states):
                expected=F(0) if j==0 else F(4,15) if y==x else F(1,15) if y in line else F(1,20)
                self.assertEqual(self.kernel[i][j],expected)
            self.assertEqual(sum(self.kernel[i]),1)
        for j in range(1,16):
            self.assertEqual(sum(self.kernel[i][j] for i in range(1,16)),1)

    def test_three_component_cone_for_every_live_state_on_small_supports(self):
        # Two independent repetition[4,1] blocks: A q repeats each state
        # coordinate four times and C sums each block, so C A=0.
        z=F(2,3)
        for support in ((),(0,),(0,1),(0,4),(0,1,4)):
            occ=(sum(i<4 for i in support),sum(i>=4 for i in support))
            measures=[]
            for old in self.states:
                measure=[F(0) for _ in self.states]
                for values in product(range(4),repeat=len(support)):
                    x=[0]*8
                    for i,value in zip(support,values):x[i]=value
                    syndrome=[0,0]
                    for i,value in enumerate(x):syndrome[i//4]^=value
                    before=tuple(a^b for a,b in zip(old,syndrome))
                    weight=sum((value^old[i//4])!=0 for i,value in enumerate(x))
                    scale=z**weight/F(4**len(support))
                    for j,probability in enumerate(self.kernel[self.index[before]]):
                        measure[j]+=scale*probability
                measures.append(np.array([float(value) for value in measure]))
            for mode in ('density','mass'):
                cap=packet_matrix(occ,float(z),q=4,h=1,flat_mode=mode)
                for measure,row in [(measures[0],0)]+[(v,1) for v in measures[1:]]+[(sum(measures[1:])/15,2)]:
                    self.assertLessEqual(measure[0],cap[row,0]+2e-14)
                    excess=np.maximum(0.,measure[1:]-cap[row,2]/15).sum()
                    self.assertLessEqual(excess,cap[row,1]+2e-14)

    def test_empty_packet_has_no_return_and_no_reused_dense_refresh_cap(self):
        z=.5
        cap=packet_matrix((0,0,0,0),z)
        np.testing.assert_array_equal(cap[0],[1.,0.,0.])
        self.assertEqual(cap[1,0],0.)
        self.assertEqual(cap[2,0],0.)
        self.assertEqual(cap[1,1],z**13/256)
        self.assertGreaterEqual(cap[2,2],z**13/256)

    def test_log_products_preserve_noncommuting_order(self):
        a=np.array([[1.,2.,0.],[0.,1.,3.],[.5,0.,1.]])
        b=np.array([[.5,0.,1.],[1.,.25,0.],[0.,2.,1.]])
        with np.errstate(divide='ignore'):
            np.testing.assert_allclose(np.exp2(log_product(np.log2(a),np.log2(b))),a@b)
            np.testing.assert_allclose(np.exp2(log_power(np.log2(a),7)),np.linalg.matrix_power(a,7))

    def test_pair_geometry_normalization_and_complete_sparse_profiles(self):
        # At z=1 the zero row sum is one. Route coefficients therefore
        # count unsigned pair placements before state-path multiplication.
        coeff=pair_coefficients(1.,1.,2)
        from math import comb
        for (s,f),logs in coeff.items():
            self.assertAlmostEqual(np.exp2(logs)[0].sum(),comb(32,s)*comb(32-s,f))
        self.assertEqual(set(exact_sparse_regions(.99,1.,2)),{(0,0),(1,0),(2,0),(2,1)})
        with self.assertRaises(ValueError):exact_sparse_regions(.99,1.,3)

    def test_arbitrary_path_obstructs_medium_weights_for_every_mark(self):
        # The universal statement is analytic: throughout this interval,
        # 1024*r-8*32768-9>0 and256*r+13*(32768-1)-734003<=0.
        for weight in (257,1024,1203):
            self.assertGreater(1024*weight-8*32768-9,0)
            self.assertLessEqual(256*weight+13*(32768-1)-734003,0)
            for z in (.01,.5,.999,1.):
                for alpha in (.001,.1,1.):
                    self.assertGreater(arbitrary_path_lower_log(weight,z,alpha),-5.)
        for occ in ((0,0,0,1),(0,0,3,16),(1,2,3,4),(16,16,16,16)):
            z=.7
            matrix=packet_matrix(occ,z)
            self.assertGreaterEqual(matrix[1,1]+1e-30,z**(sum(occ)+13)/256)
            self.assertGreaterEqual(matrix[0,1]+1e-30,(255/256)*z**sum(occ)/256)


if __name__=='__main__':
    unittest.main()

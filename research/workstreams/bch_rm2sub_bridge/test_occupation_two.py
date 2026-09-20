"""Exact field, epoch-placement, support-pair, and Q1-axis checks for Q2."""
import itertools
import math
import unittest
from collections import Counter
from fractions import Fraction as F

import numpy as np
from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import occupation_two as q2
from test_activation_bridge import field_mul,word


class OccupationTwoTests(unittest.TestCase):
    def test_weight_two_exact_field_paths(self):
        rows=[0xff,0xaa,0xcc,0xf0]
        columns=[sum(((r>>j)&1)<<i for i,r in enumerate(rows)) for j in range(8)]
        self.assertEqual(len(set(columns)),8)
        self.assertTrue(all(columns[j]^columns[k] for j,k in itertools.combinations(range(8),2)))
        spectrum=Counter(word(q,rows).bit_count() for q in range(1,16))
        checks=0
        for z in (F(1,3),F(3,4),F(99,100)):
            exact=[]
            for weight in range(3):
                matrix=[[F(0) for _ in range(16)] for _ in range(16)]
                for positions in itertools.combinations(range(8),weight):
                    x=sum(1<<j for j in positions); syndrome=0
                    for j in positions: syndrome^=columns[j]
                    for q in range(16):
                        factor=z**(x^word(q,rows)).bit_count()/(15*math.comb(8,weight))
                        for alpha in range(1,16):
                            matrix[q][field_mul(alpha,q)^syndrome]+=factor
                exact.append(matrix)
            envelope=q2.epochs(8,4,spectrum,z,F)
            starts=[([F(int(q==j)) for q in range(16)],int(j!=0)) for j in range(16)]
            starts += [([F(0)]+[F(1,15)]*15,2)]
            starts += [([F(0)]+[F(int(q!=omit),14) for q in range(1,16)],2) for omit in range(1,16)]
            for state,kind in starts:
                def visit(v,u,depth):
                    nonlocal checks
                    self.assertLessEqual(sum(v),sum(u)); checks+=1
                    if not depth:return
                    for weight in range(3):
                        nxt=[sum((v[q]*exact[weight][q][j] for q in range(16)),F(0)) for j in range(16)]
                        visit(nxt,q1.positive_vector_mul(u,envelope[weight]),depth-1)
                visit(state,tuple(F(int(j==kind)) for j in range(3)),3)
        self.assertEqual(checks,3840)

    def test_all_region_placements_and_support_pairs(self):
        t,s,spectrum,z,length=8,4,{4:14,8:1},F(3,4),16
        epoch=q2.epochs(t,s,spectrum,z,F); region=q2.regions(t,s,spectrum,z,F,length)
        identity=tuple(F(int(i==j)) for i in range(3) for j in range(3))
        for weight in range(3):
            total=[F(0)]*9
            for selected in itertools.combinations(range(length),weight):
                result=identity
                for start in range(0,length,t):
                    count=sum(start<=position<start+t for position in selected)
                    result=q1.positive_mul(result,epoch[count])
                total=[a+b for a,b in zip(total,result)]
            self.assertEqual(tuple(v/math.comb(length,weight) for v in total),region[weight])
        rounded=np.array([[q2.upper_float(arb(v.numerator)/v.denominator) for v in matrix] for matrix in region])
        upper=q2.pair_counts(rounded,4)
        for a in range(5):
            for b in range(5):
                total=F(0)
                for first in itertools.combinations(range(4),a):
                    for second in itertools.combinations(range(4),b):
                        v=(F(1),F(0),F(0))
                        for j in range(4):
                            v=q1.positive_vector_mul(v,region[int(j in first)+int(j in second)])
                        total+=sum(v)
                total/=math.comb(4,a)*math.comb(4,b)
                self.assertGreaterEqual(F.from_float(float(upper[a,b])),total)
                self.assertLess(float(F.from_float(float(upper[a,b]))-total),1e-12)
        tiny=np.array([np.nextafter(0.,np.inf)])
        self.assertGreater(q2.mul_up(tiny,np.array([0.5]))[0],0.)

    def test_actual_q1_axes_and_saved_aggregation(self):
        ctx.prec=256
        name='t128_s15'; t,s,spectrum=base.load_map(name)
        path=base.HERE/'generated'/f'{name}_q2_j-75.json'; receipt=base.read(path)
        for filename,digest in receipt['local_sha256'].items():
            self.assertEqual(base.sha(base.HERE/filename),digest)
        for filename,digest in receipt['outer_sha256'].items():
            self.assertEqual(base.sha(base.BCH/filename),digest)
        lam=(arb(-75)/10).exp(); z=(-lam).exp()
        region=q2.regions(t,s,spectrum,z,arb)
        first=q1.positive_regions(t,s,spectrum,z,arb,8192//t)
        for matrix,other in zip(region,first):
            self.assertTrue(all(x.overlaps(y) for x,y in zip(matrix,other)))
        axis=q1.positive_coefficients(*first,arb,256)
        with np.load(path.with_suffix('.npz')) as saved:
            conditional=saved['conditional_pair_upper']; rounded=saved['region_upper']
        replay=q2.pair_counts(rounded)
        replay=np.minimum(1.,q2.mul_up(replay,q2.upper_float((209716*lam).exp())))
        np.testing.assert_array_equal(replay,conditional)
        for w,value in enumerate(axis):
            expected=value*(209716*lam).exp()
            if expected>=1:
                self.assertEqual(conditional[0,w],1.)
                self.assertEqual(conditional[w,0],1.)
            else:
                self.assertGreaterEqual(arb(float(conditional[0,w])),expected)
                self.assertGreaterEqual(arb(float(conditional[w,0])),expected)
        total,_=q2.aggregate(conditional)
        self.assertEqual(total,base.decode(receipt['Q2_upper']))
        self.assertLess(total,F(1,1<<91))
        first=base.read(base.HERE/'generated'/'t128_s15_activation_q1_outward.json')
        self.assertLess(total+base.decode(first['Q1_upper']),F(1,1<<49))
        self.assertFalse(receipt['all_occupations_certified'])


if __name__=='__main__':
    unittest.main()

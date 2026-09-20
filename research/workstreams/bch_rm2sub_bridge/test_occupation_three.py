"""Exact Q3 epoch and auxiliary-conditioning identities on small instances."""
import itertools
import math
import unittest
from collections import Counter
from fractions import Fraction as F

from flint import arb,ctx
import activation_bridge as q1
import occupation_three as q3
from test_activation_bridge import field_mul,word
from certify_q3_compact import dyadic_upper
from verify_q3_compact import independent_regions


class Q3Tests(unittest.TestCase):
    def test_full_field_prefixes(self):
        rows=[0xff,0xaa,0xcc,0xf0]
        columns=[sum(((r>>j)&1)<<i for i,r in enumerate(rows)) for j in range(8)]
        spectrum=Counter(word(q,rows).bit_count() for q in range(1,16));checks=0
        for z in (F(1,3),F(3,4),F(99,100)):
            exact=[]
            for j in range(4):
                kernel=[[F(0) for _ in range(16)] for _ in range(16)]
                for positions in itertools.combinations(range(8),j):
                    x=sum(1<<i for i in positions);b=0
                    for i in positions:b^=columns[i]
                    self.assertEqual(b==0,j==0)
                    for q in range(16):
                        factor=z**(x^word(q,rows)).bit_count()/(15*math.comb(8,j))
                        for alpha in range(1,16):kernel[q][field_mul(alpha,q)^b]+=factor
                exact.append(kernel)
            envelope=q3.region_matrices(8,4,spectrum,z,F,length=8)
            starts=[([F(int(q==j)) for q in range(16)],int(j!=0)) for j in range(16)]
            starts += [([F(0)]+[F(1,15)]*15,2)]
            starts += [([F(0)]+[F(int(q!=omit),14) for q in range(1,16)],2) for omit in range(1,16)]
            for state,kind in starts:
                def visit(v,u,depth):
                    nonlocal checks
                    self.assertLessEqual(sum(v),sum(u));checks+=1
                    if depth==0:return
                    for j in range(4):
                        nxt=[sum((v[q]*exact[j][q][k] for q in range(16)),F(0)) for k in range(16)]
                        visit(nxt,q1.positive_vector_mul(u,envelope[j]),depth-1)
                visit(state,tuple(F(int(j==kind)) for j in range(3)),2)
        self.assertEqual(checks,2016)

    def test_region_and_conditioning_identity(self):
        z=F(3,4);spectrum={4:14,8:1}
        epoch=q3.region_matrices(8,4,spectrum,z,F,length=8)
        region=q3.region_matrices(8,4,spectrum,z,F,length=16)
        for degree in range(4):
            total=[F(0)]*9
            for support in itertools.combinations(range(16),degree):
                a=sum(j<8 for j in support)
                result=q1.positive_mul(epoch[a],epoch[degree-a])
                total=[x+y for x,y in zip(total,result)]
            self.assertEqual(tuple(x/math.comb(16,degree) for x in total),region[degree])
        ctx.prec=256
        independent=independent_regions(8,4,spectrum,arb(3)/4,16)
        for exact,matrix in zip(region,independent):
            for value,bound in zip(exact,itertools.chain.from_iterable(matrix)):
                self.assertTrue(bound.contains(arb(value.numerator)/value.denominator))
        totals={weights:F(0) for weights in itertools.product(range(4),repeat=3)}
        supports=[set(j for j in range(3) if (mask>>j)&1) for mask in range(8)]
        for triple in itertools.product(supports,repeat=3):
            state=(F(1),F(0),F(0))
            for j in range(3):
                state=q1.positive_vector_mul(state,region[sum(j in support for support in triple)])
            totals[tuple(map(len,triple))]+=sum(state)
        for weights in totals:totals[weights]/=math.prod(math.comb(3,w) for w in weights)
        for ps in ((F(1,3),F(1,2),F(3,4)),(F(3,4),F(3,4),F(3,4))):
            mixed=q3.coefficient(region,ps,F,positions=3)
            expected=F(0)
            for weights,value in totals.items():
                probability=math.prod(math.comb(3,w)*p**w*(1-p)**(3-w) for w,p in zip(weights,ps))
                expected+=probability*value
                self.assertLessEqual(value,mixed/probability)
            self.assertEqual(expected,mixed)

    def test_dyadic_bounds_and_multiplicities(self):
        for numerator in (1,17,1<<300):
            for denominator in (1,19,1<<500):
                value=F(numerator,denominator);upper=dyadic_upper(value)
                self.assertLessEqual(value,upper)
                self.assertLessEqual(upper,value*(1+F(1,1<<190)))
        for box in itertools.combinations_with_replacement(range(13),3):
            self.assertEqual(q3.multiplicity(box),len(set(itertools.permutations(box))))


if __name__=='__main__':unittest.main()

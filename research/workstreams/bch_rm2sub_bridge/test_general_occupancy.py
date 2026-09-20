"""Exact tests of all epoch weights, including nonzero kernel inputs."""
import itertools
import math
import unittest
from collections import Counter
from fractions import Fraction as F

import general_occupancy as model
import activation_bridge as q1
from test_activation_bridge import field_mul,word


class GeneralTests(unittest.TestCase):
    def test_all_epoch_weights_and_prefixes(self):
        rows=[0xff,0xaa,0xcc,0xf0]
        columns=[sum(((r>>j)&1)<<i for i,r in enumerate(rows)) for j in range(8)]
        spectrum=Counter(word(q,rows).bit_count() for q in range(1,16))
        kernel=Counter()
        for x in range(256):
            syndrome=0
            for j in range(8):
                if x>>j&1:syndrome^=columns[j]
            if syndrome==0:kernel[x.bit_count()]+=1
        self.assertEqual(dict(kernel),{0:1,4:14,8:1})
        checks=0
        for z in (F(1,3),F(3,4),F(99,100)):
            exact=[]
            for j in range(9):
                matrix=[[F(0) for _ in range(16)] for _ in range(16)]
                for selected in itertools.combinations(range(8),j):
                    x=sum(1<<k for k in selected);b=0
                    for k in selected:b^=columns[k]
                    for q in range(16):
                        factor=z**(x^word(q,rows)).bit_count()/(15*math.comb(8,j))
                        for alpha in range(1,16):matrix[q][field_mul(alpha,q)^b]+=factor
                exact.append(matrix)
            envelope=model.epoch_matrices(8,4,spectrum,kernel,z,F,8)
            starts=[([F(int(q==j)) for q in range(16)],int(j!=0)) for j in range(16)]
            starts += [([F(0)]+[F(1,15)]*15,2)]
            starts += [([F(0)]+[F(int(q!=omit),14) for q in range(1,16)],2) for omit in range(1,16)]
            for state,kind in starts:
                def visit(v,u,depth):
                    nonlocal checks
                    self.assertLessEqual(sum(v),sum(u));checks+=1
                    if depth==0:return
                    for j in range(9):
                        nxt=[sum((v[q]*exact[j][q][k] for q in range(16)),F(0)) for k in range(16)]
                        visit(nxt,q1.positive_vector_mul(u,envelope[j]),depth-1)
                visit(state,tuple(F(int(j==kind)) for j in range(3)),2)
        self.assertEqual(checks,8736)

    def test_all_region_weights_and_compositions(self):
        kernel={0:1,4:14,8:1};spectrum={4:14,8:1};z=F(3,4)
        epoch=model.epoch_matrices(8,4,spectrum,kernel,z,F,8)
        region=model.regions(8,4,spectrum,kernel,z,F,16,length=16)
        for j in range(17):
            total=[F(0)]*9
            for a in range(max(0,j-8),min(j,8)+1):
                value=q1.positive_mul(epoch[a],epoch[j-a]);count=math.comb(8,a)*math.comb(8,j-a)
                total=[x+count*y for x,y in zip(total,value)]
            self.assertEqual(tuple(x/math.comb(16,j) for x in total),region[j])
        for q in range(1,17):
            rows=list(model.compositions(q))
            self.assertEqual(len(rows),math.comb(q+4,4))
            self.assertEqual(sum(model.multiplicity(c) for c in rows),5**q)
        counts=(2,1,0,0,1);ps=[F(1,4),F(1,2),F(1,3),F(2,3),F(3,4)]
        dist=model.distribution(counts,ps,F);self.assertEqual(sum(dist),1)
        expanded=[p for count,p in zip(counts,ps) for _ in range(count)]
        expected=[F(0)]*5
        for bits in itertools.product((0,1),repeat=4):
            expected[sum(bits)]+=math.prod(p if b else 1-p for p,b in zip(expanded,bits))
        self.assertEqual(dist,expected)

    def test_density_majorant(self):
        # Unnormalized counting measure, not a normalized random-code law.
        n=4;groups=((1,2),(3,4));caps={1:2,2:3,3:1,4:1};ps=(F(1,3),F(3,4))
        for band,p in zip(groups,ps):
            gamma=max(F(caps[w],math.comb(n,w))/(p**w*(1-p)**(n-w)) for w in band)
            for mask in range(1<<n):
                w=mask.bit_count()
                original=F(caps[w],math.comb(n,w)) if w in band else F(0)
                reference=p**w*(1-p)**(n-w)
                self.assertLessEqual(original,gamma*reference)


if __name__=='__main__':unittest.main()

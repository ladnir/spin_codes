from fractions import Fraction as Q
from itertools import product,permutations
from math import comb,factorial
import unittest

import mixture
import boundary
import density_tangent


class MixtureTests(unittest.TestCase):
    def test_group_expansion_preserves_total_coefficient(self):
        rows=[(Q(1),Q(0)),(Q(1),Q(1)),(Q(7),Q(1,2)),(Q(3),Q(1,3)),(Q(3),Q(2,3))]
        groups=mixture.group_components(rows)
        self.assertEqual(len(groups),70)
        self.assertEqual(sum(c for _,c,_,_ in groups),sum(c for c,_ in rows)**4)
        for name,coefficient,law,active in groups:
            self.assertEqual(sum(law),1)
            self.assertEqual(active,int(name!='0000'))
            self.assertGreater(coefficient,0)

    def test_packet_law_after_uniform_lane_shuffle(self):
        rows=[(Q(1),Q(0)),(Q(1),Q(1)),(Q(1),Q(1,3)),(Q(1),Q(2,5))]
        groups=mixture.group_components(rows)
        for name,_,law,_ in groups:
            p=[rows[int(i)][1] for i in name]
            expected=[Q(0)]*16
            for bits in product((0,1),repeat=4):
                probability=Q(1)
                for b,r in zip(bits,p):probability*=r if b else 1-r
                for order in permutations(range(4)):
                    word=sum(bits[i]<<j for j,i in enumerate(order))
                    expected[word]+=probability/24
            for x,value in enumerate(expected):
                self.assertEqual(value,law[x.bit_count()]/comb(4,x.bit_count()))

    def test_log_power_against_direct_matrix_power(self):
        import numpy as np
        a=np.array([[.7,.2],[.1,.5]])
        for n in (1,2,3,19,64):
            self.assertAlmostEqual(mixture.log_power(a,n),np.log(np.linalg.matrix_power(a,n)[0].sum()),places=12)

    def test_complete_boundary_lists_against_exhaustive_counts(self):
        features=[(Q(0),Q(0)),(Q(1,10),Q(0)),(Q(1,4),Q(1,10)),(Q(0),Q(1,2))]
        active=[0,1,1,1];slots=6;q_min=2
        for lows in product(range(4),repeat=2):
            cell=tuple(v for k in lows for v in (Q(k,8),Q(k+1,8)))
            expected=[]
            for c1,c2,c3 in product(range(slots+1),repeat=3):
                q=c1+c2+c3
                if not q_min<=q<=slots:continue
                counts=(slots-q,c1,c2,c3)
                point=[sum(c*f[j] for c,f in zip(counts,features))/slots for j in range(2)]
                if all(lo<=x<=hi for x,lo,hi in zip(point,cell[::2],cell[1::2])):expected.append(counts)
            actual=boundary.compositions(features,active,cell,slots,q_min,term_limit=100,node_limit=100000,list_limit=1000)
            self.assertEqual(actual,tuple(sorted(expected)))

    def test_boundary_limits_never_return_partial_list(self):
        features=[(Q(0),),(Q(1,10),),(Q(1,5),)]
        result=boundary.compositions(features,[0,1,1],(Q(0),Q(1)),6,1,node_limit=1,term_limit=100)
        self.assertIsNone(result)

    def test_five_category_density_tangents_exhaustively(self):
        n=4
        for anchors in ((4,0,0,0,0),(1,1,1,1,0),(3,1,0,0,0)):
            base,ratios=density_tangent.exact(n,anchors)
            for head in product(range(n+1),repeat=4):
                if sum(head)>n:continue
                counts=(*head,n-sum(head));actual=Q(n**n,factorial(n));upper=base
                for count,ratio in zip(counts,ratios):
                    actual*=Q(factorial(count),count**count);upper*=ratio**count
                self.assertLessEqual(actual,upper)


if __name__=='__main__':unittest.main()

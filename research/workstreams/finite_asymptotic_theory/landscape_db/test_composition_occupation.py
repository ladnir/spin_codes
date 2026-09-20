"""Independent counting checks for the composition-preserving sparse bound."""
import itertools
import math
import unittest
from fractions import Fraction as F

import numpy as np

import activation_occupation as general
import composition_occupation as sparse
from test_activation_q1 import mul,identity


class SparseCompositionTest(unittest.TestCase):
    def test_multiplicities_cover_every_ordered_assignment(self):
        for groups in (1,3,6):
            for q in (1,2,3,4):
                indices=sparse.compositions(groups,q)
                self.assertEqual(len(indices),math.comb(groups+q-1,q))
                self.assertAlmostEqual(float(np.exp(sparse.log_multiplicities(indices)).sum()),groups**q,places=8)

    def test_actual_support_counting_measure_is_dominated(self):
        # One fixed [4,2] code with nonzero spectrum 2*y^2+y^4.
        # Enumerate the independent shuffled supports, retaining the same code.
        counts={2:2,4:1}; q=3; block=4
        regions=[[[F(1+i+2*j+d,32) for j in range(3)] for i in range(3)] for d in range(q+1)]
        supports=[(sum(1<<i for i in bits),F(count,math.comb(block,w)))
                  for w,count in counts.items() for bits in itertools.combinations(range(block),w)]
        actual=F(0)
        for rows in itertools.product(supports,repeat=q):
            product=identity(3)
            for position in range(block):
                degree=sum(bool(word&(1<<position)) for word,mass in rows)
                product=mul(product,regions[degree])
            actual+=math.prod(mass for word,mass in rows)*sum(product[0])
        engine=sparse.SparseComposition(counts,block,q)
        logs=np.log(np.array(regions,dtype=float))
        for shift in (-1.,0.,1.):
            bound=math.exp(engine.aggregate(engine.components(logs,0,.2,shift),q))
            self.assertGreaterEqual(bound+1e-12,float(actual))

    def test_preserving_compositions_improves_adaptive_maximum(self):
        counts={1:3,2:4,4:1}; block=4; q=3
        rng=np.random.default_rng(531)
        regions=np.log(rng.uniform(.005,.7,size=(q+1,3,3)))
        engine=sparse.SparseComposition(counts,block,q)
        roots,active,inactive=general.density_roots(counts,block,engine.bands,.5)
        final=list(general.adaptive_logs(regions,roots,active,inactive))[-1][1]
        adaptive=q*math.log(len(engine.bands))+general.terminal_log(final,block)
        actual=engine.aggregate(engine.components(regions,0,.2,.5),q)
        self.assertLessEqual(actual,adaptive+1e-12)


if __name__=='__main__':
    unittest.main()

"""Exact regressions for the faster graded-basis / quadratic-tower maps."""
import random
import unittest

from hd_blockwise_shared_transfer import mm,power,multiply,dot,rank
from hd_sumrank_maps import construction,mask_ranks,Q


class TowerBasisTests(unittest.TestCase):
    def test_moore_butterfly(self):
        m=construction(16,True)
        self.assertEqual(m['basis'],(1,79,256,20224))
        self.assertEqual(rank(m['moore'][:4],Q),4)
        rng=random.Random(98212)
        for _ in range(500):
            x=[rng.randrange(Q) for _ in range(4)]
            a=x[0]^multiply(79,x[1],Q)
            b=x[2]^multiply(79,x[3],Q)
            c=a^x[1];d=b^x[3]
            y=(a^multiply(256,b,Q),c^multiply(335,d,Q),
               a^multiply(257,b,Q),c^multiply(334,d,Q))
            self.assertEqual(y,tuple(dot(row,x,Q) for row in m['moore'][:4]))

    def test_nilpotence(self):
        m=construction(16,True)
        self.assertEqual(mm(m['C'],m['A'],Q),((0,)*16,)*16)
        self.assertEqual(mm(m['C'],mm(m['F'],m['A'],Q),Q),((0,)*16,)*16)
        self.assertEqual(m['F'],construction(16)['F'])

    def test_native_map_hash(self):
        m=construction(16,True)
        h=0
        for name in ('F','A','C'):
            for row in m[name]:
                for x in row:h=((h^x)*0x100000001b3)&((1<<64)-1)
        # Independently generated native setup checks the same exact maps.
        self.assertEqual(h,14806900258954508744)

    def test_ranks(self):
        rng=random.Random(726182)
        for _ in range(120):
            a=sorted(rng.sample(range(64),rng.randrange(65)))
            b=sorted(rng.sample(range(64),rng.randrange(65)))
            r=mask_ranks(a,b,16,True)
            for kind in ('quiet','flat','input'):
                self.assertEqual(r[kind+'_actual'],r[kind+'_expected'],(a,b,r))


if __name__=='__main__':unittest.main()

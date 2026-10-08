"""Exact map and residual-rank checks for the sixteen-state FFT refinement."""
import random
import unittest
from hd_sumrank_maps import construction,mask_ranks,Q
from hd_blockwise_shared_transfer import mm,power,rank


class FftWeightTests(unittest.TestCase):
    def test_maps(self):
        m=construction(16,True,True)
        self.assertEqual(m['F'],construction(16,True)['F'])
        self.assertEqual(len({power(a,4,Q) for a in m['representatives']}),8)
        self.assertNotIn(0,m['representatives'])
        self.assertEqual(mm(m['C'],m['A'],Q),((0,)*16,)*16)
        self.assertEqual(mm(m['C'],mm(m['F'],m['A'],Q),Q),((0,)*16,)*16)
        self.assertEqual(rank(m['C'],Q),16)
        h=0
        for name in ('F','A','C'):
            for row in m[name]:
                for x in row:h=((h^x)*0x100000001b3)&((1<<64)-1)
        self.assertEqual(h,17847517350586406002)

    def test_masks(self):
        rng=random.Random(5136)
        for _ in range(160):
            a=sorted(rng.sample(range(64),rng.randrange(65)))
            b=sorted(rng.sample(range(64),rng.randrange(65)))
            r=mask_ranks(a,b,16,True,True)
            for kind in ('quiet','flat','input'):
                self.assertEqual(r[kind+'_actual'],r[kind+'_expected'],(a,b,r))


if __name__=='__main__':unittest.main()

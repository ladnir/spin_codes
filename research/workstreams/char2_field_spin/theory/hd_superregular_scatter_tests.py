"""Exact finite-map tests for the proposed, untimed Lagrange-basis variant."""
from itertools import combinations
import random
import unittest
import numpy as np
from hd_common_refresh_enumerator import extend
from hd_blockwise_shared_transfer import mm
from hd_scattered_feedback import count_small_returns
from hd_superregular_scatter import local_map,maps,rank_data,return_caps


def rank(rows):
    basis=()
    for row in rows:basis=extend(basis,row)
    return len(basis)


class SuperregularScatterTests(unittest.TestCase):
    def test_every_local_square_minor(self):
        a=local_map();tested=0
        for size in range(1,5):
            for rows in combinations(range(16),size):
                for cols in combinations(range(4),size):
                    self.assertEqual(rank(tuple(tuple(a[r][c] for c in cols) for r in rows)),size)
                    tested+=1
        self.assertEqual(tested,4844)

    def test_scattered_annihilation(self):
        a,c=maps()
        self.assertFalse(any(v for row in mm(c,a,65536) for v in row))

    def test_rank_cap_on_concrete_restrictions(self):
        a,c=maps();ac=mm(a,c,65536);occ,zero_patterns,bounds,_=rank_data()
        profiles={tuple(row):i for i,row in enumerate(occ)}
        patterns={tuple(row):i for i,row in enumerate(zero_patterns)}
        rng=random.Random(4844)
        for _ in range(600):
            counts=[rng.randrange(17) for _ in range(4)]
            supports=[rng.sample(range(16),n) for n in counts]
            zeros=[rng.sample([i for i in range(16) if i not in u],rng.randrange(17-len(u))) for u in supports]
            u=[16*b+i for b in range(4) for i in supports[b]]
            e=[16*b+i for b in range(4) for i in zeros[b]]
            actual=rank(tuple(tuple(ac[i][j] for j in u) for i in e))
            bound=bounds[profiles[tuple(sorted(counts))],patterns[tuple(min(4,len(v)) for v in zeros)]]
            self.assertLessEqual(bound,actual)

    def test_caps_cover_small_exact_averages(self):
        occ,_,_,_=rank_data()
        for z in (.07,.5,.7):
            caps=return_caps(z)
            for pattern,(counts,denominator) in count_small_returns(*maps()).items():
                exact=sum(count*z**weight for weight,count in counts.items())/denominator
                got=caps[np.all(occ==pattern,axis=1)][0]
                self.assertGreaterEqual(got,exact*(1-1e-12))


if __name__=='__main__':unittest.main()

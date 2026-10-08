"""Independent checks for joint return-rank bounds."""
from fractions import Fraction as F
from itertools import combinations,product
import random
import unittest
import numpy as np
from hd_common_refresh_enumerator import extend
from hd_scattered_feedback import maps,small_return_counts
from hd_scattered_return_rank import grid_rank_lower,rank_data,inactive_caps,averaged_inactive_caps
from hd_blockwise_shared_transfer import mm


def rank(rows):
    basis=()
    for row in rows:basis=extend(basis,row)
    return len(basis)


class ScatteredReturnRankTests(unittest.TestCase):
    def test_degree_order_counterexample(self):
        # An arbitrary RS column order is NOT superregular: degree-three
        # basis has zeros at 0,1. Best-assignment pivoting was invalid.
        a,c=maps();ac=mm(a,c,65536)
        support=(48,);zeros=(0,1)
        self.assertEqual(rank(tuple(tuple(ac[i][j] for j in support) for i in zeros)),0)
        self.assertEqual(grid_rank_lower((0,0,0,1),(2,0,0,0)),0)

    def test_random_concrete_input_and_output_restrictions(self):
        a,c=maps();ac=mm(a,c,65536);rng=random.Random(99035)
        for _ in range(600):
            counts=[rng.randrange(17) for _ in range(4)]
            supports=[rng.sample(range(16),n) for n in counts]
            empty=[rng.sample([i for i in range(16) if i not in u],rng.randrange(17-len(u))) for u in supports]
            u=[16*b+i for b in range(4) for i in supports[b]]
            e=[16*b+i for b in range(4) for i in empty[b]]
            actual=rank(tuple(tuple(ac[i][j] for j in u) for i in e))
            self.assertLessEqual(grid_rank_lower(counts,[min(4,len(v)) for v in empty]),actual)

    def test_caps_majorize_exact_one_and_two_input_averages(self):
        occ,_,_,_=rank_data()
        for z in (.07,.5,.7):
            caps=inactive_caps(z)
            for pattern,(counts,denominator) in small_return_counts().items():
                exact=sum(count*z**weight for weight,count in counts.items())/denominator
                got=caps[np.all(occ==pattern,axis=1)][0]
                self.assertGreaterEqual(got,exact*(1-1e-12))

    def test_full_syndrome_rank_matches_uniform_inactive_moment(self):
        from hd_block_emission_screen import uniform_mark
        occ,_,_,_=rank_data();q=65536.
        for z in (.07,.5,.7):
            caps=inactive_caps(z)
            for profile in ((4,4,4,4),(4,5,6,7),(6,6,6,6)):
                inactive=np.prod([uniform_mark(16-l,4,z) for l in profile])-q**(-16)
                expected=inactive*(z+(1-z)/q)**(sum(profile)-16)
                actual=caps[np.all(occ==profile,axis=1)][0]
                self.assertAlmostEqual(actual/expected,1,places=10)

    def test_labelled_average_refines_worst_bound(self):
        occ,_,_,_=rank_data()
        for z in (.07,.5,.7):
            averaged=averaged_inactive_caps(z);worst=inactive_caps(z)
            self.assertTrue(np.all(averaged<=worst*(1+1e-12)))
            for pattern,(counts,denominator) in small_return_counts().items():
                exact=sum(count*z**weight for weight,count in counts.items())/denominator
                got=averaged[np.all(occ==pattern,axis=1)][0]
                self.assertGreaterEqual(got,exact*(1-1e-12))


if __name__=='__main__':unittest.main()

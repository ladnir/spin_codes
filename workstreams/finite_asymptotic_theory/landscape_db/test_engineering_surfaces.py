"""Check ensemble normalization and the explanatory onset identity."""
import math
import unittest

import bch_growth_model as model
from study_engineering_surfaces import random_mean_counts


class EngineeringSurfaceTest(unittest.TestCase):
    def test_random_mean_mass(self):
        for b in (8,32,128,512):
            self.assertAlmostEqual(sum(random_mean_counts(b).values())/((1<<(b//2))-1),1.,places=13)

    def test_mean_matches_all_two_dimensional_binary_subspaces(self):
        # Every independent pair generates {0,a,b,a^b}; deduplicate bases.
        subspaces={tuple(sorted((a,b,a^b))) for a in range(1,16) for b in range(a+1,16)}
        actual={w:sum(sum(v.bit_count()==w for v in code) for code in subspaces)/len(subspaces)
                for w in range(1,5)}
        expected=random_mean_counts(4)
        for w in actual: self.assertAlmostEqual(actual[w],expected[w],places=14)

    def test_random_onset_closed_form(self):
        for b in (8,32,64,128,256,512):
            counts=random_mean_counts(b);tail=.2*b;m=math.floor(tail);f=tail-m
            literal=sum(n*model.onset_probability(b,w) for w,n in counts.items())
            closed=(((1<<(b//2))-1)/((1<<b)-1))*((1+f)*(1<<m)-1)
            self.assertAlmostEqual(literal/closed,1.,places=12)


if __name__=='__main__': unittest.main()

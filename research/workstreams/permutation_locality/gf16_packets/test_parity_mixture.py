import unittest
from itertools import product
from math import prod
from fractions import Fraction as Q
import parity_mixture as pm


class ParityMixtureTests(unittest.TestCase):
    def test_exact_small_activity_patterns(self):
        for probabilities in ([Q(1,2)]*4,[Q(0),Q(1,3),Q(2,3),Q(1)],
                              [Q(0),Q(0),Q(1,3),Q(2,3)],[Q(0),Q(0),Q(0),Q(1,3)]):
            activity=1-prod(1-p for p in probabilities)
            for n in (1,2,3):
                counts=[Q(0)]*(1<<n)
                for packets in product(range(16),repeat=n):
                    parity=0;support=0;mass=Q(1)
                    for position,value in enumerate(packets):
                        parity^=value;support|=int(value!=0)<<position
                        mass*=prod(p if value>>i&1 else 1-p for i,p in enumerate(probabilities))
                    if not parity:counts[support]+=mass
                for floor in range(1,n+1):
                    factor=pm.parity_factor(probabilities,floor)
                    for support,mass in enumerate(counts):
                        u=support.bit_count()
                        if u>=floor:self.assertLessEqual(mass,factor*activity**u*(1-activity)**(n-u))

    def test_single_row_keeps_its_parity_loss(self):
        self.assertEqual(pm.parity_factor([0,0,0,Q(2,5)],38),1)
        factor=pm.parity_factor([Q(1,2)]*4,38)
        self.assertEqual(factor,(1+15*Q(1,15)**38)/16)
        self.assertEqual(pm.parity_factor([0,0,0,0],38),1)

    def test_inactive_component_is_not_rescaled(self):
        rows=[(Q(1),Q(0)),(Q(1),Q(1)),(Q(100),Q(1,2))]
        original=pm.sc.group_components(rows);refined=pm.components(rows,2)
        for before,after in zip(original,refined):
            self.assertEqual(before[0],after[0]);self.assertEqual(before[2:],after[2:])
            self.assertLessEqual(after[1],before[1])
            if not before[3]:self.assertEqual(before,after)

    def test_reject_invalid_input(self):
        for probabilities,floor in [([0]*3,1),([0,0,0,2],1),([0]*4,0)]:
            with self.assertRaises(ValueError):pm.parity_factor(probabilities,floor)


if __name__=='__main__':unittest.main()

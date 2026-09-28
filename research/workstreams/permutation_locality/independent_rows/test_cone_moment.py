"""Clamping is checked against exact tilted measures, not just loose matrices."""
from fractions import Fraction as Q
from itertools import product
from math import log
import unittest

from flint import arb, arb_mat, ctx
from cone_moment import clamp, moment, log_moment
from occupancy_memory import Z,M,C
from mature_tail import L48,L56


class ConeMoment(unittest.TestCase):
    def test_valid_uppers(self):
        for masses in product((Q(0),Q(1,8),Q(3,4)),repeat=3):
            actual = [Q(0)]*11
            actual[M] = sum(masses)
            actual[L48] = masses[0]
            actual[L56] = masses[0]+masses[1]
            actual[C] = max(masses)
            for slack in product((Q(0),Q(1,8),Q(2)),repeat=4):
                upper = actual[:]
                for i,s in zip((M,C,L48,L56),slack):
                    upper[i] += s
                result = clamp(upper)
                self.assertTrue(all(a<=b<=c for a,b,c in zip(actual,result,upper)))

    def test_concrete_weighted_process(self):
        ctx.prec = 192
        transitions = [[Q(1,8),Q(1,16),Q(1,4)],
                       [Q(1,16),Q(1,8),Q(1,8)],
                       [Q(1,32),Q(1,16),Q(1,8)]]
        t = arb_mat(11,11)
        t[Z,M] = 1
        for target in (C,L48,L56):
            t[Z,target] = 2
        total = max(map(sum,transitions))
        t[M,M] = arb(total.numerator)/total.denominator
        for target in (C,L48,L56):
            t[M,target] = 2*t[M,M]
        actual = [Q(1,3)]*3
        for steps in range(1,12):
            if steps>1:
                actual = [sum(actual[i]*transitions[i][j] for i in range(3)) for j in range(3)]
            upper = moment(t,steps)
            point = upper.fmpq()
            self.assertGreaterEqual(Q(int(point.p),int(point.q)),sum(actual))
            raw = t**steps
            self.assertLessEqual(upper,raw[Z,M])
            proposal = log_moment([[float(t[i,j]) for j in range(11)] for i in range(11)],steps)
            self.assertAlmostEqual(proposal,log(float(upper)),places=12)


if __name__ == '__main__':
    unittest.main()

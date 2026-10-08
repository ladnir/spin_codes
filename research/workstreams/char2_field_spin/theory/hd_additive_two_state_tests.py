"""Exhaustive GF4 checks independent of the derived counting formulas."""
from fractions import Fraction as F
from itertools import product
import math
import unittest

import numpy as np

from hd_additive_two_state import (exact_probability, final_log_moment,
                                   local_transfer, support_law)
from hd_rational_post_tests import mul


def coefficients():
    return [(c,d1,d2,f) for c,d1,d2,f in product(range(1,4),repeat=4)
            if c ^ d1 ^ mul(f,d2)]


def encode(x, setup):
    s1 = s2 = 0
    output = []
    for a,(c,d1,d2,f) in zip(x,setup):
        output.append(s1 ^ mul(f,s2) ^ mul(c,a))
        s1 ^= mul(d1,a)
        s2 ^= mul(d2,a)
    return tuple(output)


class AdditiveTwoStateTests(unittest.TestCase):
    def test_all_gf4_local_values(self):
        tuples = coefficients()
        self.assertEqual(len(tuples),(4-1)**2*(4**2-3*4+3))
        for s1,s2,x in product(range(4),repeat=3):
            counts = {}
            for c,d1,d2,f in tuples:
                u = s1 ^ mul(d1,x)
                v = s2 ^ mul(d2,x)
                y = s1 ^ mul(f,s2) ^ mul(c,x)
                key = 2*int(bool(u))+int(bool(v)), int(bool(y))
                counts[key] = counts.get(key,F(0))+F(1,len(tuples))
            self.assertEqual(counts,support_law(2*int(bool(s1))+int(bool(s2)),
                                               int(bool(x)),4))

    def test_all_length_two_setups_and_messages(self):
        tuples = coefficients()
        counts = [[0]*3 for _ in range(3)]
        for setup in product(tuples,repeat=2):
            outputs = set()
            for x in product(range(4),repeat=2):
                y = encode(x,setup)
                outputs.add(y)
                counts[sum(v != 0 for v in x)][sum(v != 0 for v in y)] += 1
            self.assertEqual(len(outputs),16)
        for w in range(3):
            denom = len(tuples)**2*math.comb(2,w)*3**w
            for cutoff in range(3):
                self.assertEqual(F(sum(counts[w][:cutoff+1]),denom),
                                 exact_probability(2,w,cutoff,4))

    def test_positive_vector_bound(self):
        for q in (4,16,65536):
            for lx,ly in product((-8.,-2.,0.,2.,8.),(-8.,-2.,0.)):
                t0,t1 = local_transfer(math.exp(ly),q)
                matrix = t0+math.exp(lx)*t1
                for n in (1,3,20):
                    actual = float(np.linalg.matrix_power(matrix,n)[0].sum())
                    self.assertGreaterEqual(final_log_moment(lx,ly,n,q)+1e-10,
                                            math.log(actual))

    def test_all_probability_masses(self):
        for q in (4,8,16,256,65536):
            for old,active in product(range(4),range(2)):
                law = support_law(old,active,q)
                self.assertEqual(sum(law.values()),1)
                self.assertTrue(all(p>0 for p in law.values()))

    def test_inverse_tail_really_can_have_rank_two(self):
        setup = ((3,3,3,3),(3,2,1,2),(1,3,1,1),(1,2,1,2))
        columns = []
        for j in range(4):
            s1 = s2 = 0
            inputs = []
            for i,(c,d1,d2,f) in enumerate(setup):
                x = mul(mul(c,c),int(i == j)^s1^mul(f,s2))
                inputs.append(x)
                s1 ^= mul(d1,x)
                s2 ^= mul(d2,x)
            self.assertEqual(encode(inputs,setup),tuple(int(i == j) for i in range(4)))
            columns.append(inputs)
        determinant = mul(columns[0][2],columns[1][3]) ^ mul(columns[0][3],columns[1][2])
        self.assertEqual(determinant,3)


if __name__ == '__main__':
    unittest.main()

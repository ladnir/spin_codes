"""Independent finite checks of the rational scalar support law."""
from fractions import Fraction as F
from itertools import product
import math
import unittest

import numpy as np

from hd_rational_post import (exact_probability, final_log_moment,
                               local_transfer, support_law)


def mul(a, b):
    # GF4 = F2[X]/(X^2+X+1).
    value = 0
    for _ in range(2):
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return value


def encode(x, coefficients):
    state = 0
    output = []
    for a, (c,d) in zip(x, coefficients):
        output.append(state ^ mul(c,a))
        state ^= mul(d,a)
    return tuple(output)


class RationalPostTests(unittest.TestCase):
    def test_all_gf4_local_values(self):
        pairs = [(c,d) for c in range(1,4) for d in range(1,4) if c != d]
        for state, x in product(range(4), repeat=2):
            counts = {}
            for c,d in pairs:
                key = int(bool(state ^ mul(d,x))), int(bool(state ^ mul(c,x)))
                counts[key] = counts.get(key, F(0))+F(1,len(pairs))
            self.assertEqual(counts, support_law(int(bool(state)),int(bool(x)),4))

    def test_all_length_three_setups_and_messages(self):
        pairs = [(c,d) for c in range(1,4) for d in range(1,4) if c != d]
        n = 3
        counts = [[0]*(n+1) for _ in range(n+1)]
        setups = list(product(pairs, repeat=n))
        for coefficients in setups:
            outputs = set()
            for x in product(range(4), repeat=n):
                y = encode(x, coefficients)
                outputs.add(y)
                counts[sum(v != 0 for v in x)][sum(v != 0 for v in y)] += 1
            self.assertEqual(len(outputs), 4**n)
        for w in range(n+1):
            denom = len(setups)*math.comb(n,w)*3**w
            for cutoff in range(n+1):
                self.assertEqual(F(sum(counts[w][:cutoff+1]),denom),
                                 exact_probability(n,w,cutoff,4))

    def test_positive_vector_moment(self):
        for q in (4,16,65536):
            for lx,ly in product((-9.,-2.,0.,2.,9.),(-9.,-2.,0.)):
                t0,t1 = local_transfer(math.exp(ly),q)
                matrix = t0+math.exp(lx)*t1
                for n in (1,3,20):
                    actual = float(np.linalg.matrix_power(matrix,n)[0].sum())
                    bound = final_log_moment(lx,ly,n,q)
                    self.assertGreaterEqual(bound+1e-10, math.log(actual))

    def test_dense_inverse_input_is_not_limited_to_twice_output(self):
        # One output at the start can have an arbitrarily long input tail.
        pairs = [(1,2)]*5
        target = (1,0,0,0,0)
        hits = [x for x in product(range(4),repeat=5) if encode(x,pairs) == target]
        self.assertEqual(len(hits),1)
        self.assertEqual(sum(v != 0 for v in hits[0]),5)


if __name__ == '__main__':
    unittest.main()

"""Exact normalization and transfer checks for sparse return enumeration."""

from collections import Counter
import itertools
import math
import unittest
import numpy as np

from gf16_return_spectrum import histogram, local, return_matrix
from gf16_sizes import local_model, singleton_weights


class ReturnSpectrumTests(unittest.TestCase):
    def test_singletons_match_existing_formula(self):
        for state in (2, 4, 6, 8):
            expected = Counter(singleton_weights(state))
            self.assertEqual(histogram(state, 1), tuple(expected[w] for w in range(17)))

    def test_unquotiented_pair_enumeration(self):
        matrix, table = return_matrix(4)
        pairs = np.array(list(itertools.product(range(1, 16), repeat=2)))
        counts = np.zeros(17, dtype=np.int64)
        for left, right in itertools.combinations(range(16), 2):
            values = table[pairs[:, 0, None], matrix[:, left]] ^ table[pairs[:, 1, None], matrix[:, right]]
            counts += np.bincount(np.count_nonzero(values, axis=1), minlength=17)
        self.assertEqual(tuple(counts), tuple(15*x for x in histogram(4, 2)))

    def test_counts_and_complete_output_moment(self):
        for state in (4, 6):
            for weight in (2, 3, 4):
                counts = histogram(state, weight)
                self.assertEqual(sum(counts), math.comb(16, weight)*15**(weight-1))
                self.assertEqual(counts[0], 0)
            for z in (.4, .9999, 1.):
                exact, upper = local(z, state), local_model(z, state)
                expected = upper.sum(axis=2)
                for weight in (2, 3, 4):
                    # The old entrywise envelope has a larger row sum:
                    # G-J_lower+J_upper. The exact replacement sums to G.
                    expected[weight, 1] = upper[weight, 1, 1]+z**16/(16**state-1)
                np.testing.assert_allclose(exact.sum(axis=2), expected, rtol=1e-14, atol=1e-20)
                self.assertTrue(np.all(exact >= 0))
                self.assertTrue(np.all(exact[2:5, 1, 0] <= upper[2:5, 1, 0]*(1+1e-13)))


if __name__ == '__main__':
    unittest.main()

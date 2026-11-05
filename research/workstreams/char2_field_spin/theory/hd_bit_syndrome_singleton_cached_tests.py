"""Small exact-reference controls; no regional provider."""
from itertools import product
from math import comb, log2
import unittest
import numpy as np
from hd_bit_syndrome_singleton_cached import column_coefficients, ordered_word_coefficients, marker_minimum


class SingletonTests(unittest.TestCase):
    def test_column_expansion_against_labelled_activity_and_atom_words(self):
        # Each directed transition has one atom; masses form a noncommuting law.
        mass = ((2, 3), (5, 7))
        atom = ((0, 1), (2, 3))
        groups, regions, q = 2, 3, 4
        def moment(k):
            with np.errstate(divide='ignore'):
                return np.log2(np.array([[mass[i][j]*(atom[i][j]**k) for j in range(2)]
                                        for i in range(2)], dtype=float))
        obtained = np.exp2(column_coefficients(moment, groups, regions, q))
        exact = np.zeros_like(obtained)
        for start in range(2):
            for path in product(range(2), repeat=regions):
                states = (start,)+path
                weight, t = 1, 1
                for i, j in zip(states, states[1:]):
                    weight *= mass[i][j]
                    t *= atom[i][j]
                for active in product((0, 1), repeat=groups):
                    k = sum(active)
                    exact[start, path[-1], k] += weight*(q-1)**k*t**k
        np.testing.assert_allclose(obtained, exact, rtol=1e-14, atol=1e-8)

    def test_ordered_word_coefficients_exact_noncommuting_enumeration(self):
        matrices = (np.array([[1, 2], [3, 4]], dtype=object),
                    np.array([[2, 1], [0, 3]], dtype=object))
        with np.errstate(divide='ignore'):
            column = np.log2(np.stack(matrices, axis=2).astype(float))
        actual = np.exp2(ordered_word_coefficients(column, 3))
        expected = [0]*4
        for degrees in product((0, 1), repeat=3):
            row = np.array([1, 0], dtype=object)
            for k in degrees:
                row = row@matrices[k]
            expected[sum(degrees)] += sum(row)
        np.testing.assert_allclose(actual, expected, rtol=1e-14)

    def test_marker_minimum_is_valid_upper_for_tail(self):
        values = np.array([log2(comb(6, k)) for k in range(7)])
        fit = marker_minimum(values, 4)
        exact = log2(sum(comb(6, k) for k in range(4, 7)))
        self.assertGreaterEqual(fit['marked_tail_log2']+1e-12, exact)
        self.assertAlmostEqual(fit['log2_marker'], 1., places=10)
        self.assertEqual(marker_minimum(values, 1)['log2_marker'], 0.)

    def test_guards(self):
        with self.assertRaises(ValueError):
            column_coefficients(lambda _: np.zeros((2, 2)), 0)
        with self.assertRaises(ValueError):
            ordered_word_coefficients(np.zeros((2, 2, 2)), 0)
        with self.assertRaises(ValueError):
            marker_minimum([0., 0.], 2)


if __name__ == '__main__':
    unittest.main()

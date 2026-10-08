"""Small exact law checks and floating-screen regression checks."""

from collections import Counter
from fractions import Fraction
import itertools
import math
import unittest

import numpy as np

import field_rate_tradeoff as screen


def gf4_mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def parity(diagonals, word, parity_length):
    out = []
    for i in range(parity_length):
        value = 0
        for j, coefficient in enumerate(word):
            value ^= gf4_mul(diagonals[i-j+len(word)-1], coefficient)
        out.append(value)
    return tuple(out)


def matrix_product(a, b):
    return [[sum(a[i][j]*b[j][k] for j in range(2)) for k in range(2)] for i in range(2)]


class RectangularLawTests(unittest.TestCase):
    def test_every_nonzero_rectangular_word_has_uniform_parity(self):
        # k>p explicitly includes last nonzero source coordinates outside
        # the parity index range. k<p checks the other rectangular shape.
        for dimension, parity_length in ((3, 1), (3, 2), (2, 3)):
            all_diagonals = tuple(itertools.product(range(4), repeat=dimension+parity_length-1))
            for word in itertools.product(range(4), repeat=dimension):
                if not any(word):
                    continue
                counts = Counter(parity(d, word, parity_length) for d in all_diagonals)
                self.assertEqual(len(counts), 4**parity_length)
                self.assertEqual(set(counts.values()), {4**(dimension-1)})

    def test_selected_diagonal_map_is_triangular(self):
        for dimension, parity_length in ((3, 1), (3, 2), (2, 3)):
            for word in itertools.product(range(4), repeat=dimension):
                if not any(word):
                    continue
                last = max(j for j, v in enumerate(word) if v)
                matrix = [[word[i-r+last] if 0 <= i-r+last < dimension else 0
                           for r in range(parity_length)] for i in range(parity_length)]
                for i in range(parity_length):
                    self.assertEqual(matrix[i][i], word[last])
                    self.assertTrue(all(matrix[i][r] == 0 for r in range(i)))

    def test_complete_value_counting_identity(self):
        dimension, parity_length = 3, 2
        diagonals = tuple(itertools.product(range(4), repeat=dimension+parity_length-1))
        words = [w for w in itertools.product(range(4), repeat=dimension) if any(w)]
        def marked(word, value):
            # A non-product function exercises the complete value measure.
            return (1+sum((j+1)*v for j, v in enumerate(word+value)))**2
        left = Fraction(sum(marked(w, parity(d, w, parity_length)) for d in diagonals for w in words), len(diagonals))
        right = Fraction(sum(marked(w, v) for w in words for v in itertools.product(range(4), repeat=parity_length)), 4**parity_length)
        self.assertEqual(left, right)


class SourceDPTests(unittest.TestCase):
    def test_seen_count_dp_matches_exhaustive_supports(self):
        active, length, dimension = 2, 4, 3
        regions = [
            [[Fraction(1, 2), Fraction(1, 7)], [Fraction(1, 9), Fraction(1, 3)]],
            [[Fraction(1, 3), Fraction(1, 8)], [Fraction(1, 10), Fraction(1, 4)]],
            [[Fraction(1, 4), Fraction(1, 9)], [Fraction(1, 11), Fraction(1, 5)]],
        ]
        exact = Fraction(0)
        for masks in itertools.product(range(1 << active), repeat=length):
            seen = 0
            for mask in masks[:dimension]:
                seen |= mask
            if seen != (1 << active)-1:
                continue
            matrix = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]]
            for mask in masks:
                h = mask.bit_count()
                matrix = matrix_product(matrix, [[15**h*v for v in row] for row in regions[h]])
            exact += sum(matrix[0])
        exact /= 16**((length-dimension)*active)
        logs = np.log(np.array(regions, dtype=float))
        actual = screen.source_moment(logs, active, length, dimension)
        self.assertAlmostEqual(actual, math.log(float(exact)), places=12)
        self.assertLessEqual(actual, screen.source_moment(logs, active, length, dimension, False))


class FormulaTests(unittest.TestCase):
    def test_geometry_padding(self):
        expected = {(32, 24): (1048704, 1398272, 43696, 139827),
                    (64, 48): (1049088, 1398784, 21856, 139878),
                    (128, 96): (1049088, 1398784, 10928, 139878)}
        for shape, result in expected.items():
            self.assertEqual(screen.geometry(20, *shape), result)
        for k_log in (16, 18, 20):
            for length, dimension in expected:
                k, n, rows, threshold = screen.geometry(k_log, length, dimension)
                self.assertEqual(rows % 16, 0)
                self.assertEqual(k, rows*dimension)
                self.assertEqual(n, rows*length)
                self.assertLess(k-(1 << k_log), dimension*16)
                self.assertEqual(threshold, (n+9)//10-1)

    def test_preserved_half_rate_formula(self):
        from rank_continuation_random_outer import screen as prior_sparse, middle_screen as prior_middle
        old = prior_sparse(20, 2, 3, 32, .9998)
        new = screen.sparse_screen(20, 2, 3, 64, 32, .9998, False)
        self.assertAlmostEqual(old, new, places=7)
        old = prior_middle(20, 1000, 3, 32, .01, .9)
        new = screen.middle_screen(20, 1000, 3, 64, 32, .01, .9)
        self.assertAlmostEqual(old, new, places=7)

    def test_binomial_log(self):
        for n, k in ((43696, 1), (43696, 16), (43696, 43695), (43696, 2000)):
            self.assertAlmostEqual(screen.log_choose(n, k), math.log(math.comb(n, k)), places=8)

    def test_mds_source_count_and_tilt(self):
        from q1_batched_screen import fixed_support_multiplicity
        counts = [fixed_support_multiplicity(16, 12, w, 16) for w in range(1, 17)]
        self.assertEqual(counts[:4], [0]*4)
        self.assertEqual(sum(math.comb(16, w)*counts[w-1] for w in range(1, 17)), 16**12-1)
        for rho in (.5, .9, 1., 1.2):
            eta = math.exp(screen.mds_log_eta(rho))
            for w in range(5, 17):
                self.assertLessEqual(counts[w-1], eta*16**(w-4)*rho**w*(1+1e-12))

    def test_uniform_input_row_sums(self):
        for state in (3, 4, 5, 6):
            for z in (.05, .2, .9, 1.):
                result = screen.uniform_input_check(z, state)
                exact = np.array(result['exact'])
                np.testing.assert_allclose(exact.sum(axis=1), [result['exact_row_sum']]*2, rtol=1e-14)
                np.testing.assert_allclose(result['envelope'][0], exact[0], rtol=1e-14)
                self.assertGreaterEqual(result['envelope_nonzero_row_ratio'], 1-1e-14)
                if z == 1.:
                    self.assertAlmostEqual(exact[0, 0], 16**-state, places=14)

    def test_aggregate_hook_normalization(self):
        from gf16_sizes import local_model
        def aggregate(v, z, state):
            local = local_model(z, state)
            return sum(math.comb(16, w)*(15*v)**w*local[w] for w in range(17))
        args = (20, 1000, 3, .01, .9, .87)
        self.assertAlmostEqual(screen.mds_middle_screen(*args),
                               screen.mds_middle_screen(*args, aggregate=aggregate), places=7)
        args = (20, 1000, 3, 32, 24, .01, .9)
        self.assertAlmostEqual(screen.middle_screen(*args),
                               screen.middle_screen(*args, aggregate=aggregate), places=7)


if __name__ == '__main__':
    unittest.main()

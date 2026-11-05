from fractions import Fraction as F
from itertools import product
from unittest import TestCase, main

import numpy as np

import hd_post_refresh_compress as compression


def multiply(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def exact_compress(matrix, weights):
    return [[matrix[0][0], max(matrix[0][j+1]/weights[j] for j in range(7))],
        [sum(weights[i]*matrix[i+1][0] for i in range(7)),
         max(sum(weights[i]*matrix[i+1][j+1] for i in range(7))/weights[j]
             for j in range(7))]]


def field_multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


class PostRefreshCompressionTests(TestCase):
    def setUp(self):
        self.weights = [F(1 << i, 127) for i in range(7)]
        self.matrices = [
            [[F((3*i+5*j+1)%11, 13) for j in range(8)] for i in range(8)],
            [[F((7*i+2*j+3)%13, 17) for j in range(8)] for i in range(8)]]

    def test_exact_embedding_and_noncommuting_chronological_words(self):
        embedding = [[F(1)]+[0]*7, [0]+self.weights]
        reduced = [exact_compress(matrix, self.weights) for matrix in self.matrices]
        self.assertNotEqual(multiply(*self.matrices), multiply(*reversed(self.matrices)))
        for matrix, small in zip(self.matrices, reduced):
            left, right = multiply(embedding, matrix), multiply(small, embedding)
            self.assertTrue(all(a <= b for x, y in zip(left, right) for a, b in zip(x, y)))
        for word in product(range(2), repeat=4):
            full, small = [[F(1)]+[0]*7], [[F(1), 0]]
            for index in word:
                full = multiply(full, self.matrices[index])
                small = multiply(small, reduced[index])
            self.assertLessEqual(sum(full[0]), sum(small[0]))

    def test_float_matches_exact_and_identity_zeros(self):
        weights = np.log2(list(map(float, self.weights)))
        for matrix in self.matrices:
            with np.errstate(divide='ignore'):
                actual = compression.compress_matrix(np.log2(np.array(matrix, dtype=float)), weights)
            np.testing.assert_allclose(np.exp2(actual), np.array(exact_compress(matrix, self.weights),
                dtype=float), rtol=2e-14)
        with np.errstate(divide='ignore'):
            identity = compression.compress_matrix(np.log2(np.eye(8)), weights)
        np.testing.assert_allclose(np.exp2(identity), np.eye(2), atol=1e-15)

    def test_fixed_predata_rule_normalized_once_and_scale_invariant(self):
        expected = (0., -128., 0., 0., -128., -128., -256.)
        self.assertEqual(compression.FIXED_LOG_WEIGHTS, expected)
        one, total = compression.normalize_live_weights(expected)
        two, total2 = compression.normalize_live_weights(tuple(x+20 for x in expected))
        np.testing.assert_allclose(one, two, atol=3e-14)
        self.assertAlmostEqual(total2-total, 20.)
        self.assertAlmostEqual(np.exp2(one).sum(), 1.)
        bank = {(m, f): np.zeros((8, 8)) for m in range(3) for f in range(min(1, m//2)+1)}
        result = compression.compress_bank(bank, maximum=2)
        self.assertEqual(result['metadata']['terminal_vector'], [1, 1])
        self.assertEqual(result['metadata']['extra_terminal_factor_log2'], 0.)
        self.assertEqual(compression.entrywise_curves(result['regions'], maximum=2).shape, (4, 3))

    def test_point_and_conditional_flat_families_close_under_positive_mixtures(self):
        masses = [F(1, 3), F(2, 5)]
        total = sum(masses)
        # P/G retain arbitrary atoms; F retains independent flat second
        # conditional on each first-coordinate atom, rather than flat first.
        distribution = {(u, v): mass/3 for u, mass in enumerate(masses, 1) for v in range(1, 4)}
        first = {u: sum(distribution[u, v] for v in range(1, 4)) for u in (1, 2)}
        self.assertEqual(sum(first.values()), total)
        self.assertNotEqual(first[1], first[2])
        for (u, v), mass in distribution.items():
            self.assertEqual(mass, total*(first[u]/total)/3)

    def test_exact_two_packet_identity_for_all_q4_old_states(self):
        for u, v, c1, c2 in product(range(4), repeat=4):
            for m1, m2 in product(range(1, 4), repeat=2):
                d1, d2 = field_multiply(m1, u ^ c1), field_multiply(m2, v ^ c2)
                u1, v1 = v, d1 ^ v
                u2, v2 = v1, field_multiply(m2, u1 ^ c2) ^ v1
                self.assertEqual((u2, v2), (v ^ d1, v ^ d1 ^ d2))
        # An arbitrary old state can return with ONE refresh hit. This is
        # not a counterexample to the first excursion from zero lemma.
        returns = sum((1 ^ field_multiply(m1, 1)) == 0
                      for m1, m2 in product(range(1, 4), repeat=2))
        self.assertEqual(F(returns, 9), F(1, 3))

    def test_fixed_word_first_return_two_hit_bound_tiny(self):
        epsilon = F(1, 3)
        for feedback in ((1, 1, 1), (1, 2, 3, 0), (1, 0, 2, 1),
                         (2, 3, 1, 2, 1), (1, 1, 1, 1, 1)):
            count = [0]*(len(feedback)+1)
            for keys in product(range(1, 4), repeat=len(feedback)):
                u = v = 0
                for n, (c, key) in enumerate(zip(feedback, keys), 1):
                    u, v = v, field_multiply(key, u ^ c) ^ v
                    if not u and not v:
                        count[n] += 1
                        break
            self.assertEqual(count[1:3], [0, 0])
            for n in range(3, len(feedback)+1):
                self.assertLessEqual(F(count[n], 3**len(feedback)), (n-2)*epsilon**2)

    def test_fail_closed_guards(self):
        for values in ((0.,)*6, (0.,)*6+(-np.inf,)):
            with self.assertRaises(ValueError):
                compression.normalize_live_weights(values)
        with self.assertRaises(ValueError):
            compression.compress_matrix(np.zeros((8, 8)), np.zeros(7))
        with self.assertRaises(ValueError):
            compression.compress_bank({}, maximum=2)
        with self.assertRaises(ValueError):
            compression.entrywise_curves({(0, 0): np.zeros((2, 2))}, maximum=2)


if __name__ == '__main__':
    main()

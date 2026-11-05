"""Independent tiny polynomial tests; construction-size work is separate."""
from itertools import product
from math import comb, log2
import unittest
import numpy as np
import hd_field128_full_quartet_coefficients_v55 as subject


class Checks(unittest.TestCase):
    def toy(self, slots):
        return np.array([[[1.+j/8, .25+j/9], [.125+j/11, .75+j/7]]
                         for j in range(2*slots+1)])

    def brute(self, matrices, slots, alpha):
        geometry = [(s, f, comb(slots, s)*comb(slots-s, f))
                    for s in range(slots+1) for f in range(slots-s+1)]
        result = np.zeros((4, 4*slots+1, 4*slots+1))
        for rows in product(geometry, repeat=4):
            raw = np.eye(2)
            weight, s, f = 1, 0, 0
            for single, full, count in rows:
                raw = raw@matrices[single+2*full]
                weight *= count
                s += single
                f += full
            result[:, s, f] += weight*(raw**alpha).reshape(4)
        return result

    def test_independent_histogram_enumeration(self):
        for slots in (1, 2):
            matrices = self.toy(slots)
            for alpha in (.125, .5, 1.):
                logs, meta = subject.all_quartet_coefficients(np.log2(matrices), alpha,
                                                             slots=slots, chunk=7)
                np.testing.assert_allclose(np.exp2(logs), self.brute(matrices, slots, alpha),
                                           rtol=3e-14, atol=3e-11)
                self.assertEqual(meta['ordered_quartets'], (2*slots+1)**4)
                self.assertEqual(meta['multiset_groups'], comb(2*slots+4, 4))
                self.assertEqual(meta['finite_raw_underflows'], 0)

    def test_repeated_multisets_and_degree_boundary(self):
        matrices = self.toy(2)
        logs, _ = subject.all_quartet_coefficients(np.log2(matrices), .25, slots=2)
        coefficients = np.exp2(logs)
        for j, position in ((0, (0, 0)), (4, (0, 8))):
            expected = np.linalg.matrix_power(matrices[j], 4)**.25
            np.testing.assert_allclose(coefficients[:, position[0], position[1]], expected.ravel(),
                                       rtol=1e-14)
        self.assertEqual(coefficients[:, 8, 8].sum(), 0.)

    def test_chunk_independence_and_point(self):
        matrices, alpha, odds = self.toy(2), .375, (log2(3), log2(5))
        first, _ = subject.all_quartet_coefficients(np.log2(matrices), alpha, slots=2, chunk=1)
        second, _ = subject.all_quartet_coefficients(np.log2(matrices), alpha, slots=2, chunk=32)
        np.testing.assert_allclose(first, second, rtol=0, atol=1e-13)
        weights = [0.]*5
        for assignment in product(range(3), repeat=2):
            s, f = assignment.count(1), assignment.count(2)
            weights[s+2*f] += 3**s*5**f
        expected = np.zeros((2, 2))
        for indices in product(range(5), repeat=4):
            raw, weight = np.eye(2), 1.
            for j in indices:
                raw = raw@matrices[j]
                weight *= weights[j]
            expected += weight*raw**alpha
        np.testing.assert_allclose(np.exp2(subject.point_logs(first, odds)), expected,
                                   rtol=2e-14)

    def test_alpha_one_point_is_fourth_matrix_power(self):
        matrices, slots, odds = self.toy(2), 2, (log2(2), log2(3))
        polynomial, _ = subject.all_quartet_coefficients(np.log2(matrices), 1., slots=slots)
        single = np.zeros((2, 2))
        for s in range(slots+1):
            for f in range(slots-s+1):
                single += comb(slots, s)*comb(slots-s, f)*2**s*3**f*matrices[s+2*f]
        np.testing.assert_allclose(np.exp2(subject.point_logs(polynomial, odds)),
                                   np.linalg.matrix_power(single, 4), rtol=2e-14)

    def test_algebraic_zeros_and_range_rejection(self):
        logs = np.full((3, 2, 2), -np.inf)
        logs[:, 0, 0] = 0.
        logs[:, 1, 1] = -1.
        coefficients, _ = subject.all_quartet_coefficients(logs, .5, slots=1)
        self.assertTrue(np.all(np.isneginf(coefficients[[1, 2]])))
        with self.assertRaises(ValueError):
            subject.all_quartet_coefficients(logs-1000, .5, slots=1)
        with self.assertRaises(ValueError):
            subject.all_quartet_coefficients(logs, .5, slots=1, chunk=0)

    def test_alias_degree128_folds_and_uses_power32(self):
        first = np.array([[.5, .25], [.125, .75]])
        second = np.array([[.25, .125], [.25, .5]])
        polynomial = np.full((4, 129, 129), -np.inf)
        polynomial[:, 0, 0] = np.log2(first).ravel()
        polynomial[:, 128, 0] = np.log2(second).ravel()
        odds = (.001, 0.)
        table = subject.quartet_alias_table(polynomial, odds)
        expected = np.linalg.matrix_power(first+2**(.001*128)*second, 32)
        np.testing.assert_allclose(np.exp2(table.positive_numerator).reshape(2, 2),
                                   expected, rtol=1e-13)
        np.testing.assert_allclose(np.exp2(table.logs[:, 0, 0]).reshape(2, 2),
                                   expected, rtol=2e-11)
        self.assertEqual(int(table.reliable.sum()), 4)
        self.assertTrue(np.all(table.reliable[:, 0, 0]))
        self.assertEqual(table.diagnostic['quartets_per_region'], 32)


if __name__ == '__main__':
    unittest.main()

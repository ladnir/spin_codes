"""Order, regression, and conditioning checks for fine-count path grouping."""
from itertools import product
from math import comb, exp
import unittest

import numpy as np

import fine_grouped_gate as fine


def matrix_product(matrices):
    result = np.eye(matrices[0].shape[0])
    for m in matrices:
        result = result@m
    return result


class FineGroupedTests(unittest.TestCase):
    def setUp(self):
        self.local = np.array([[[.8, .0], [.0, .3]],
            [[.15, .25], [.03, .08]], [[.08, .02], [.09, .04]]])

    def test_explicit_fine_counts_and_marker(self):
        for g in (1, 2, 4):
            tuples = fine.tuple_products(self.local, g, cap=1)
            actual = fine.fine_operators(tuples, .4, nu=.7)
            expected = np.zeros_like(actual)
            for counts in product(range(3), repeat=g):
                j = sum(counts)
                choices = np.prod([comb(2, k) for k in counts])
                mark = sum(min(k, 1) for k in counts)
                expected[j] += (choices/comb(2*g, j)*exp(.7*mark)
                    *matrix_product([self.local[k] for k in counts])**.4)
            np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=1e-15)

    def test_fine_group_no_worse_than_physical_or_coarse(self):
        for g in (2, 4):
            tuples = fine.tuple_products(self.local, g)
            coarse, _ = fine.coarse.grouped_operators(self.local, g)
            for alpha in (.3, .4, .5, 1.):
                actual = fine.fine_operators(tuples, alpha)
                physical, _ = fine.coarse.grouped_operators(self.local**alpha, g)
                self.assertTrue(np.all(actual <= physical+2e-14))
                self.assertTrue(np.all(actual <= coarse**alpha+2e-14))

    def test_alpha1_regresses_marked_physical_product(self):
        for g in (1, 2, 4):
            tuples = fine.tuple_products(self.local, g, cap=1)
            for nu in (0., .7, 2.):
                actual = fine.fine_operators(tuples, 1., nu=nu)
                marked = self.local*np.exp(nu*np.minimum(np.arange(3), 1))[:, None, None]
                expected, _ = fine.coarse.grouped_operators(marked, g)
                np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=1e-15)

    def test_pointwise_clipping_with_h2_marker(self):
        alpha, nu, beta, h = .4, .7, 20., 5
        for counts in product(range(3), repeat=4):
            matrices = [self.local[k] for k in counts]
            cost = sum(min(k, 2) for k in counts)
            moment = matrix_product(matrices)[0].sum()
            clipped = min(1., beta*moment) if cost >= h else 0.
            path = ((matrices[0]@matrices[1])**alpha
                @(matrices[2]@matrices[3])**alpha)[0].sum()
            bound = beta**alpha*exp(nu*(cost-h))*path
            self.assertLessEqual(clipped, bound+2e-13)

    def test_zero_entries_remain_zero(self):
        local = np.array([np.diag([.3, .4]), np.diag([.1, .2])])
        actual = fine.fine_operators(fine.tuple_products(local, 2, cap=1), .4)
        self.assertTrue(np.all(actual[:, 0, 1] == 0))
        self.assertTrue(np.all(actual[:, 1, 0] == 0))


if __name__ == '__main__':
    unittest.main()

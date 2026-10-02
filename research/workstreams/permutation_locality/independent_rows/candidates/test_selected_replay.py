import unittest
from fractions import Fraction as Q
from itertools import combinations
from math import comb

from flint import arb, arb_mat, fmpq_mat, ctx

from selected_replay import region_action, point
from test_fresh_history import endpoint


class SelectedReplayTests(unittest.TestCase):
    def test_noncommuting_placement_matches_exact_enumeration(self):
        ctx.prec = 192
        rows = [[[2, 1], [1, 3]], [[1, 3], [2, 1]], [[1, 0], [2, 2]]]
        operators = [arb_mat(t) for t in rows]
        exact = [fmpq_mat(t) for t in rows]
        potential = [point(Q(1, 3)), point(Q(2, 7))]
        vector = fmpq_mat([[1], [1]])
        vector[0, 0] /= 3
        vector[1, 0] *= Q(2, 7).numerator
        vector[1, 0] /= Q(2, 7).denominator
        for epochs in (1, 2, 3):
            result = region_action(operators, potential, 2*epochs, epochs=epochs, windows=2)
            for r, values in enumerate(result):
                total = fmpq_mat(2, 1)
                for chosen in combinations(range(2*epochs), r):
                    matrix = fmpq_mat([[1, 0], [0, 1]])
                    for epoch in range(epochs):
                        matrix *= exact[sum(x//2 == epoch for x in chosen)]
                    total += matrix*vector
                total /= comb(2*epochs, r)
                for i in range(2):
                    expected = Q(str(total[i, 0]))
                    self.assertGreaterEqual(endpoint(values[i]), expected)
                    self.assertLess(endpoint(values[i])-expected, Q(1, 2**170))

    def test_invalid_geometry_or_potential(self):
        operators = [arb_mat([[1]])]*3
        for degree, epochs, windows in ((7, 3, 2), (-1, 3, 2), (1, 0, 2)):
            with self.assertRaises(ValueError):
                region_action(operators, [arb(1)], degree, epochs=epochs, windows=windows)
        with self.assertRaises(ValueError):
            region_action(operators, [arb(0)], 1)


if __name__ == '__main__':
    unittest.main()

"""Independent exact finite tests for return-pivot cover evaluation."""

from fractions import Fraction as F
import math
import unittest

import return_cover
import return_cover_certificate as cert


class ReturnCoverTests(unittest.TestCase):
    def test_alive_return_potential(self):
        self.assertEqual(return_cover.exact_return_checks(), 720)

    def test_exact_weight_intersection_enumerator(self):
        witness = cert.Witness(F(7, 3), F(1, 4), F(5, 2), F(2, 5))
        x, y, tau = witness.x, witness.y, witness.tau
        expected = [[F(1), F(0)], [F(0), F(0)]]
        for wi in range(17):
            for wo in range(17):
                for both in range(max(0, wi+wo-16), min(wi, wo)+1):
                    count = math.comb(16, wi)*math.comb(wi, both)*math.comb(16-wi, wo-both)
                    if wi == wo == both and wi:
                        value = count*(x*y)**wi
                        expected[0][1] += value
                        if wi >= 5:
                            expected[0][0] += value/tau**4
                    if wi+wo-both >= 13:
                        value = count*x**wi*y**wo/tau**(wi-both)
                        expected[1][1] += value
                        if wi and wo and (both or wo != 1):
                            expected[1][0] += value/tau**min(both, 4)
        actual = cert.transfer(witness)
        for i in range(2):
            for j in range(2):
                self.assertLessEqual(F(actual[i][j].lo), expected[i][j])
                self.assertGreaterEqual(F(actual[i][j].hi), expected[i][j])

    def test_interval_terminal_power(self):
        matrix = [[F(3, 2), F(2, 7)], [F(5, 4), F(1, 3)]]
        current = [F(1), F(0)]
        for _ in range(7):
            current = [sum(current[i]*matrix[i][j] for i in range(2)) for j in range(2)]
        actual = cert.terminal_power([[cert.interval(value) for value in row] for row in matrix], 7)
        self.assertLessEqual(F(actual.lo), sum(current))
        self.assertGreaterEqual(F(actual.hi), sum(current))

    def test_density_endpoint_bound(self):
        for witness in (cert.WITNESSES[0][1], cert.Witness(F(2), F(1, 2), F(2), F(1, 10))):
            low, high = F(19, 100), F(191, 1000)
            upper = cert.row_upper(low, high, witness)
            for a in (low, (low+high)/2, high):
                exact = (1-a)**32+witness.tau**(-10)*(1+(witness.tau*witness.tilt-1)*a)**32
                self.assertGreaterEqual(F(upper.hi), exact)

    def test_integer_cells_are_contiguous(self):
        n = 96256
        cells = [((i*n+9999)//10000, ((i+1)*n+9999)//10000-1) for i in range(1900, 2400)]
        self.assertEqual(cells[0][0], 18289)
        self.assertEqual(cells[-1][1], 23101)
        for left, right in zip(cells, cells[1:]):
            self.assertEqual(left[1]+1, right[0])


if __name__ == '__main__':
    unittest.main()

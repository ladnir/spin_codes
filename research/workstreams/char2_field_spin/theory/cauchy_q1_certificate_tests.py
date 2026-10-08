"""Exact and outward-rounding checks for the Q1-only certificate."""

from fractions import Fraction as F
from itertools import product
import math
import unittest

import cauchy_q1_certificate as certificate


def multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


def power(matrix, exponent):
    result = [[F(1), F(0)], [F(0), F(1)]]
    for _ in range(exponent):
        result = multiply(result, matrix)
    return result


def exact_marginals(y):
    z = (1-y)/255
    m0 = sum(math.comb(16, j)*y**(16-j)*z**j for j in range(4))
    m1 = (y+z)*sum(math.comb(15, j)*y**(15-j)*z**j for j in range(4))
    return m0, m1


class Q1CertificateTests(unittest.TestCase):
    def assert_encloses(self, interval, exact):
        self.assertLessEqual(F(interval.lo), exact)
        self.assertGreaterEqual(F(interval.hi), exact)

    def test_positive_geometric_sum(self):
        for value in (F(1), F(2, 3), F(9997, 10000)):
            for count in (0, 1, 2, 3, 7, 16, 31):
                powered, summed = certificate.positive_power_sum(certificate.interval(value), count)
                self.assert_encloses(powered, value**count)
                self.assert_encloses(summed, sum((value**j for j in range(count)), F(0)))

    def test_marginal_enclosures(self):
        for mark in (F(1), F(2, 3), F(9997, 10000)):
            m0, m1 = exact_marginals(mark)
            self.assert_encloses(certificate.marginal(mark, 0), m0)
            self.assert_encloses(certificate.marginal(mark, 1), m1)

    def test_support_count_inclusion_exclusion(self):
        for parity in (2, 10):
            for weight in range(parity+1, 33):
                all_nonzero = sum((-1)**j*math.comb(weight, j)
                                  *256**max(0, weight-parity-j)
                                  for j in range(weight+1))
                self.assertEqual(all_nonzero % 255, 0)
                self.assertEqual(certificate.projective_support_count(weight, parity),
                                 all_nonzero//255)

    def test_complete_small_fraction_expression(self):
        mark, length, parity, threshold = F(2, 3), 5, 2, 9
        m0, m1 = exact_marginals(mark)
        t0 = [[F(1), F(0)], [F(0), m0]]
        t1 = [[F(0), mark], [mark**16/255**4, m1]]
        for packets in (1, 2, 3):
            rows = 16*packets
            empty = power(t0, packets)
            active = [[F(0), F(0)], [F(0), F(0)]]
            for position in range(packets):
                matrix = multiply(multiply(power(t0, position), t1),
                                  power(t0, packets-position-1))
                for i in range(2):
                    for j in range(2):
                        active[i][j] += matrix[i][j]/packets
            exact = F(0)
            for mask in product((0, 1), repeat=length):
                weight = sum(mask)
                if weight <= parity:
                    continue
                matrix = [[F(1), F(0)], [F(0), F(1)]]
                for bit in mask:
                    matrix = multiply(matrix, active if bit else empty)
                count = sum((-1)**j*math.comb(weight, j)*256**max(0, weight-parity-j)
                            for j in range(weight+1))//255
                exact += count*sum(matrix[0])
            exact *= rows/mark**threshold
            actual = certificate.q1_expression(rows, mark, length=length,
                                               parity=parity, threshold=threshold)
            self.assert_encloses(actual, exact)

    def test_actual_pairing_and_return_weight(self):
        from cauchy_diagonal_packet_reference import LABELS, emission, injection, step
        from rate_design_coset import mul, power as gf_power
        for phase in range(4):
            for lane in range(16):
                unit = tuple(int(i == lane) for i in range(16))
                column = injection(unit, phase)
                self.assertTrue(all(column))
                paired = emission(column, phase)
                self.assertEqual(paired, tuple(gf_power(u ^ LABELS[phase][lane], 3)
                                               for u in LABELS[phase]))
                for scalar in (1, 2, 255):
                    values = tuple(mul(scalar, bit) for bit in unit)
                    state = tuple(mul(scalar, value) for value in column)
                    outputs, following = step(values, state, (1, 1, 1, 1), phase)
                    self.assertEqual(following, (0, 0, 0, 0))
                    self.assertTrue(all(outputs))
                    self.assertEqual(outputs[lane], scalar)

    def test_frozen_certificates(self):
        expected = {16:107, 18:105, 20:103}
        for k_log, margin in expected.items():
            result = certificate.certify(k_log)
            self.assertEqual(result['class_margin_bits_lower'], margin)
            self.assertLessEqual(F(result['failure_upper']), F(2)**(-margin))
            self.assertEqual(result['active_outer_rows'], 1)
            self.assertFalse(result['whole_code'])


if __name__ == '__main__':
    unittest.main()

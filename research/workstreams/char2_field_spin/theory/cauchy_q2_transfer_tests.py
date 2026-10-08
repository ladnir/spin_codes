"""Independent finite-route, source-count, and directed Q2 checks."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import math
import unittest

import numpy as np

import cauchy_q2_transfer as certificate


def multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]


def route_average(local, rows, width, active):
    total = [[F(0), F(0)], [F(0), F(0)]]
    for slots in combinations(range(rows), active):
        chosen = set(slots)
        matrix = [[F(1), F(0)], [F(0), F(1)]]
        for start in range(0, rows, width):
            count = sum(slot in chosen for slot in range(start, start+width))
            matrix = multiply(matrix, local[count])
        for i in range(2):
            for j in range(2):
                total[i][j] += matrix[i][j]/math.comb(rows, active)
    return total


def enumerate_supports(regions, length):
    result = [[[F(0), F(0)] for _ in range(length+1)]
              for _ in range(length+1)]
    for first in product((0, 1), repeat=length):
        for second in product((0, 1), repeat=length):
            matrix = [[F(1), F(0)], [F(0), F(1)]]
            for a, b in zip(first, second):
                matrix = multiply(matrix, regions[a+b])
            for end in range(2):
                result[sum(first)][sum(second)][end] += matrix[0][end]
    return result


class Q2TransferTests(unittest.TestCase):
    def assert_encloses(self, interval, exact):
        self.assertLessEqual(F(interval.lo), exact)
        self.assertGreaterEqual(F(interval.hi), exact)

    def test_distinct_slot_average_exhaustive(self):
        local = [[[F(1, 2), F(1, 3)], [F(1, 5), F(1, 7)]],
                 [[F(2, 3), F(3, 5)], [F(4, 7), F(5, 11)]],
                 [[F(3, 7), F(4, 9)], [F(5, 13), F(6, 17)]]]
        for rows, width in ((4, 2), (6, 2), (6, 3), (9, 3)):
            actual = certificate.regional_coefficients(local, rows, width)
            for active in range(3):
                self.assertEqual(actual[active], route_average(local, rows, width, active))

    def test_same_packet_collision_mass(self):
        for rows in (32, 48, 64):
            total = sum(a//16 == b//16 for a, b in combinations(range(rows), 2))
            self.assertEqual(F(total, math.comb(rows, 2)), F(15, rows-1))

    def test_two_markers_against_support_enumeration(self):
        regions = [[[F(1), F(0)], [F(0), F(1, 2)]],
                   [[F(0), F(2, 3)], [F(1, 5), F(3, 7)]],
                   [[F(0), F(4, 9)], [F(2, 11), F(5, 13)]]]
        expected = enumerate_supports(regions, 4)
        actual = certificate.two_row_polynomial(regions, 4)
        self.assertEqual(actual, expected)
        for a in range(5):
            for b in range(5):
                self.assertEqual(actual[a][b], actual[b][a])

    def test_one_global_projective_quotient_gf4(self):
        # Evaluation of a+b*x at the four GF4 points: a [4,2] MDS code.
        def mul(a, b):
            value = 0
            while b:
                if b & 1:
                    value ^= a
                a <<= 1
                if a & 4:
                    a ^= 7
                b >>= 1
            return value

        words = [tuple(a ^ mul(b, x) for x in range(4))
                 for a, b in product(range(4), repeat=2) if a or b]
        counts = Counter(tuple(i for i, x in enumerate(word) if x) for word in words)
        self.assertEqual(len(words), 15)
        per_support = {}
        for support, count in counts.items():
            n = certificate.exact_mds_projective_support(len(support), parity=2, q=4)
            self.assertEqual(count, 3*n)
            per_support[support] = n
        actual = Counter((tuple(i for i, x in enumerate(a) if x),
                          tuple(i for i, x in enumerate(b) if x))
                         for a, b in product(words, repeat=2))
        self.assertEqual(sum(actual.values())//3, 75)
        for (a, b), count in actual.items():
            self.assertEqual(F(count, 3), 3*per_support[a]*per_support[b])

    def test_packet_outward_enclosures(self):
        for mark in (F(1), F(2, 3), F(24993, 25000)):
            exact = certificate.packet_matrices(mark)
            directed = certificate.directed_packets(mark)
            for c, i, j in product(range(3), range(2), range(2)):
                self.assert_encloses(directed[c][i][j], exact[c][i][j])

    def test_complete_small_expression(self):
        mark, length, parity = F(2, 3), 4, 2
        for rows in (16, 32):
            local = certificate.packet_matrices(mark)
            regions = [route_average(local, rows, 16, c) for c in range(3)]
            polynomial = enumerate_supports(regions, length)
            counts = [sum((-1)**j*math.comb(w, j)*256**max(0, w-parity-j)
                          for j in range(w+1))//255 if w > parity else 0
                      for w in range(length+1)]
            exact = sum(sum(polynomial[a][b])*counts[a]*counts[b]
                        for a in range(length+1) for b in range(length+1))
            threshold = (length*rows+9)//10-1
            exact *= math.comb(rows, 2)*255/mark**threshold
            self.assert_encloses(certificate.directed_expression(rows, mark, length, parity), exact)

    def test_float_log_path_matches_positive_reference(self):
        mark, rows = F(2, 3), 32
        regions = certificate.regional_coefficients(certificate.packet_matrices(mark), rows)
        log_regions = certificate.log_regions(float(mark), rows)
        np.testing.assert_allclose(np.exp(log_regions), np.asarray(regions, dtype=float), rtol=1e-13)
        exact = certificate.two_row_polynomial(regions, 4)
        log_result = certificate.log_two_row_polynomial(log_regions, 4)
        np.testing.assert_allclose(np.exp(log_result), np.asarray(exact, dtype=float), rtol=1e-13)

    def test_frozen_certificates(self):
        for k_log, margin in ((16, 65), (18, 72), (20, 77)):
            result = certificate.certify(k_log)
            self.assertEqual(result['class_margin_bits_lower'], margin)
            self.assertLessEqual(F(result['failure_upper']), F(2)**(-margin))
            self.assertEqual(result['active_rows'], 2)
            self.assertFalse(result['whole_code'])

    def test_union_source_recurrence_by_enumeration(self):
        length, parity, theta = 4, 2, F(2, 7)
        for active in (0, 1, 2, 3):
            actual = certificate.log_union_source(active, math.log(theta), length, parity)
            for size in range(length+1):
                masks = [mask for mask in product((0, 1), repeat=size) if sum(mask) > parity]
                expected = F(0)
                for chosen in product(masks, repeat=active):
                    if all(any(mask[j] for mask in chosen) for j in range(size)):
                        term = F(1)
                        for mask in chosen:
                            weight = sum(mask)
                            term *= theta**weight*certificate.exact_mds_projective_support(weight, parity)
                        expected += term
                if expected:
                    self.assertAlmostEqual(actual[size], math.log(expected), places=11)
                else:
                    self.assertEqual(actual[size], -math.inf)

    def test_union_envelope_contains_exact_two_row_expression(self):
        # Compare only the positive source/route moment on a small b=4,p=2
        # example; geometry and global factors are common to both methods.
        length, parity, mark, rows = 4, 2, F(2, 3), 32
        regions = certificate.log_all_regions(float(mark), rows, 2)
        exact_regions = certificate.regional_coefficients(certificate.packet_matrices(mark), rows)
        polynomial = certificate.two_row_polynomial(exact_regions, length)
        counts = [certificate.exact_mds_projective_support(w, parity) for w in range(length+1)]
        exact = sum(sum(polynomial[a][b])*counts[a]*counts[b]
                    for a in range(length+1) for b in range(length+1))
        for theta in (F(1, 17), F(2, 7), F(2)):
            logtheta = math.log(theta)
            nonempty = np.max(regions[1:]-np.arange(1, 3)[:, None, None]*logtheta, axis=0)
            source = certificate.log_union_source(2, logtheta, length, parity)
            routed = certificate.log_region_union_polynomial(regions[0], nonempty, length)
            upper = float(certificate.logsumexp(source+routed))
            self.assertGreaterEqual(upper+1e-12, math.log(exact))

    def test_word_moment_obstruction_is_not_a_failure_probability(self):
        witness = certificate.word_moment_obstruction()
        self.assertEqual(witness['window'], 4752)
        self.assertEqual(witness['dimension'], 228)
        self.assertEqual(witness['active_rows_at_most'], 16)
        self.assertGreater(witness['expected_projective_lines_lower'], 2**111)


if __name__ == '__main__':
    unittest.main()

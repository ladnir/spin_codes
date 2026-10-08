"""Partial-certificate replay, interval checks, and first-moment obstruction."""

from copy import deepcopy
from fractions import Fraction
from itertools import product
import unittest

from subfield_certificate import (Interval, certify_sparse, interval_regions,
                                  verify_sparse, window_obstruction)
from batched_transfer import multiply


def gf4_rank(matrix):
    matrix = [list(row) for row in matrix]
    if not matrix:
        return 0
    row = 0
    for column in range(len(matrix[0])):
        pivot = next((i for i in range(row, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        inverse = next(x for x in range(1, 4) if multiply(matrix[row][column], x, 7) == 1)
        matrix[row] = [multiply(x, inverse, 7) for x in matrix[row]]
        for i in range(len(matrix)):
            if i != row and matrix[i][column]:
                factor = matrix[i][column]
                matrix[i] = [x ^ multiply(factor, y, 7) for x, y in zip(matrix[i], matrix[row])]
        row += 1
        if row == len(matrix):
            break
    return row


def gf4_dot(left, right):
    out = 0
    for x, y in zip(left, right):
        out ^= multiply(x, y, 7)
    return out


def sl2_gf4():
    return [[[a, b], [c, d]] for a, b, c, d in product(range(4), repeat=4)
            if multiply(a, d, 7) ^ multiply(b, c, 7) == 1]


def gf4_product(left, right):
    columns = list(zip(*right))
    return [[gf4_dot(row, column) for column in columns] for row in left]


class SubfieldCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = certify_sparse()

    def test_sparse_certificate_replays(self):
        self.assertTrue(verify_sparse(self.record))
        self.assertFalse(self.record['complete'])
        self.assertEqual(self.record['covered_active_rows'], [1, 2])
        self.assertGreaterEqual(self.record['certified_integer_bits'], 43)

    def test_whole_certificate_claim_rejected(self):
        changed = deepcopy(self.record)
        changed['complete'] = True
        with self.assertRaises(ValueError):
            verify_sparse(changed)

    def test_missing_class_rejected(self):
        changed = deepcopy(self.record)
        changed['entries'].pop()
        with self.assertRaises(ValueError):
            verify_sparse(changed)

    def test_understated_bound_rejected(self):
        changed = deepcopy(self.record)
        changed['entries'][1]['failure_upper'] = '1e-20'
        with self.assertRaises(ValueError):
            verify_sparse(changed)

    def test_exaggerated_margin_rejected(self):
        changed = deepcopy(self.record)
        changed['certified_integer_bits'] = 80
        with self.assertRaises(ValueError):
            verify_sparse(changed)

    def test_interval_arithmetic_encloses_exact(self):
        for a in [Fraction(1, 3), Fraction(-7, 11), Fraction(511, 512)]:
            ia = Interval.rational(a.numerator, a.denominator)
            for b in [Fraction(2, 7), Fraction(11, 2)]:
                ib = Interval.rational(b.numerator, b.denominator)
                for result, exact in [(ia+ib, a+b), (ia-ib, a-b), (ia*ib, a*b),
                                      (ia/ib, a/b), (ia**9, a**9)]:
                    self.assertLessEqual(Fraction(result.lo), exact)
                    self.assertGreaterEqual(Fraction(result.hi), exact)

    def test_region_normalization_at_one(self):
        for matrix in interval_regions(Interval.integer(1), 2):
            for row in matrix:
                total = row[0]+row[1]
                self.assertLessEqual(total.lo, 1)
                self.assertGreaterEqual(total.hi, 1)

    def test_first_moment_obstruction(self):
        result = window_obstruction()
        self.assertEqual(result['subspace_dimension_lower'], 226)
        self.assertEqual(result['maximum_output_weight'], 13056)
        self.assertGreater(result['expected_bad_lines_lower'], 2**82)
        # The rare routing event is not itself a likely code failure.
        self.assertLess(result['route_event_probability'], Fraction(1, 2**1700))

    def test_diagonal_rank_potential(self):
        forms_list = [((1, 0), (0, 1), (1, 1)), ((1, 0), (1, 0), (1, 0)),
                      ((1, 0), (2, 0), (0, 1)), ((1, 1), (1, 2), (1, 3))]
        matrices = [((1, 0, 1), (0, 1, 1)), ((1, 1, 1), (1, 1, 1)),
                    ((1, 2, 3), (0, 1, 2)), ((0, 1, 1), (1, 0, 2)),
                    ((1, 0, 0), (0, 1, 0)), ((0, 0, 0), (0, 0, 0))]
        messages = list(product(range(4), repeat=2))
        for tau in [2, 7]:
            g = Fraction(1, 3)+Fraction(2, 3*tau)
            for forms in forms_list:
                alive_inputs = [m for m in messages if all(gf4_dot(f, m) for f in forms)]
                for matrix in matrices:
                    total = Fraction(0)
                    for coefficients in product(range(1, 4), repeat=3):
                        scaled = [[multiply(c, value, 7) for value in f]
                                  for c, f in zip(coefficients, forms)]
                        equations = gf4_product(matrix, scaled)
                        if any(all(gf4_dot(eq, m) == 0 for eq in equations) for m in alive_inputs):
                            total += Fraction(1, tau**gf4_rank(equations))
                    self.assertLessEqual(total/27, g**gf4_rank(matrix))

    def test_sl2_return_rank_factors(self):
        matrices = sl2_gf4()
        self.assertEqual(len(matrices), 60)
        for h, forms in [(1, [[1, 0], [0, 0]]), (2, [[1, 0], [0, 1]])]:
            images = [gf4_product(m, forms) for m in matrices]
            for target_entries in product(range(4), repeat=4):
                target = [target_entries[:2], target_entries[2:]]
                ranks = [gf4_rank([[x ^ y for x, y in zip(row, other)]
                                   for row, other in zip(image, target)]) for image in images]
                for tau in [2, 7]:
                    actual = sum((Fraction(1, tau**r) for r in ranks), Fraction(0))/60
                    if h == 1:
                        upper = Fraction(1, tau)+(1-Fraction(1, tau))/15
                    else:
                        upper = (Fraction(1, tau*tau)+(Fraction(1, tau)-Fraction(1, tau*tau))/3
                                 +(1-Fraction(1, tau))/60)
                    self.assertLessEqual(actual, upper)

    def test_empty_input_rank_operator(self):
        matrices = sl2_gf4()
        q, width, tau = 4, 4, 3
        z = Fraction(2, 3)
        messages = list(product(range(4), repeat=2))
        for h, forms in [(1, [[1, 0], [0, 0]]), (2, [[1, 0], [0, 1]])]:
            totals = [Fraction(0)]*3
            for matrix in matrices:
                refreshed = gf4_product(matrix, forms)
                emitted = gf4_product([[1, j] for j in range(width)], refreshed)
                for mask in range(1 << width):
                    zero_rows = [emitted[j] for j in range(width) if (mask >> j) & 1]
                    nonzero_rows = [emitted[j] for j in range(width) if not ((mask >> j) & 1)]
                    alive = any(any(gf4_dot(f, m) for f in forms)
                                and all(gf4_dot(f, m) == 0 for f in zero_rows)
                                and all(gf4_dot(f, m) != 0 for f in nonzero_rows)
                                for m in messages)
                    if not alive:
                        continue
                    rank = gf4_rank(zero_rows)
                    next_rank = gf4_rank(zero_rows+forms)-rank
                    totals[next_rank] += z**len(nonzero_rows)/tau**rank/len(matrices)
            if h == 1:
                expected = (1-Fraction(width, q+1))*z**width+Fraction(width, q+1)*z**(width-1)
                self.assertEqual(totals, [0, expected, 0])
            else:
                self.assertEqual(totals, [0, Fraction(width, tau)*z**(width-1), z**width])


if __name__ == '__main__':
    unittest.main()

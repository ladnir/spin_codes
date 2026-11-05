"""Exact rational and finite-scope tests for fresh directed mixed replay."""
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest

from flint import arb, ctx, fmpq

from hd_xor_fixed_mask_floor import comparison_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_mask_certificate import (EXPECTED_SCOPE, validate_scope,
    offset_resource_arb, mixed_vectors_arb, fixed_branch_caps_arb,
    column_coefficients_arb, positive_powers_arb, bivariate_tail_arb, receipt)


def directed(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def multiply(left, right):
    return tuple(sum(left[2*i+k]*right[2*k+j] for k in range(2))
                 for i in range(2) for j in range(2))


class MixedMaskCertificateTests(unittest.TestCase):
    def assertContainsRational(self, interval, value):
        value = F(value)
        rational = fmpq(value.numerator, value.denominator)
        self.assertLessEqual(interval.lower().fmpq(), rational)
        self.assertGreaterEqual(interval.upper().fmpq(), rational)

    def assertUpper(self, upper, value):
        value = F(value)
        self.assertGreaterEqual(upper.lower().fmpq(), fmpq(value.numerator, value.denominator))

    def test_noncommuting_offset_words_and_reachability_pruning(self):
        exact = [[tuple(F(1+7*i+3*m+j, 29+j) for j in range(4)) for m in range(4)]
                 for i in range(4)]
        with ctx.workprec(192):
            tables = [[tuple(directed(value) for value in matrix) for matrix in table] for table in exact]
            potential = F(3, 11)
            full = offset_resource_arb(tables, directed(potential), 0, 2, range(16))
            pruned = offset_resource_arb(tables, directed(potential), 0, 2, (0, 3, 5, 15), (1, 4, 7))
            for mask in pruned:
                self.assertEqual(set(pruned[mask]), {1, 4, 7})
                for target in pruned[mask]:
                    self.assertEqual(pruned[mask][target], full[mask][target])
            for mask, word in product(range(16), product(range(3), repeat=4)):
                matrix = (F(1), F(0), F(0), F(1))
                for i, occupancy in enumerate(word):
                    matrix = multiply(matrix, exact[i][occupancy+((mask >> i) & 1)])
                first, second = matrix[0]+matrix[1]*potential, matrix[2]+matrix[3]*potential
                self.assertUpper(full[mask][sum(word)][0], first)
                self.assertUpper(full[mask][sum(word)][1], second)

    def test_fixed_occupancy_is_exact_word_upper(self):
        with ctx.workprec(192):
            exact = [[tuple(F(2+i+m+j, 7+j) for j in range(4)) for m in range(3)] for i in range(4)]
            tables = [[tuple(directed(value) for value in matrix) for matrix in table] for table in exact]
            actual = offset_resource_arb(tables, directed(F(1, 4)), 1, 1, range(16))
            for mask in range(16):
                matrix = (F(1), F(0), F(0), F(1))
                for i in range(4):
                    matrix = multiply(matrix, exact[i][1+((mask >> i) & 1)])
                for value, expected in zip(actual[mask][0], (matrix[0]+matrix[1]/4, matrix[2]+matrix[3]/4)):
                    self.assertUpper(value, expected)
                    self.assertLess(abs(float(value/directed(expected))-1), 1e-14)

    def test_cancellation_probabilities_and_dense_mask_normalization(self):
        with ctx.workprec(192):
            identity = (arb(1), arb(0), arb(0), arb(1))
            tables = [[identity]*4 for _ in range(8)]
            potential = directed(F(1, 8))
            masks = [mask for mask, count in enumerate(comparison_counts(195)) if count]
            vectors = mixed_vectors_arb(tables, potential, masks, groups=2, mu=2)
            for table in vectors.values():
                for first, second in table.values():
                    self.assertContainsRational(first, 1)
                    self.assertContainsRational(second, F(1, 8))
            empty, nonempty = fixed_branch_caps_arb(vectors, 195, potential)
            for kind in categories(2):
                self.assertUpper(empty[kind], F(1, 256**4))
                self.assertUpper(nonempty[kind], 1-F(1, 256**4))
                self.assertLess(abs(float(empty[kind]*256**4)-1), 1e-14)

    def test_dense_mixture_precedes_the_quotient(self):
        with ctx.workprec(192):
            exact = {mask: (F(1+mask**2), F(1, 4)+((mask+7) % 16)**2) for mask in range(16)}
            vectors = {(0, 0): {mask: tuple(directed(x) for x in value) for mask, value in exact.items()}}
            empty, nonempty = fixed_branch_caps_arb(vectors, 7, directed(F(1, 4)), group=4, q=4)
            counts = comparison_counts(7, 4, 4)
            means = [sum(F(counts[mask], 64)*exact[mask][j] for mask in range(1, 16)) for j in range(2)]
            value = max(means[0], 4*means[1])
            self.assertUpper(nonempty[0, 0], value)
            self.assertLess(abs(float(nonempty[0, 0]/directed(value))-1), 1e-14)
            premature = sum(F(counts[mask], 64)*max(exact[mask][0], 4*exact[mask][1]) for mask in range(1, 16))
            self.assertLess(value, premature)

    def test_positive_bivariate_tail_matches_rational_word_enumeration(self):
        empty, nonempty = [F(2, 3), F(1, 5), F(3, 7)], [F(1, 11), F(2, 5)]
        choices = [(flag, degree, value) for flag, poly in enumerate((empty, nonempty)) for degree, value in enumerate(poly)]
        expected = {}
        for word in product(choices, repeat=4):
            x, y, weight = 0, 0, F(1)
            for flag, degree, value in word:
                x += degree
                y += flag
                weight *= value
            expected[x, y] = expected.get((x, y), 0)+weight
        with ctx.workprec(192):
            first, second = list(map(directed, empty)), list(map(directed, nonempty))
            for xf, yf in product(range(5), range(5)):
                actual, by_y = bivariate_tail_arb(first, second, columns=4, x_floor=xf, y_floor=yf)
                self.assertContainsRational(actual, sum(value for (x, y), value in expected.items() if x >= xf and y >= yf))
                for y, value in enumerate(by_y, yf):
                    self.assertContainsRational(value, sum(v for (x, yy), v in expected.items() if x >= xf and yy == y))

    def test_empty_scaling_and_independent_sparse_normalizer(self):
        with ctx.workprec(192):
            base = column_coefficients_arb({kind: arb(1) for kind in categories(1)}, groups=1, q=4)
            probability = F(1, 64)
            empty = [value*directed(probability) for value in base]
            nonempty = [value*directed(1-probability) for value in base]
            tail, _ = bivariate_tail_arb(empty, nonempty, 3, 4, 2,
                empty_powers=positive_powers_arb(base, 1), empty_scale=directed(probability))
            sparse = sum(comb(6, degree)*3**degree for degree in range(4, 7))
            dense = sum(comb(3, degree)*(1-probability)**degree*probability**(3-degree) for degree in range(2, 4))
            self.assertContainsRational(tail, sparse*dense)
            self.assertContainsRational(tail/3**2, F(sparse, 3**2)*dense)

    def test_scope_and_masks_are_checked_before_expensive_work(self):
        validate_scope()
        for key, value in (('cohort_counts', (1, 15, 1)), ('field_order', 255), ('forbidden_weight', 734004),
                           ('total_columns', 31), ('source_matching_counts', (0, 8)), ('mu_range', (1, 16)),
                           ('sparse_support_floor', 543), ('dense_nonempty_column_floor', 16)):
            scope = dict(EXPECTED_SCOPE, **{key: value})
            with self.assertRaises(ValueError):
                receipt(scope=scope)
        for masks in ((3,), (195, 195), (195.0,), ()):
            with self.assertRaises(ValueError):
                receipt(active_masks=masks)
        with self.assertRaises(ValueError):
            receipt(precision=64)


if __name__ == '__main__':
    unittest.main()

"""Tiny exact/positive checks; no provider or production common-table build."""
from collections import Counter
from fractions import Fraction
from itertools import product
import unittest

import numpy as np

import hd_lambda_shear_shared_ratio as ratio
from hd_lambda_scheduled_relations import row_forms
from hd_quaternary_outer_relations import multiply_gf4
from hd_pair_state_tail_general import entrywise_power


def logs(values):
    with np.errstate(divide='ignore'):
        return np.log2(np.asarray(values, dtype=float))


def mask(values):
    return sum(bool(x) << i for i, x in enumerate(values))


def monomial(left, right, marks):
    answer = Fraction(1)
    rows = len(marks)//2
    for side, support in enumerate((left, right)):
        for i in range(rows):
            if support >> i & 1:
                answer *= marks[side*rows+i]
    return answer


def shear(left, right):
    return (tuple(a ^ b if i % 2 == 0 else a for i, (a, b) in enumerate(zip(left, right))),
            tuple(b if i % 2 == 0 else a ^ b for i, (a, b) in enumerate(zip(left, right))))


def evaluate(forms, values):
    output = []
    for form in forms:
        value = 0
        for i, x in enumerate(values):
            value ^= multiply_gf4(form >> (2*i) & 3, x)
        output.append(value)
    return tuple(output)


def tiny_words(arity):
    # Invertible two-row map at arity2; both coefficients nonzero at arity1.
    forms = (2, 3) if arity == 1 else (2 | 3 << 2, 3 | 2 << 2)
    return tuple((sum(bool(x) for x in values), evaluate(forms, values))
                 for values in product(range(4), repeat=arity))


def tiny_histogram(words):
    counter = Counter((weight, mask(word)) for weight, word in words)
    return tuple((weight, support, count) for (weight, support), count in sorted(counter.items()))


class SharedRatioTests(unittest.TestCase):
    def test_physical_or_xor_masks(self):
        for rows in (1, 2, 4):
            for first, second in product(range(1 << rows), repeat=2):
                got = ratio.output_masks(np.array([first]), np.array([second]), rows=rows)
                union, difference = first | second, first ^ second
                expected = [0]*4
                for i in range(rows):
                    expected[0] |= ((union if i % 2 == 0 else first) >> i & 1) << i
                    expected[1] |= ((second if i % 2 == 0 else union) >> i & 1) << i
                    expected[2] |= ((difference if i % 2 == 0 else first) >> i & 1) << i
                    expected[3] |= ((second if i % 2 == 0 else difference) >> i & 1) << i
                self.assertEqual(tuple(int(x[0]) for x in got), tuple(expected))

    def test_exact_GF4_ratio_average_including_zero_marks(self):
        mark_sets = ((Fraction(0), Fraction(1, 3), Fraction(1), Fraction(2, 3)),
                     (Fraction(1, 2), Fraction(2, 3), Fraction(3, 4), Fraction(1, 4)))
        for left, right in product(product(range(4), repeat=2), repeat=2):
            supports = ratio.output_masks(np.array([mask(left)]), np.array([mask(right)]), rows=2)
            for marks in mark_sets:
                exact = Fraction(0)
                for rho in (1, 2, 3):
                    scaled = tuple(multiply_gf4(rho, value) for value in right)
                    u, v = shear(left, scaled)
                    exact += monomial(mask(u), mask(v), marks)/3
                upper = (2*monomial(int(supports[0][0]), int(supports[1][0]), marks)
                         +monomial(int(supports[2][0]), int(supports[3][0]), marks))/3
                self.assertLessEqual(exact, upper)

    def test_minkowski_probabilities_are_outside_roots(self):
        matrices = np.array([[[2., 1.], [.2, .5]], [[.5, .1], [1., 2.]]])
        prepared = ratio.prepare_atoms(logs([.1, .9]), logs(matrices), 3, rows=2)
        first, second = np.array([3, 1]), np.array([3, 2])
        masks = ratio.output_masks(first, second, rows=2)
        jo = ratio.branch_log_matrices(masks[0], masks[1], prepared.table)
        jx = ratio.branch_log_matrices(masks[2], masks[3], prepared.table)
        expected = Fraction(2, 3)*np.exp2(jo/3)+Fraction(1, 3)*np.exp2(jx/3)
        got = ratio.mixture_log_norms(first, second, prepared, q=4)
        np.testing.assert_allclose(np.exp2(got), np.asarray(expected, float), rtol=3e-15)

    def test_true_GF4_a1_a2_moment_polynomials_noncommuting_measure(self):
        atoms = np.array([.2, .8])
        matrices = np.array([[[.7, .2], [.1, .6]], [[.2, .5], [.3, .1]]])
        self.assertFalse(np.allclose(matrices[0] @ matrices[1], matrices[1] @ matrices[0]))
        for arity, groups in ((1, 3), (2, 2)):
            words = tiny_words(arity)
            histogram = tiny_histogram(words)
            joint = Counter()
            for (wa, left), (wb, right) in product(words, repeat=2):
                u, v = shear(left, right)
                joint[wa+wb, mask(u), mask(v)] += 1
            prepared = ratio.prepare_atoms(logs(atoms), logs(matrices), groups, rows=2)
            coefficients = ratio.pair_coefficients(histogram, histogram, prepared, q=4, chunk_size=3)
            upper = np.exp2(entrywise_power(coefficients, groups))
            exact = np.zeros_like(upper)
            for choices in product(range(2), repeat=4):
                marks = atoms[list(choices)]
                polynomial = np.zeros(2*arity+1)
                for (weight, left, right), count in joint.items():
                    value = float(count)
                    for support, offset in ((left, 0), (right, 2)):
                        for i in range(2):
                            if support >> i & 1:
                                value *= marks[offset+i]
                    polynomial[weight] += value
                powered = np.array([1.])
                for _ in range(groups):
                    powered = np.convolve(powered, polynomial)
                matrix = np.eye(2)
                for choice in choices:
                    matrix = matrix @ matrices[choice]
                exact += powered[:, None, None]*matrix
            self.assertTrue(np.all(upper+1e-10*np.maximum(1, exact) >= exact))
            np.testing.assert_allclose(upper[0], exact[0], rtol=4e-15)

    def test_actual_scheduled_a1_a2_scalar_sanity(self):
        marks = tuple(Fraction((i % 5)+1, 6) for i in range(32))
        for source, column in (((3,), 0), ((0, 9), 8)):
            a = len(source)
            left_forms = row_forms(source, ratio.SCHEDULE, column)
            right_forms = row_forms(source, ratio.SCHEDULE, column+1)
            exact = [Fraction(0)]*(2*a+1)
            for left, right in product(product(range(4), repeat=a), repeat=2):
                u, v = shear(evaluate(left_forms, left), evaluate(right_forms, right))
                exact[sum(bool(x) for x in left+right)] += monomial(mask(u), mask(v), marks)
            upper = [Fraction(0)]*(2*a+1)
            first = ratio.column_histogram(source, column, q=4)
            second = ratio.column_histogram(source, column+1, q=4)
            for (wa, sa, na), (wb, sb, nb) in product(first, second):
                m = ratio.output_masks(np.array([sa]), np.array([sb]))
                upper[wa+wb] += na*nb*(2*monomial(int(m[0][0]), int(m[1][0]), marks)
                    +monomial(int(m[2][0]), int(m[3][0]), marks))/3
            self.assertTrue(all(x <= y for x, y in zip(exact, upper)))

    def test_a4_histograms_exact_mass_and_bounded_products(self):
        for source in ((0, 1, 2, 3), (0, 1, 4, 8)):
            first = ratio.column_histogram(source, 0)
            second = ratio.column_histogram(source, 1)
            self.assertLessEqual(len(first)*len(second), 529**2)
            self.assertEqual(sum(n for _, _, n in first), 256**4)
            self.assertEqual(sum(na*nb for _, _, na in first for _, _, nb in second), 256**8)
            self.assertEqual(first[0], (0, 0, 1))

    def test_log_scale_and_chunk_invariance(self):
        matrixlogs = logs([[[1., 2.], [3., 1.]], [[2., 1.], [1., 2.]]])
        base = ratio.prepare_atoms(logs([0., 1.]), matrixlogs, 2, rows=2)
        tiny = ratio.prepare_atoms(logs([0., 1.]), matrixlogs-1000, 2, rows=2)
        histogram = ((0, 0, 1), (1, 3, 3))
        expected = ratio.pair_coefficients(histogram, histogram, base, q=4)
        actual = ratio.pair_coefficients(histogram, histogram, tiny, q=4, chunk_size=1)
        np.testing.assert_allclose(actual, expected-2000, rtol=0, atol=2e-12)

    def test_dependencies_scope_and_guards(self):
        self.assertIn('hd_lambda_shear_shared_ratio.py', ratio.dependencies())
        matrices = np.zeros((1, 2, 2))
        with self.assertRaises(ValueError):
            ratio.prepare_atoms([.01], matrices, 2, rows=2)
        prepared = ratio.prepare_atoms([0.], matrices, 2, rows=2)
        with self.assertRaises(ValueError):
            ratio.pair_coefficients(((0, 0, 1),), ((0, 0, 1),), prepared, q=8)
        with self.assertRaises(ValueError):
            ratio.pair_coefficients(((0, 0, 0),), ((0, 0, 1),), prepared)
        with self.assertRaises(ValueError):
            ratio.pair_polynomial((0,), 1, prepared)
        with self.assertRaises(ValueError):
            ratio.column_histogram((0, 1, 2, 3, 4), 0)
        for name in ('density_or_renyi_upper', 'whole_code_certificate',
                     'coefficientwise_support_domination', 'rs_normalization_included'):
            self.assertFalse(ratio.metadata()[name])


if __name__ == '__main__':
    unittest.main()

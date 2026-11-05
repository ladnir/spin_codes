"""Tiny exact local tests only; no all-mask preparation or provider run."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest

import numpy as np

from hd_bit_syndrome_cross_fibers import cross_fiber, inside_support_counts
from hd_bit_syndrome_flat_geometry import outside_image_support_counts
import hd_bit_syndrome_scalar_allarity as scalar


def mul(x, y, q):
    result = 0
    modulus = {2: 3, 4: 7}[q]
    while y:
        if y & 1:
            result ^= x
        y >>= 1
        x <<= 1
        if x & q:
            x ^= modulus
    return result


def dot(x, y, q):
    result = 0
    for a, b in zip(x, y):
        result ^= mul(a, b, q)
    return result


def histogram(words, length):
    counts = Counter(sum(x != 0 for x in word) for word in words)
    return tuple(counts[j] for j in range(length+1))


def evaluate(coefficients, t):
    return sum(value*t**j for j, value in enumerate(coefficients))


class ScalarAllarityTests(unittest.TestCase):
    def test_macwilliams_against_small_field_enumeration(self):
        for q in (2, 4):
            rows = ((1, 0, 1, q-1), (0, 1, q-1, 1))
            code = {tuple(dot(values, column, q) for column in zip(*rows))
                    for values in product(range(q), repeat=2)}
            kernel = {word for word in product(range(q), repeat=4)
                      if all(dot(row, word, q) == 0 for row in rows)}
            dual = histogram(code, 4)
            self.assertEqual(scalar.macwilliams_kernel_counts(dual, q), histogram(kernel, 4))
            self.assertEqual(scalar.macwilliams_kernel_counts((1, 0, 0), q),
                             (1, 2*(q-1), (q-1)**2))

    def test_translated_fiber_and_positive_nonempty_toy(self):
        # H=I+11^T in characteristic2, length4: a symmetric involution.
        # Exhaust every input for each of15 labelled source sets over GF4.
        q, size = 4, 4
        h = tuple(tuple(1 ^ int(i == j) for j in range(size)) for i in range(size))
        for mask in range(1, 1 << size):
            a = tuple(i for i in range(size) if mask >> i & 1)
            o = tuple(i for i in range(size) if not mask >> i & 1)
            fibers = {}
            actual = Counter()
            for x in product(range(q), repeat=len(a)):
                y = tuple(dot(tuple(row[j] for j in a), x, q) for row in h)
                fibers.setdefault(tuple(y[i] for i in o), []).append(tuple(y[i] for i in a))
                actual[sum(value != 0 for value in y)] += 1
            c = histogram(fibers, len(o))
            z = histogram(fibers[(0,)*len(o)], len(a))
            nonempty = scalar.positive_nonempty_coefficients(c, z)
            for t in (F(0), F(1, 7), F(1, 2), F(1)):
                homogeneous = evaluate(z, t)
                for fiber in fibers.values():
                    self.assertLessEqual(evaluate(histogram(fiber, len(a)), t), homogeneous)
                actual_nonempty = sum(count*t**w for w, count in actual.items() if w)
                self.assertGreaterEqual(evaluate(nonempty, t), actual_nonempty)
                self.assertEqual(nonempty[0], 0)

    def test_actual_image_low_arity_and_endpoint_counts(self):
        for mask in (1, 9, 0x7f):
            supports = outside_image_support_counts(mask)
            expected = [0]*(17-mask.bit_count())
            for support, count in supports.items():
                expected[support.bit_count()] += count
            self.assertEqual(scalar.outside_weight_counts(mask), tuple(expected))
        self.assertEqual(scalar.outside_weight_counts(0), (1,)+(0,)*16)
        self.assertEqual(scalar.outside_weight_counts(0xffff), (1,))
        self.assertEqual(scalar.outside_weight_counts(0xfffe), (1, 255))

    def test_actual_eight_nine_exceptions_against_inverse_elimination(self):
        # No affine transport: these are independently labelled actual fibers.
        for mask in (0x00ff, 0x017f, 0xff80):
            source = tuple(i for i in range(16) if mask >> i & 1)
            fiber = cross_fiber(source)
            exact = inside_support_counts(fiber, max_subsets=512)
            expected = [0]*(len(source)+1)
            for support, count in exact.items():
                expected[support.bit_count()] += count
            dual = scalar.outside_weight_counts(0xffff ^ mask)
            self.assertEqual(scalar.macwilliams_kernel_counts(dual), tuple(expected))
            self.assertEqual(sum(dual), 256**fiber['cross_rank'])
        self.assertEqual(cross_fiber(tuple(range(7, 16)))['zero_outside_excess_dimension'], 1)

    def test_prepare_compression_fixed_word_and_all_arities_small_sample(self):
        masks = tuple((1 << a)-1 for a in range(1, 17))
        prepared = scalar.prepare_geometry(source_masks=masks, max_masks=16)
        self.assertFalse(prepared['metadata']['all_labelled_nonempty_rowsets'])
        self.assertEqual(prepared['source_mask_count'], 16)
        self.assertEqual(prepared['counts_by_arity'], {a: 1 for a in range(1, 17)})
        exact = scalar.group_logs(-3., prepared)
        coarse = scalar.group_logs(-3., prepared, fiber_bound='dimension')
        for a in range(1, 17):
            self.assertLessEqual(exact[a], coarse[a]+1e-11)
        at_one = scalar.group_logs(0., prepared)
        for a in range(1, 17):
            tail = sum(comb(16, j)*(256**a-1)**j for j in range(9, 17))
            self.assertAlmostEqual(at_one[a], log2(tail)-8*a*log2(255), places=10)
        self.assertTrue(all(np.isneginf(value) for value in scalar.group_logs(-np.inf, prepared).values()))

    def test_class_multiplicity_is_outside_sixteen_column_tail(self):
        prepared = scalar.prepare_geometry(source_masks=(1, 2), max_masks=2)
        self.assertEqual(len(prepared['classes']), 1)
        item = prepared['classes'][0]
        self.assertEqual(item['multiplicity'], 2)
        t = F(1, 2)
        nonempty = evaluate(item['nonempty_coefficients'], t)
        expected = 2*sum(comb(16, j)*nonempty**j for j in range(9, 17))/255**8
        self.assertAlmostEqual(scalar.group_logs(-1., prepared)[1], log2(expected), places=11)

    def test_positive_binomial_tail_exact_small_and_zero(self):
        for x in (F(1, 16), F(1), F(5)):
            expected = sum(comb(5, j)*x**j for j in range(3, 6))
            self.assertAlmostEqual(scalar.nonzero_column_tail_log(log2(x), columns=5, minimum=3),
                                   log2(expected), places=13)
        self.assertEqual(scalar.nonzero_column_tail_log(-np.inf), -np.inf)

    def test_scalar_envelope_rows_prefix_and_scope(self):
        curves = np.array([[0., -1., -2.], [-np.inf, -4., -5.],
                           [-np.inf, -6., -7.], [-1., -2., -3.]])
        envelope = scalar.scalar_envelope(curves, -2., -1.)
        potential = np.array([1., .25])
        for m in range(3):
            matrix = np.exp2(curves[:, m].reshape(2, 2))
            self.assertTrue(np.all(matrix @ potential <= 2**envelope['log2_eta']*.5**m*potential*(1+1e-14)))
        self.assertEqual(envelope['endpoint_log2'], 2.)
        partial = scalar.prepare_geometry(source_masks=(1,), max_masks=1)
        with self.assertRaises(ValueError):
            scalar.allarity_bound(partial, envelope, groups=2, source_count=3,
                                   source_scope='toy', mark=F(9, 10))

    def test_guards(self):
        for bad in (-1, True, 65536):
            with self.assertRaises(ValueError):
                scalar.outside_weight_counts(bad)
        for masks in ((), (1, 1), (0,), (True,)):
            with self.assertRaises(ValueError):
                scalar.prepare_geometry(source_masks=masks)
        with self.assertRaises(ValueError):
            scalar.prepare_geometry(source_masks=(1, 2), max_masks=1)
        for bad in (.01, np.nan):
            with self.assertRaises(ValueError):
                scalar.log_polynomial((0., 0.), bad)
        with self.assertRaises(ValueError):
            scalar.macwilliams_kernel_counts((1, 1, 0), 4)
        with self.assertRaises(ValueError):
            scalar.positive_nonempty_coefficients((2,), (1,))


if __name__ == '__main__':
    unittest.main()

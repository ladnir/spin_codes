from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, log2
import unittest

import numpy as np

from hd_bit_syndrome_counted_general import ENDPOINTS
from hd_bit_syndrome_scalar_allarity import positive_nonempty_coefficients
from hd_local_shuffle_control import column_moment, path_measure, shuffled_coefficients
from hd_local_shuffle_maclaurin_control import (
    HELPER_PINS, coefficient_classes, column_moments, combine_words,
    kernel_logs, scalar_nonempty_logs)
from hd_post_refresh_quiet_rows import moment_atoms

MEASURE = dict(log2_probability=-3., log2_atom00=-2., log2_atom01=-1., log2_atom10=-.5,
               log2_scale01=-2., log2_scale10=-4., log2_scale11=-1.)


def elementary_exact(atoms):
    result = [F(1)]+[F(0)]*len(atoms)
    for atom in atoms:
        for degree in range(len(atoms), 0, -1):
            result[degree] += result[degree-1]*atom
    return result


def elementary_table(atoms):
    return np.array([[log2(v) if v else -np.inf for v in elementary_exact(atoms)]])


class MaclaurinControlTests(unittest.TestCase):
    def test_exact_maclaurin_and_translated_fiber_toy(self):
        # Binary H=I+11^T, A={0,1}; actual enumerator 1+t^2+2t^3.
        old = shuffled_coefficients((1, 0, 1), (1, 0, 1), 1, size=4, q=2)
        active = positive_nonempty_coefficients((1, 0, 1), (1, 0, 1))
        for atoms in ((F(0), F(1, 5), F(4, 5), F(1)), (F(1, 10),)*4,
                      (F(0),)*4, (F(1),)*4):
            elementary = elementary_exact(atoms)
            mean = sum(atoms)/4
            for w in range(5):
                self.assertLessEqual(elementary[w]/comb(4, w), mean**w)
            truth = F(0)
            for x0, x1 in product(range(2), repeat=2):
                y = tuple(v ^ x0 ^ x1 for v in (x0, x1, 0, 0))
                for perm in permutations(range(4)):
                    mass = F(1)
                    for i, v in enumerate(y):
                        if v:
                            mass *= atoms[perm[i]]
                    truth += mass/24
            shuffle = sum(c*e for c, e in zip(old, elementary))
            scalar = 1+sum(c*mean**w for w, c in enumerate(active))
            self.assertLessEqual(truth, min(shuffle, scalar))
            self.assertLessEqual(truth-1, min(shuffle-1, scalar-1))
            kernels = kernel_logs(old, active, elementary_table(atoms))
            self.assertAlmostEqual(2**kernels['full'][0], float(min(shuffle, scalar)), places=13)
            self.assertAlmostEqual(2**kernels['active'][0], float(min(shuffle-1, scalar-1)), places=13)

    def test_zero_and_one_endpoints(self):
        old = shuffled_coefficients((1, 0, 1), (1, 0, 1), 1, size=4, q=2)
        active = (0, 0, 2, 0, 1)
        for atoms, full, nonempty in (((F(0),)*4, 1, 0), ((F(1),)*4, 4, 3)):
            got = kernel_logs(old, active, elementary_table(atoms))
            self.assertFalse(any(np.isnan(v).any() for v in got.values()))
            self.assertAlmostEqual(2**got['full'][0], full)
            self.assertAlmostEqual(2**got['active'][0], nonempty)

    def test_noncommuting_count_integration_matches_ordered_paths(self):
        old = shuffled_coefficients((1, 0, 1), (1, 0, 1), 1, size=4, q=2)
        active = (0, 0, 2, 0, 1)
        paths, elementary, _ = path_measure(MEASURE, length=4)
        got = column_moments(old, active, paths, elementary, 3)
        masses, _ = moment_atoms(MEASURE)
        atoms = (0., 2**MEASURE['log2_atom00'], 2**MEASURE['log2_atom01'],
                 2**MEASURE['log2_atom10'], 1.)
        expected = [np.zeros((2, 2)), np.zeros((2, 2))]
        for labels in product(range(5), repeat=4):
            edges = [ENDPOINTS[i] for i in labels]
            if any(edges[i][1] != edges[i+1][0] for i in range(3)):
                continue
            kernel = kernel_logs(old, active, elementary_table(tuple(F(atoms[i]) for i in labels)))
            for index, name in enumerate(('full', 'active')):
                expected[index][edges[0][0], edges[-1][1]] += 2**(
                    sum(masses[i] for i in labels)+3*kernel[name][0])
        for actual, wanted in zip(got, expected):
            np.testing.assert_allclose(2**actual, wanted, rtol=2e-13, atol=1e-14)
        self.assertNotAlmostEqual(expected[0][0, 1], expected[0][1, 0])
        # No pointwise repair can worsen either inherited positive branch.
        for actual, nonzero in zip(got, (False, True)):
            self.assertTrue(np.all(actual <= column_moment(old, paths, elementary, 3, nonzero=nonzero)+1e-12))

    def test_class_key_retains_both_polynomials(self):
        # Algebraic fixtures with the same old bound, different C*Z products.
        prepared = dict(classes=[
            dict(arity=2, multiplicity=3, counts_by_outputweight=(1, 0, 1), zero_fiber_counts=(1, 1), d=1),
            dict(arity=2, multiplicity=5, counts_by_outputweight=(1, 1, 2), zero_fiber_counts=(1,), d=0)])
        classes = coefficient_classes(prepared, size=4, q=2)
        self.assertEqual(len(classes), 2)
        keys = list(classes)
        self.assertEqual(keys[0][1], keys[1][1])
        self.assertNotEqual(keys[0][2], keys[1][2])
        self.assertEqual(sum(classes.values()), 8)
        prepared['classes'].append(dict(prepared['classes'][0], multiplicity=7))
        self.assertEqual(sorted(coefficient_classes(prepared, size=4, q=2).values()), [5, 10])

    def test_wholeword_norm_and_pergroup_normalization(self):
        classes = Counter({(1, (), ()): 3, (2, (), ()): 5})
        words = np.array([6., 9.])
        per, mixed = combine_words(classes, words, groups=3, prefix=11., row_dimension=2, q=4)
        b1, b2 = 3*2**(6/3)/3**2, 5*2**(9/3)/3**4
        self.assertAlmostEqual(per[1], 11+3*log2(b1), places=12)
        self.assertAlmostEqual(per[2], 11+3*log2(b2), places=12)
        self.assertAlmostEqual(mixed, 11+3*log2(b1+b2), places=12)

    def test_scope_guards_and_direct_dependency_pins(self):
        with self.assertRaises(ValueError):
            path_measure(dict(MEASURE, log2_atom01=.1))
        with self.assertRaises(ValueError):
            scalar_nonempty_logs((1, 2), np.array([-.1]))
        with self.assertRaises(ValueError):
            scalar_nonempty_logs((0, 2), np.array([.1]))
        with self.assertRaises(ValueError):
            column_moments((1, 1), (0, 1), np.zeros((1, 2, 2)), np.zeros((1, 2)), 0)
        from hashlib import sha256
        from hd_local_shuffle_maclaurin_control import DIRECTORY
        for name, digest in HELPER_PINS.items():
            self.assertEqual(sha256((DIRECTORY/name).read_bytes()).hexdigest(), digest, name)
        self.assertIn('hd_bit_syndrome_counted_general.py', HELPER_PINS)
        self.assertIn('hd_inverse_column_control.py', HELPER_PINS)


if __name__ == '__main__':
    unittest.main()

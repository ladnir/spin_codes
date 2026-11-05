"""Exact tiny controls for the positive h2 feature interface; no large bank."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, prod
import unittest

import numpy as np

from hd_xor_dense_holder import uniform_mask_count
from hd_xor_dense_holder_classes import class_moments
from hd_xor_dense_two_feature import (feature_convolution, ordered_box_vectors,
    feature_caps, positive_bivariate_tail, frozen_type_and_witness, CONTROLS)
from hd_xor_joint_columns import categories


def matrices(length=4, maximum=3):
    return [[[F(1+i+m), F(1+(i+2*m) % 3)],
             [F(2+(2*i+m) % 4), F(1+i+2*m)]]
            for i in range(length) for m in range(maximum+1)]


def exact_tables(length=4, maximum=3):
    flat = matrices(length, maximum)
    return [flat[i*(maximum+1):(i+1)*(maximum+1)] for i in range(length)]


def action(tables, occupancy, potential):
    vector = (F(1), F(potential))
    for i in reversed(range(len(tables))):
        matrix = tables[i][occupancy[i]]
        vector = tuple(sum(matrix[row][col]*vector[col] for col in range(2))
                       for row in range(2))
    return vector


def logs(tables):
    return np.log2(np.asarray(tables, dtype=float))


class DenseTwoFeatureTests(unittest.TestCase):
    def test_labelled_two_mask_convolution_matches_exact_pair_enumeration(self):
        weights = [F((i*7+3) % 11, 13) for i in range(16)]
        expected = {}
        for left, right in product(range(16), repeat=2):
            key = (int(left != 0)+int(right != 0),
                   left.bit_count()+right.bit_count(),
                   int(left == 15)+int(right == 15))
            expected[key] = expected.get(key, 0)+weights[left]*weights[right]
        expected = {key: value for key, value in expected.items() if value}
        actual = feature_convolution(weights, 2)
        self.assertEqual(actual, expected)
        self.assertEqual(sum(actual.values()), sum(weights)**2)
        self.assertEqual(actual[0, 0, 0], weights[0]**2)
        # A nonuniform labelled input must not be replaced by C(4,w)*one weight.
        self.assertNotEqual(sum(weights[1 << i] for i in range(4)), 4*weights[1])
        with self.assertRaises(ValueError):
            feature_convolution([F(1), F(-1)], 2)

    def test_holder_tensorizes_without_extra_type_multiplicity(self):
        active = (7, 15)
        nums, denominator = class_moments(active, q=4, width=4, columns=2, dimension=1)
        moments = [F(n, denominator) for n in nums]
        word_sums = {}
        for first, second in product(range(16), repeat=2):
            exact = sum(F(uniform_mask_count(4, 4, a, first)
                          * uniform_mask_count(4, 4, a, second), 3**a.bit_count())
                        for a in active)
            self.assertLessEqual(exact**2, moments[first]*moments[second])
            word_sums[first, second] = exact
        # Separate fixed source groups give an ordered A-pair sum; no 2!.
        for first, second in product(((0, 0), (3, 15), (15, 15), (7, 11)), repeat=2):
            tensor = word_sums[first]*word_sums[second]
            explicit = sum(
                F(prod(uniform_mask_count(4, 4, a, t) for t in first), 3**a.bit_count())
                * F(prod(uniform_mask_count(4, 4, b, t) for t in second), 3**b.bit_count())
                for a, b in product(active, repeat=2))
            self.assertEqual(tensor, explicit)
            self.assertLessEqual(tensor**2,
                prod(moments[t] for t in first+second))

    def test_ordered_resource_dominates_every_exact_noncommuting_word(self):
        tables, potential = exact_tables(), F(2, 3)
        bound = np.exp2(ordered_box_vectors(logs(tables), np.log2(float(potential)), 0, 2))
        for word in product(range(3), repeat=4):
            exact = action(tables, word, potential)
            for state in range(2):
                self.assertGreaterEqual(bound[sum(word), state]*(1+1e-12), float(exact[state]))
        fixed = np.exp2(ordered_box_vectors(logs(tables), np.log2(float(potential)), 1, 1))[0]
        exact = action(tables, (1,)*4, potential)
        np.testing.assert_allclose(fixed, np.array(exact, dtype=float), rtol=1e-13)
        reverse = action(list(reversed(tables)), (1,)*4, potential)
        self.assertNotEqual(exact, reverse)

    def test_features_dominate_fixed_sparse_geometry_and_whole_dense_masks(self):
        # Width four, q=4, one sparse group and two dense groups. Enumerate
        # the actual fixed constituent pair or singleton, all dense masks,
        # and all cancellation coins before the common potential quotient.
        tables, potential = exact_tables(), F(2, 3)
        weights = [F((i*3+1) % 7+1, 11) for i in range(16)]
        feature = feature_convolution(weights, 2)
        caps, _ = feature_caps(logs(tables), np.log2(float(potential)), feature, 1, 2, q=4)
        exact_cases = {(0, 0): [[(F(1), (0,)*4)]],
                       (1, 0): [[(F(1), tuple(int(i != missing) for i in range(4)))]
                                for missing in range(4)],
                       (0, 1): [[(F(2, 3), (1,)*4),
                                 (F(1, 3), tuple(int(i in edge) for i in range(4)))]
                                for edge in combinations(range(4), 2)]}
        for index, kind in enumerate(categories(1)):
            for fixed_geometry in exact_cases[kind]:
                values = [[F(0), F(0)] for _ in range(3)]
                for left, right in product(range(16), repeat=2):
                    e = int(left != 0)+int(right != 0)
                    for probability, sparse in fixed_geometry:
                        occupancy = tuple(sparse[i]+((left >> i) & 1)+((right >> i) & 1)
                                          for i in range(4))
                        vector = action(tables, occupancy, potential)
                        for state in range(2):
                            values[e][state] += weights[left]*weights[right]*probability*vector[state]
                for e in range(3):
                    exact = max(values[e][0], values[e][1]/potential)
                    self.assertGreaterEqual(2**caps[index, e]*(1+1e-12), float(exact))

    def test_positive_two_floor_tail_matches_integer_column_words(self):
        polynomial = [[1, 2, 0], [3, 1, 4], [2, 0, 1]]
        terms = [(x, y, coefficient) for y, row in enumerate(polynomial)
                 for x, coefficient in enumerate(row) if coefficient]
        exact = 0
        for word in product(terms, repeat=4):
            if sum(t[0] for t in word) >= 4 and sum(t[1] for t in word) >= 5:
                exact += prod(t[2] for t in word)
        with np.errstate(divide='ignore'):
            result, _ = positive_bivariate_tail(np.log2(polynomial), 4, 4, 5)
        self.assertAlmostEqual(result, np.log2(exact), places=12)
        self.assertNotEqual(exact, sum(t[2] for t in terms)**4)
        self.assertEqual(positive_bivariate_tail([[0.]], 2, 0, 1)[0], -np.inf)

    def test_frozen_scope_and_source_coloring(self):
        self.assertEqual(len(CONTROLS), 3)
        for control, (h, label) in CONTROLS.items():
            members, z, potential = frozen_type_and_witness(label)
            self.assertEqual(len(members), 36)
            self.assertTrue(all(mask.bit_count() == 5 for mask in members))
            self.assertTrue(0 < z < 1 and potential < 0)
            self.assertEqual(comb(17, h), 17 if h == 1 else 136)
            self.assertEqual((34*(17-h), 17*h, 32*(17-h)),
                             (544, 17, 512) if h == 1 else (510, 34, 480))


if __name__ == '__main__':
    unittest.main()

"""Tiny exact MacWilliams and positive-measure controls; no production bank."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import prod
import unittest

import numpy as np

import hd_lambda_dual_kernel as dual
from hd_pair_identity_kernel import mask_table_from_entries
from hd_pair_state_kernel import common_moment_atoms
from hd_quaternary_outer_relations import multiply_gf4


def logs(values):
    with np.errstate(divide='ignore'):
        return np.log2(np.asarray(values, float))


def orthogonal_matrix(rows):
    # I + 2J is a symmetric involution over GF4 when rows is even.
    return tuple(tuple(2 ^ (i == j) for j in range(rows)) for i in range(rows))


def words(matrix, source):
    result = []
    for values in product(range(4), repeat=len(source)):
        word = []
        for row in matrix:
            value = 0
            for j, x in zip(source, values):
                value ^= multiply_gf4(row[j], x)
            word.append(value)
        result.append(tuple(word))
    return result


def enumerator(word_list, atoms):
    return sum((prod(atoms[i] for i, x in enumerate(word) if x)
                for word in word_list), F(0))


def support_counts(word_list):
    return Counter(sum(bool(x) << i for i, x in enumerate(word)) for word in word_list)


class DualKernelTests(unittest.TestCase):
    def test_raw_macwilliams_exact_tiny_orthogonal_codes(self):
        for rows in (2, 4):
            matrix = orthogonal_matrix(rows)
            for i, j in product(range(rows), repeat=2):
                dot = 0
                for k in range(rows):
                    dot ^= multiply_gf4(matrix[k][i], matrix[k][j])
                self.assertEqual(dot, int(i == j))
            atoms = (F(0), F(1), F(1, 4), F(2, 3))[:rows]
            for arity in range(rows+1):
                for source in combinations(range(rows), arity):
                    complement = tuple(i for i in range(rows) if i not in source)
                    primal, dual_words = words(matrix, source), words(matrix, complement)
                    transformed = sum((prod(1-atoms[i] if x else 1+3*atoms[i]
                                            for i, x in enumerate(word))
                                       for word in dual_words), F(0))/4**len(complement)
                    self.assertEqual(enumerator(primal, atoms), transformed)

    def test_noncommuting_atom_products_and_minkowski(self):
        rows, source = 4, (0, 2)
        matrix = orthogonal_matrix(rows)
        primal = words(matrix, source)
        counts = support_counts(words(matrix, (1, 3)))
        masks, raw = tuple(counts), tuple(counts.values())
        atoms = np.array((0., .25, 1.))
        matrices = np.array(([[2., 1.], [0., 3.]], [[1., 0.], [2., 1.]],
                             [[1., 3.], [1., 0.]]))
        for groups in (1, 2, 3):
            entries = dual.entry_logs_from_atoms(logs(atoms), logs(matrices), groups, q=4)
            table = mask_table_from_entries(*entries, rows=rows)
            got = np.exp2(dual.column_log_bound(masks, raw, 2, groups, table, q=4))
            direct = np.zeros((2, 2))
            for choices in product(range(len(atoms)), repeat=rows):
                weight = np.eye(2)
                for choice in choices:
                    weight = weight @ matrices[choice]
                column = float(enumerator(primal, atoms[list(choices)]))
                direct += weight*column**groups
            if groups == 1:
                np.testing.assert_allclose(got, direct, rtol=2e-14)
            else:
                self.assertTrue(np.all(got >= direct*(1-2e-14)))
            for mask in range(1 << rows):
                expected = np.eye(2)
                for i in range(rows):
                    expected = expected @ np.exp2(entries[mask >> i & 1])
                np.testing.assert_allclose(np.exp2(table[mask]), expected, rtol=2e-14)

    def test_endpoint_atoms_and_mass_not_powered(self):
        rows, groups, b, q = 4, 3, 2, 4
        counts = support_counts(words(orthogonal_matrix(rows), (0, 1)))
        masks, raw = tuple(counts), tuple(counts.values())
        mass = np.array([[2., 3.], [1., 4.]])
        mass_word = np.linalg.matrix_power(mass, rows)
        for atom in (0., 1.):
            entries = dual.entry_logs_from_atoms(logs([atom]), logs([mass]), groups, q=q)
            table = mask_table_from_entries(*entries, rows=rows)
            got = np.exp2(dual.column_log_bound(masks, raw, b, groups, table, q=q))
            expected = mass_word*(q**((rows-b)*groups) if atom else 1)
            np.testing.assert_allclose(got, expected, rtol=1e-14)
            if atom:
                self.assertTrue(np.all(np.isneginf(entries[1])))
                self.assertTrue(np.all(np.isneginf(table[1:])))
            else:
                np.testing.assert_array_equal(entries[0], entries[1])

    def test_unmerged_categories_batch_and_zero_counts(self):
        entries = dual.entry_logs_from_atoms(logs([.5]), logs([[[2, 1], [0, 3]]]), 2, q=4)
        table = mask_table_from_entries(*entries, rows=2)
        masks = np.array([[0, 3, 3, 1], [0, 1, 1, 3]])
        counts = (1, 1, 2, 0)
        batch = dual.column_log_bound(masks, counts, 1, 2, table, q=4)
        for i in range(2):
            expected = dual.column_log_bound(masks[i, :2], (1, 3), 1, 2, table, q=4)
            np.testing.assert_allclose(batch[i], expected, rtol=0, atol=1e-14)
        np.testing.assert_array_equal(batch, dual.column_log_bound(
            masks, np.array([counts, counts], dtype=object), 1, 2, table, q=4))

    def test_common_measure_api_and_tiny_full_table(self):
        measure = dict(log2_probability=-2., log2_atom00=-3., log2_atom01=-2.,
                       log2_atom10=-1., log2_scale01=-4., log2_scale10=-5., log2_scale11=-6.)
        for groups in (1, 18):
            entries = dual.entry_log_matrices(measure, groups)
            expected = dual.entry_logs_from_atoms(*common_moment_atoms(measure), groups)
            for got, want in zip(entries, expected):
                np.testing.assert_array_equal(got, want)
            np.testing.assert_array_equal(dual.mask_log_matrices(measure, groups, rows=2),
                                          mask_table_from_entries(*entries, rows=2))

    def test_extreme_logs_and_zero_measure(self):
        atoms = np.array([-np.inf, -1e-12, 0.])
        matrices = logs([[[2, 1], [1, 3]]]*3)-2000
        low = dual.entry_logs_from_atoms(atoms, matrices, 18)
        ordinary = dual.entry_logs_from_atoms(atoms, matrices+2000, 18)
        for a, b in zip(low, ordinary):
            np.testing.assert_allclose(a, b-2000, rtol=0, atol=5e-13)
        table = mask_table_from_entries(*low, rows=2)
        self.assertTrue(np.all(np.isfinite(dual.column_log_bound([0, 3], [1, 255], 1, 18, table))))
        empty = dual.entry_logs_from_atoms([0.], np.full((1, 2, 2), -np.inf), 2)
        table = mask_table_from_entries(*empty, rows=2)
        self.assertTrue(np.all(np.isneginf(dual.column_log_bound([0, 3], [1, 255], 1, 2, table))))

    def test_guards_and_scope(self):
        matrices = np.zeros((1, 2, 2))
        for atom in ([.01], [np.nan], [np.inf]):
            with self.assertRaises(ValueError):
                dual.entry_logs_from_atoms(atom, matrices, 2)
        for groups, q in ((0, 4), (2, 3), (True, 4)):
            with self.assertRaises(ValueError):
                dual.entry_logs_from_atoms([-1.], matrices, groups, q=q)
        table = mask_table_from_entries(np.zeros((2, 2)), np.zeros((2, 2)), rows=2)
        for counts in ((1, 2), (1., 3.), (True, 3), (5, -1)):
            with self.assertRaises(ValueError):
                dual.column_log_bound([0, 3], counts, 1, 2, table, q=4)
        for masks in ([0, 4], [0., 3.], [0, -1]):
            with self.assertRaises(ValueError):
                dual.column_log_bound(masks, [1, 3], 1, 2, table, q=4)
        scope = dual.metadata()
        for key in ('atom_masses_powered', 'original_input_floor_included',
                    'rs_normalization_included', 'whole_code_certificate'):
            self.assertFalse(scope[key])


if __name__ == '__main__':
    unittest.main()

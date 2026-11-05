"""Tiny noncommuting identity-layer controls; no full mask table or bank."""
from itertools import product
import unittest

import numpy as np

import hd_pair_identity_kernel as identity
from hd_pair_state_kernel import common_moment_atoms


def logs(values):
    with np.errstate(divide='ignore'):
        return np.log2(np.asarray(values, float))


def ordered(values):
    result = np.eye(2, dtype=object)
    for value in values:
        result = result @ np.asarray(value, object)
    return np.asarray(result, float)


class IdentityKernelTests(unittest.TestCase):
    def test_all_tiny_masks_with_distinct_noncommuting_rows(self):
        for rows in (1, 2, 4):
            zero = [np.array([[1+i, 2], [1, 1]]) for i in range(rows)]
            nonzero = [np.array([[2, 1], [3+i, 2]]) for i in range(rows)]
            table = identity.mask_table_from_entries(logs(zero), logs(nonzero), rows=rows)
            for mask in range(1 << rows):
                direct = ordered([nonzero[i] if mask >> i & 1 else zero[i] for i in range(rows)])
                np.testing.assert_allclose(np.exp2(table[mask]), direct, rtol=6e-15)

    def test_common_entries_and_column_join(self):
        zero, nonzero = np.array([[1, 2], [0, 1]]), np.array([[2, 0], [3, 1]])
        table = identity.mask_table_from_entries(logs(zero), logs(nonzero), rows=4)
        first, second = zip(*product(range(16), repeat=2))
        got = identity.branch_log_matrices(first, second, table)
        for index, (a, b) in enumerate(zip(first, second)):
            sequence = [nonzero if mask >> i & 1 else zero for mask in (a, b) for i in range(4)]
            np.testing.assert_allclose(np.exp2(got[index]), ordered(sequence), rtol=8e-15)

    def test_five_atom_law_and_mass_not_powered(self):
        measure = dict(log2_probability=-2., log2_atom00=-3., log2_atom01=-2.,
                       log2_atom10=-1., log2_scale01=-4., log2_scale10=-5., log2_scale11=-6.)
        atoms, matrices = common_moment_atoms(measure)
        for groups in (1, 3, 18):
            zero, nonzero = identity.entry_log_matrices(measure, groups)
            np.testing.assert_allclose(np.exp2(zero), np.exp2(matrices).sum(axis=0))
            expected = (np.exp2(matrices)*np.exp2(groups*atoms)[:, None, None]).sum(axis=0)
            np.testing.assert_allclose(np.exp2(nonzero), expected, rtol=2e-15)
            table = identity.mask_log_matrices(measure, groups, rows=2)
            np.testing.assert_array_equal(table, identity.mask_table_from_entries(zero, nonzero, rows=2))

    def test_tiny_atoms_direct_full_chronology(self):
        measure = dict(log2_probability=-1., log2_atom00=-2., log2_atom01=-1.,
                       log2_atom10=-3., log2_scale01=-2., log2_scale10=-3., log2_scale11=-4.)
        atomlogs, matrixlogs = common_moment_atoms(measure)
        atoms, matrices = np.exp2(atomlogs), np.exp2(matrixlogs)
        table = identity.mask_log_matrices(measure, 2, rows=2)
        got = identity.branch_log_matrices([1], [2], table)[0]
        expected = np.zeros((2, 2))
        for choices in product(range(5), repeat=4):
            factor = atoms[choices[0]]**2*atoms[choices[3]]**2
            expected += factor*ordered(matrices[list(choices)])
        np.testing.assert_allclose(np.exp2(got), expected, rtol=5e-15)

    def test_log_underflow_and_exact_zeros(self):
        base = logs([[1, 2], [3, 1]])
        table = identity.mask_table_from_entries(base-1000, base-2000, rows=4)
        actual = identity.branch_log_matrices([5], [3], table)[0]
        expected = np.log2(ordered([np.exp2(base)]*8))-12000
        np.testing.assert_allclose(actual, expected, rtol=0, atol=4e-12)
        zero = identity.mask_table_from_entries(logs(np.eye(2)), np.full((2, 2), -np.inf), rows=2)
        self.assertTrue(np.all(np.isneginf(zero[1:])))

    def test_guards_and_scope(self):
        for rows in (0, 17):
            with self.assertRaises(ValueError):
                identity.mask_table_from_entries(np.zeros((2, 2)), np.zeros((2, 2)), rows=rows)
        table = identity.mask_table_from_entries(np.zeros((2, 2)), np.zeros((2, 2)), rows=2)
        for first in ([4], [-1], [1.5], [0]*8193):
            with self.assertRaises(ValueError):
                identity.branch_log_matrices(first, [0]*len(first), table)
        self.assertFalse(identity.metadata()['outer_crossed_pair_layer'])
        self.assertTrue(identity.metadata()['global_paired_route_unchanged'])
        self.assertFalse(identity.metadata()['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()

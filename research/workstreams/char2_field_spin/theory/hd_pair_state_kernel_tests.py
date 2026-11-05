"""Tiny exact/direct controls for the FLOAT two-column state kernel."""
from fractions import Fraction as F
from itertools import product
import unittest

import numpy as np

import hd_pair_state_kernel as kernel


def logs(values):
    values = np.asarray(values, dtype=float)
    with np.errstate(divide='ignore'):
        return np.log2(values)


def ordered(matrices, *, exact=False):
    answer = np.array([[1, 0], [0, 1]], dtype=object if exact else float)
    for matrix in matrices:
        answer = answer @ matrix
    return answer


def pair_value(t, u, kind, q=4, mode=kernel.EXACT_RATIO):
    if kind == 0:
        return F(1)
    if kind == 1:
        return t*u
    if mode == kernel.EXACT_RATIO:
        return (t+u+(q-3)*t*u)/(q-1)
    return (t+u+(q-4)*t*u)/(q-2)


def direct_measure(atoms, matrices, rows, delta, branches, groups, *, exact=False):
    """Independent enumeration of both natural atom/state chains.

    branches consists of (S,T,coefficient). The branch sum is taken BEFORE
    its groups power, unlike the separate-branch Minkowski majorant.
    """
    answer = np.zeros((2, 2), dtype=object if exact else float)
    for choices in product(range(len(atoms)), repeat=2*rows):
        first, second = choices[:rows], choices[rows:]
        scalar = F(0) if exact else 0.
        for left, right, coefficient in branches:
            value = coefficient
            for i in range(rows):
                kind = ((left >> i) & 1)+((right >> (i ^ delta)) & 1)
                value *= pair_value(atoms[first[i]], atoms[second[i ^ delta]], kind)
            scalar += value
        chain = ordered([matrices[i] for i in choices], exact=exact)
        answer += scalar**groups*chain
    return np.asarray(answer, dtype=float)


class PairStateKernelTests(unittest.TestCase):
    def setUp(self):
        self.atoms = [F(1, 3), F(3, 4)]
        self.matrices = np.array([
            [[F(1, 2), F(1, 5)], [F(1, 7), F(1, 3)]],
            [[F(1, 4), F(1, 3)], [F(1, 2), F(1, 9)]],
        ], dtype=object)

    def test_pair_functions_and_symmetric_domain(self):
        for t, u in product((F(0), F(1, 4), F(1)), repeat=2):
            for mode in (kernel.EXACT_RATIO, kernel.SYMMETRIC_MAJORANT):
                got = np.exp2(kernel.support_pair_logs(logs(t), logs(u), q=4, mode=mode))
                expected = [pair_value(t, u, k, mode=mode) for k in range(3)]
                np.testing.assert_allclose(got, np.asarray(expected, float), atol=1e-15)
                if mode == kernel.SYMMETRIC_MAJORANT:
                    # Equal nonzero amplitudes have both outputs nonzero.
                    self.assertGreaterEqual(float(expected[2]), float(t*u))
        with self.assertRaises(ValueError):
            kernel.support_pair_logs(.1, -1, mode=kernel.SYMMETRIC_MAJORANT)

    def test_leaf_fraction_enumeration(self):
        for groups in (1, 3):
            got = np.exp2(kernel.leaf_tensors(logs(self.atoms), logs(self.matrices), groups, q=4))
            expected = np.zeros((3, 2, 2, 2, 2), dtype=object)
            for kind, i, j, a, b, c, d in product(range(3), range(2), range(2),
                                                 range(2), range(2), range(2), range(2)):
                expected[kind, a, b, c, d] += (pair_value(self.atoms[i], self.atoms[j], kind)**groups
                    * self.matrices[i, a, b]*self.matrices[j, c, d])
            np.testing.assert_allclose(got, np.asarray(expected, float), rtol=2e-15, atol=1e-16)

    def test_masses_not_powered(self):
        mats = np.asarray(self.matrices, float)
        total = mats.sum(axis=0)
        expected = total[:, :, None, None]*total[None, None, :, :]
        for groups in (1, 18):
            got = np.exp2(kernel.leaf_tensors([0., 0.], logs(mats), groups))
            for kind in range(3):
                np.testing.assert_allclose(got[kind], expected, rtol=2e-14)

    def test_noncommuting_tree_all_xors(self):
        for rows in (1, 2, 4, 8):
            first = [np.array([[1+i, 2], [1, 2+i]], dtype=object) for i in range(rows)]
            second = [np.array([[2, 1+i], [3+i, 1]], dtype=object) for i in range(rows)]
            p, q = ordered(first, exact=True), ordered(second, exact=True)
            expected = np.asarray(p[:, :, None, None]*q[None, None, :, :], float)
            for delta in range(rows):
                leaves = [first[i][:, :, None, None]*second[i ^ delta][None, None, :, :]
                          for i in range(rows)]
                got = kernel.tensor_tree(logs(leaves), delta)
                np.testing.assert_allclose(np.exp2(got), expected, rtol=1e-14)
                np.testing.assert_allclose(np.exp2(kernel.join_columns(got)),
                                           np.asarray(p @ q, float), rtol=1e-14)

    def test_two_rows_fraction_atom_enumeration(self):
        for groups in (1, 3):
            leaves = kernel.leaf_tensors(logs(self.atoms), logs(self.matrices), groups, q=4)
            for delta in range(2):
                for left, right in ((0, 0), (1, 2), (3, 1), (3, 3)):
                    direct = direct_measure(self.atoms, self.matrices, 2, delta,
                                            [(left, right, F(1))], groups, exact=True)
                    got = kernel.branch_log_matrix(left, right, delta, leaves, rows=2)
                    np.testing.assert_allclose(np.exp2(got), direct, rtol=8e-15, atol=1e-18)

    def test_four_rows_direct_all_xors(self):
        atoms, matrices = list(map(float, self.atoms)), np.asarray(self.matrices, float)
        leaves = kernel.leaf_tensors(logs(atoms), logs(matrices), 2, q=4)
        for delta in range(4):
            direct = direct_measure(atoms, matrices, 4, delta, [(5, 11, 1.)], 2)
            got = kernel.branch_log_matrix(5, 11, delta, leaves, rows=4)
            np.testing.assert_allclose(np.exp2(got), direct, rtol=1e-14, atol=1e-19)

    def test_branch_minkowski_counts_stay_outside_roots(self):
        branches = [(1, 2, F(2)), (3, 1, F(3)), (0, 0, F(1, 4))]
        for groups in (1, 3):
            leaves = kernel.leaf_tensors(logs(self.atoms), logs(self.matrices), groups, q=4)
            matrices = kernel.branch_log_matrices([b[0] for b in branches], [b[1] for b in branches],
                                                 1, leaves, rows=2)
            terms = matrices/groups+logs([b[2] for b in branches])[:, None, None]
            bound = np.exp2(groups*np.logaddexp2.reduce(terms, axis=0))
            exact = direct_measure(self.atoms, self.matrices, 2, 1, branches, groups, exact=True)
            if groups == 1:
                np.testing.assert_allclose(bound, exact, rtol=8e-15)
            else:
                self.assertTrue(np.all(bound >= exact*(1-1e-14)))
                one = branches[:1]
                exact_one = direct_measure(self.atoms, self.matrices, 2, 1, one, groups, exact=True)
                np.testing.assert_allclose(np.exp2(matrices[0])*float(one[0][2])**groups,
                                           exact_one, rtol=8e-15)

    def test_common_five_atom_measure(self):
        measure = dict(log2_probability=-2., log2_atom00=-3., log2_atom01=-2.,
                       log2_atom10=-1., log2_scale01=-4., log2_scale10=-5., log2_scale11=-6.)
        marks, mats = kernel.common_moment_atoms(measure)
        p = .25
        for m in (0, 1, 3):
            factors = np.ones(5) if m == 0 else np.exp2(m*marks)
            got = (np.exp2(mats)*factors[:, None, None]).sum(axis=0)
            expected = np.array([[(1-p if m == 0 else 0)+p*2.**(-3*m), 2.**(-4-2*m)],
                                 [2.**(-5-m), 2.**-6]])
            np.testing.assert_allclose(got, expected)
        leaves = kernel.leaves_from_measure(measure, 2, q=4)
        direct = direct_measure(np.exp2(marks), np.exp2(mats), 2, 1, [(1, 2, 1.)], 2)
        np.testing.assert_allclose(np.exp2(kernel.branch_log_matrix(1, 2, 1, leaves, rows=2)),
                                   direct, rtol=1e-14, atol=1e-19)

    def test_batch_matches_single(self):
        leaves = kernel.leaf_tensors(logs(self.atoms), logs(self.matrices), 2, q=4)
        first, second = [0, 1, 5, 15], [8, 3, 11, 15]
        for delta in range(4):
            got = kernel.branch_log_matrices(first, second, delta, leaves, rows=4)
            expected = [kernel.branch_log_matrix(s, t, delta, leaves, rows=4)
                        for s, t in zip(first, second)]
            np.testing.assert_array_equal(got, expected)

    def test_log_underflow_and_exact_zeros(self):
        base = kernel.leaf_tensors(logs(self.atoms), logs(self.matrices), 18, q=4)
        tiny = kernel.leaf_tensors(logs(self.atoms), logs(self.matrices)-1000, 18, q=4)
        left = kernel.branch_log_matrix(5, 11, 3, base, rows=4)
        right = kernel.branch_log_matrix(5, 11, 3, tiny, rows=4)
        self.assertTrue(np.all(np.isfinite(right)))
        np.testing.assert_allclose(right, left-8000, rtol=0, atol=2e-12)
        zero = kernel.leaf_tensors([-np.inf], logs([np.eye(2)]), 18)
        self.assertTrue(np.all(np.isneginf(zero[1:])))
        result = kernel.branch_log_matrix(0, 0, 3, zero, rows=4)
        np.testing.assert_array_equal(result, logs(np.eye(2)))

    def test_validation_and_scope(self):
        leaves = kernel.leaf_tensors([0.], logs([np.eye(2)]), 1)
        for left, right, delta, rows in (([16], [0], 0, 4), ([1.5], [0], 0, 4),
                                         ([0], [0], 4, 4), ([0], [0], 0, 3),
                                         ([0]*8193, [0]*8193, 0, 4)):
            with self.assertRaises(ValueError):
                kernel.branch_log_matrices(left, right, delta, leaves, rows=rows)
        with self.assertRaises(ValueError):
            kernel.leaf_tensors([np.nan], logs([np.eye(2)]), 1)
        info = kernel.metadata()
        self.assertFalse(info['whole_code_certificate'])
        self.assertFalse(info['atom_matrix_masses_powered'])
        self.assertTrue(info['both_chronological_state_chains_retained'])


if __name__ == '__main__':
    unittest.main()

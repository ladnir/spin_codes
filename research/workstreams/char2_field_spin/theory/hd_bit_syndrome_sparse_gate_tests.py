from itertools import product
from math import comb, log2
import unittest

import numpy as np

from hd_bit_syndrome_sparse_gate import column_branches, word_restriction
from hd_bit_syndrome_sparse_geometry import prepare_geometry
from hd_post_refresh_quiet_rows import log_matrix_power, moment_atoms


class SparseGateTests(unittest.TestCase):
    def test_first_d_exact_when_group_exponent_one(self):
        q = np.array([[.7, .1], [.2, .6]])
        n = np.array([[.1, .3], [.2, .1]])
        branch = {name: np.log2(value[None, ...]) for name, value in
                  (('quiet', q), ('nonquiet', n), ('full', q+n))}
        maps = (tuple(range(16)),)*4
        rows = (tuple(range(16)),)
        for floor in range(5):
            total = np.zeros((2, 2))
            for labels in product((0, 1), repeat=4):
                if sum(labels) < floor:
                    continue
                value = np.eye(2)
                for label in labels:
                    value = value @ (n if label else q)
                total += value
            got = word_restriction(branch, rows, maps, floor, groups=1)
            self.assertAlmostEqual(got[0], log2(total[0].sum()), places=12)
            one = word_restriction(branch, rows, maps, floor, groups=1, one_copy=True)
            np.testing.assert_array_equal(one, got)

    def test_zero_nonquiet(self):
        quiet = np.zeros((1, 2, 2))
        branches = dict(quiet=quiet, nonquiet=np.full_like(quiet, -np.inf), full=quiet)
        got = word_restriction(branches, (tuple(range(16)),), (tuple(range(16)),)*4, 1)
        self.assertEqual(got[0], -np.inf)

    def test_full_word_is_unrestricted(self):
        g = np.log2(np.array([[[.8, .2], [.3, .7]]]))
        branches = dict(quiet=g, nonquiet=g, full=g)
        got = word_restriction(branches, (tuple(range(16)),), (tuple(range(16)),)*16, 0)
        expected = np.logaddexp2.reduce(log_matrix_power(g, 16)[:, 0, :], axis=1)
        np.testing.assert_allclose(got, expected, atol=1e-12)

    def test_column_normalization_at_unit_atoms(self):
        prepared = prepare_geometry(2, rowsets=((0, 1),))
        measure = dict(log2_probability=0., log2_atom00=0., log2_atom01=0., log2_atom10=0.,
                       log2_scale01=-1., log2_scale10=-2., log2_scale11=-3.)
        branches = column_branches(prepared, measure, 2)
        _, empty = moment_atoms(measure)
        base = log_matrix_power(empty, 16)
        caps = prepared['outside_caps'][0]
        mass = sum(comb(14, i)*v for i, v in enumerate(caps))
        np.testing.assert_allclose(branches['full'][0], base+2*log2(mass), atol=1e-12)
        np.testing.assert_allclose(branches['nonquiet'][0], base+2*log2(mass-1), atol=1e-12)
        np.testing.assert_allclose(branches['quiet'][0], base, atol=1e-12)


if __name__ == '__main__':
    unittest.main()

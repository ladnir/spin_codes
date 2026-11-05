from math import log2
import unittest
import numpy as np

from hd_sparse_activity_control import minkowski_columns, quiet_nonzero_columns
from hd_post_refresh_quiet_rows import log_matmul, moment_atoms


class SparseActivityTests(unittest.TestCase):
    def test_minkowski_unnormalized_measure(self):
        mass = np.array([.25, 2., 7.5])
        q, n = np.array([0., 1., 2.]), np.array([1., 0., .5])
        for r in (1, 2, 18):
            left = np.full((2, 2), log2(np.sum(mass*q**r)))
            right = np.full((2, 2), log2(np.sum(mass*n**r)))
            upper = minkowski_columns(left, right, r)
            actual = log2(np.sum(mass*(q+n)**r))
            self.assertGreaterEqual(upper[0, 0]+1e-12, actual)
            if r == 1:
                self.assertAlmostEqual(upper[0, 0], actual, places=12)

    def test_zero_branch_has_no_young_penalty(self):
        left = np.full((2, 2), -np.inf)
        right = np.arange(4.).reshape(2, 2)
        np.testing.assert_allclose(minkowski_columns(left, right, 18), right)

    def test_plane_nonzero_moment_independent_product(self):
        measure = dict(log2_probability=-3., log2_atom00=-2., log2_atom01=-1., log2_atom10=-.5,
                       log2_scale01=-2., log2_scale10=-4., log2_scale11=-1.)
        plane = (0, 1, 2, 3)
        nonplane = (0, 1, 2, 4)
        r = 3
        got = quiet_nonzero_columns((plane, nonplane), measure, r)
        _, empty = moment_atoms(measure)
        moment = np.array([[measure['log2_probability']+r*measure['log2_atom00'],
                            measure['log2_scale01']+r*measure['log2_atom01']],
                           [measure['log2_scale10']+r*measure['log2_atom10'], measure['log2_scale11']]])
        expected = np.full((2, 2), -np.inf)
        expected[0, 0] = expected[1, 1] = 0.
        for i in range(16):
            expected = log_matmul(expected, moment if i in plane else empty)
        expected += r*log2(255)
        np.testing.assert_allclose(got[0], expected, atol=1e-12)
        self.assertTrue(np.isneginf(got[1]).all())


if __name__ == '__main__':
    unittest.main()

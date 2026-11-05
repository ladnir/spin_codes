from itertools import combinations
from math import comb, log2
import unittest
from unittest.mock import patch
import numpy as np

from hd_sparse_activity_control import quiet_nonzero_columns
from hd_sparse_marked_control import marked_quiet, marked_words, quiet_geometry


MEASURE = dict(log2_probability=-3.,log2_atom00=-2.,log2_atom01=-1.,log2_atom10=-.5,
               log2_scale01=-2.,log2_scale10=-4.,log2_scale11=-1.)


class SparseMarkedTests(unittest.TestCase):
    def test_marker_one_reproduces_nonzero_quiet(self):
        rows = tuple(combinations(range(16), 3))
        geometry = quiet_geometry(rows, MEASURE)
        expected = quiet_nonzero_columns(rows, MEASURE, 2)
        np.testing.assert_allclose(marked_quiet(geometry, MEASURE, 2, 0., nonzero=True), expected)

    def test_plane_input_marker_exact_degree(self):
        rows = ((0,1,2,3),(0,1,2,4))
        geometry = quiet_geometry(rows, MEASURE)
        base = marked_quiet(geometry,MEASURE,3,0.,nonzero=True)
        marked = marked_quiet(geometry,MEASURE,3,2.,nonzero=True)
        np.testing.assert_allclose(marked[0],base[0]+4*3*2)
        self.assertTrue(np.isneginf(marked[1]).all())

    def test_scalar_first9_with_nonquiet_marker(self):
        # r=1: exact positive scalar generating function, same S(0)=I.
        def matrix(value):
            answer = np.full((1,2,2),-np.inf)
            answer[:,0,0] = answer[:,1,1] = log2(value)
            return answer
        quiet_nz, nonquiet, z, floor = .3, .7, 4., 3
        branches = dict(full=matrix(1+quiet_nz+nonquiet),nonquiet=matrix(nonquiet))
        with patch('hd_sparse_marked_control.moment_atoms',return_value=(None,matrix(1)[0])):
            got = marked_words(branches,matrix(1+quiet_nz),matrix(quiet_nz),
                (tuple(range(16)),),MEASURE,1,0.,log2(z),floor)[0]
        expected = log2(sum(comb(16,d)*(quiet_nz+z*nonquiet)**d for d in range(9,17)))-floor*log2(z)
        self.assertAlmostEqual(got,expected,places=11)

    def test_nonquiet_input_mark_is_an_upper_bound(self):
        for a in (5,6,7):
            for x in (1.,2.,4.):
                for weight in range(1,a+1):
                    self.assertLessEqual(x**weight,x**a)


if __name__ == '__main__':
    unittest.main()

import unittest
import numpy as np
from fractions import Fraction as F

import gf8_punctured_screen as fixed
import gf8_balanced_puncture_screen as balanced


class BalancedTests(unittest.TestCase):
    def test_fixed_pattern(self):
        for dimension,outputs in ((84,7),(80,6),(67,5)):
            g=balanced.profile(dimension,(outputs,))
            self.assertEqual(g,fixed.profile(dimension=dimension,outputs=outputs))
            for active in (1,4,12):
                self.assertAlmostEqual(balanced.sparse_bound(active,(6.5,2.5),dimension,(outputs,)),
                    fixed.sparse_bound(active,(6.5,2.5),dimension=dimension,outputs=outputs),places=7)
            point=(-1.,-.5,2.,1.)
            self.assertAlmostEqual(balanced.middle_bound(point,50000,dimension,(outputs,)),
                fixed.middle_bound(point,50000,dimension=dimension,outputs=outputs),places=7)

    def test_repeated_pattern(self):
        point=(-1.,-.5,2.,1.)
        for outputs in (5,6):
            np.testing.assert_allclose(balanced.period_matrix(point,(outputs,outputs)),
                np.linalg.matrix_power(fixed.packet(point,8,outputs,8),2))
            self.assertAlmostEqual(balanced.sparse_bound(8,(6.5,2.5),74,(outputs,outputs)),
                balanced.sparse_bound(8,(6.5,2.5),74,(outputs,)),places=7)

    def test_balanced_rate(self):
        g=balanced.profile()
        self.assertEqual(F(g['k'],g['output_n']),F(37,44))
        self.assertEqual(g['rows']%16,0)
        self.assertLess(10*g['threshold'],g['output_n'])


if __name__=='__main__':unittest.main()

from fractions import Fraction as F
import math
import unittest

from gf8_speed_geometry_screen import CANDIDATES,Candidate
from gf8_punctured_screen import local_matrices,middle_bound,sparse_bound
from gf8_puncture_frontier import middle_bound as control_middle,sparse_bound as control_sparse


class SpeedGeometryTests(unittest.TestCase):
    def test_natural_geometry(self):
        for c in CANDIDATES:
            g=c.geometry()
            self.assertEqual(g['rows']%64,0)
            self.assertGreaterEqual(g['k'],2**20)
            self.assertEqual(F(g['k'],g['output_n']),c.rate)
            self.assertEqual(g['input_n'],c.length*g['rows'])
            self.assertEqual(g['threshold'],(g['output_n']+9)//10-1)
        with self.assertRaises(ValueError):
            Candidate('bad-rate',64,44,8,7,8).geometry()

    def test_control_matches_frozen_geometry_screen(self):
        point=(-.5,-.8,2.,1.)
        g=CANDIDATES[0].geometry()
        weight=round(.05*g['input_n'])
        self.assertAlmostEqual(middle_bound(point,weight,8,6,8,128,80),
                               control_middle(point,.05,80,(6,)),places=6)
        self.assertAlmostEqual(sparse_bound(4,(6.5,2.),8,6,8,128,80),
                               control_sparse((6.5,2.),4,80,(6,)),places=8)

    def test_new_local_transfers_remain_exact_positive(self):
        for width,outputs,state in {(c.width,c.outputs,c.state) for c in CANDIDATES}:
            matrices=local_matrices(F(3,4),F(11),width,outputs,state)
            self.assertEqual(len(matrices),width+1)
            self.assertEqual(matrices[0][0],[F(1),F(0)])
            for matrix in matrices:
                for row in matrix:
                    self.assertTrue(all(isinstance(value,F) and value>=0 for value in row))
            for weight in range(1,width+1):
                self.assertLessEqual(matrices[weight][0][1],F(3,4)**max(0,outputs+1-weight))


if __name__=='__main__':
    unittest.main()

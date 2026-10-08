from fractions import Fraction as F
import math
import unittest

from gf8_feedthrough_mds_screen import weighted_matrices
from gf8_mds_geometry_screen import local_matrices, packet


class MdsGeometryTests(unittest.TestCase):
    def test_agreement_with_independent_width_eight_transfer(self):
        for state in (4,8):
            for y,tau in ((F(9,10),F(6)),(F(1,2),F(255)),(F(1),F(1))):
                got = local_matrices(y,tau,8,state)
                ref = weighted_matrices(y,tau,False,False,state,False)
                for w in range(9):
                    for i in range(2):
                        for j in range(2):
                            scaled = got[w][i][j]*math.comb(8,w)
                            if (i,j)==(0,1):
                                self.assertLessEqual(scaled,ref[w][i][j])
                            else:
                                self.assertEqual(scaled,ref[w][i][j])

    def test_full_output_zero_exclusions(self):
        for width,state in ((8,16),(16,8),(16,16),(32,16)):
            values = local_matrices(F(1,10),F(6),width,state)
            self.assertEqual(values[0][0][0],1)
            self.assertEqual(values[0][0][1],0)
            self.assertEqual(values[0][1][0],0)
            self.assertTrue(all(x>=0 for matrix in values for row in matrix for x in row))
            # A fully live input of weight <=state cannot cancel into state0
            # when the incoming state is zero.
            self.assertTrue(all(values[w][0][0]==0 for w in range(1,min(width,state)+1)))


if __name__=='__main__':
    unittest.main()

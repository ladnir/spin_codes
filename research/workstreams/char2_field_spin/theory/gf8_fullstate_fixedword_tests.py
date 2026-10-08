from fractions import Fraction as F
import unittest
import gf8_fullstate_fixedword as target
from gf8_fullstate_screen import local_matrix as coarse


class FixedWordTests(unittest.TestCase):
    def test_birth_and_return(self):
        y=F(9,10)
        for state in (8,16):
            self.assertEqual(target.local_matrix(0,y,state)[0],[1,0])
            for w in range(1,9):
                matrix=target.local_matrix(w,y,state)
                self.assertEqual(matrix[0],[0,y**w])
                self.assertLessEqual(matrix[1][0],F(255)**(-state))

    def test_fixed_word_bound_dominates_no_rank_expansion(self):
        for state in (8,16):
            for w in range(9):
                y=F(3,4)
                direct=target.local_matrix(w,y,state)
                rank=coarse(w,y,F(255),state)
                for i in range(2):
                    for j in range(2):
                        self.assertLessEqual(direct[i][j],rank[i][j])


if __name__=='__main__':
    unittest.main()

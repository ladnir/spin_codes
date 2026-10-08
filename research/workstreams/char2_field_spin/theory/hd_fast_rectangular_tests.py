from fractions import Fraction as F
from itertools import combinations
import unittest

from gf8_packet_maps import MUL,INV,rank
from high_distance_fields import local_matrices as square
from hd_fast_rectangular import local_matrices,bound


class RectangularTests(unittest.TestCase):
    def test_existing_local_formula_agreement(self):
        for width,outputs,state in ((4,2,4),(4,4,8),(8,6,8)):
            self.assertEqual(local_matrices(F(3,4),F(17),width,outputs,state),
                             square(F(3,4),F(17),width,outputs,state,16))

    def test_tall_cauchy_selected_output_and_all_state_ranks(self):
        # Five raw rows, four columns; normalize the two bottom state rows.
        w,s,o=2,2,3
        raw=[[INV[i^(128+j)] for j in range(w+s)] for i in range(o+s)]
        work=[row[:] for row in raw]
        for j in range(s):
            pivot=o+j
            scale=INV[work[pivot][w+j]]
            work[pivot]=[MUL[scale][x] for x in work[pivot]]
            for i in range(o,o+s):
                if i!=pivot:
                    factor=work[i][w+j]
                    work[i]=[a^MUL[factor][b] for a,b in zip(work[i],work[pivot])]
        for count in range(o+1):
            for outputs in combinations(range(o),count):
                for include_state in (False,True):
                    rows=outputs+tuple(range(o,o+s)) if include_state else outputs
                    for size in range(w+s+1):
                        for columns in combinations(range(w+s),size):
                            matrix=tuple(tuple(work[i][j] for j in columns) for i in rows)
                            self.assertEqual(rank(matrix),min(len(rows),len(columns)))

    def test_saved_rectangular_witnesses(self):
        self.assertAlmostEqual(-bound((6.422157028430665,2.28363327218258),16),
                               -1380.0058453103907,places=5)
        self.assertAlmostEqual(-bound((-2.2528233247350977,-.10458061152706508,
                                      2.3462809389297363,1.8090850960628637),59474,False),
                               -80572.86627308716,places=5)


if __name__=='__main__':unittest.main()

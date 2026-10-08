from fractions import Fraction as F
import math
import unittest

from gf8_mds_geometry_screen import local_matrices as full_matrices
from gf8_punctured_screen import local_matrices,profile
from gf8_mds_packet import fullstate_matrices
from gf8_packet_maps import MUL,rank,matmul


class PuncturedPacketTests(unittest.TestCase):
    def test_unpunctured_agrees_exactly(self):
        for width,state in ((8,4),(8,8),(8,16),(16,8)):
            for y,tau in ((F(9,10),F(6)),(F(1,2),F(255)),(F(1),F(1))):
                self.assertEqual(local_matrices(y,tau,width,width,state),
                                 full_matrices(y,tau,width,state))

    def test_rate_and_all_output_zero_branch(self):
        for outputs,dimension in ((8,96),(7,84),(6,72),(4,48)):
            g=profile(outputs=outputs)
            self.assertEqual(g['dimension'],dimension)
            self.assertEqual(F(g['k'],g['output_n']),F(3,4))
            self.assertEqual(g['input_n']//8*outputs,g['output_n'])
            # Puncturing can create an all-zero emitted packet even when
            # the state returns to zero. The bound must retain this event.
            if outputs<8:
                m=local_matrices(F(1,10**12),F(6),8,outputs,8)
                self.assertGreater(m[8][1][0],F(1,6**30))

    def test_restricted_feedthrough_and_global_rank(self):
        f,a,c,eta=fullstate_matrices()
        for outputs in (4,6,7,8):
            # Two chronological packets, zero initial state. This checks
            # the full-row-rank property used by the preimage-count gate.
            acf=matmul(a,c)
            top=[tuple(row)+(0,)*8 for row in f[:outputs]]
            bottom=[tuple(MUL[eta][v] for v in acf[i])+tuple(f[i]) for i in range(outputs)]
            self.assertEqual(rank(top+bottom),2*outputs)
            self.assertEqual(rank(f[:outputs]),outputs)


if __name__=='__main__':
    unittest.main()

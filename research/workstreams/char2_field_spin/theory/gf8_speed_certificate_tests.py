from copy import deepcopy
from fractions import Fraction as F
from random import Random
import unittest

from gf8_packet_maps import MUL,matvec,xor_sum
from gf8_speed_maps import packet,F as FEED,A_EFFECTIVE,C,FACTORED_LAMBDA
from gf8_speed_profile import metadata
from gf8_speed_certificate import combine,validate_maps
from gf8_speed_middle_certificate import certify as middle,cell_upper,prefactor,WITNESSES
from gf8_speed_dense_certificate import certify as dense
from return_cover_certificate import interval


class SpeedCertificateTests(unittest.TestCase):
    def test_exact_map_and_factored_adjoint(self):
        checked=validate_maps()
        self.assertEqual(checked['selected_rank_cases'],8192)
        self.assertEqual(checked['factored_lambda'],107)
        rng=Random(7040)
        def dot(a,b):
            return xor_sum(MUL[x][y] for x,y in zip(a,b))
        for _ in range(100):
            x,r,y,h=([rng.randrange(256) for _ in range(4)] for j in range(4))
            dx,ds=([rng.randrange(1,256) for _ in range(4)] for j in range(2))
            out,following=packet(x,r,dx,ds)
            back,previous=packet(y,h,dx,ds,True)
            self.assertEqual(dot(out,y)^dot(following,h),dot(x,back)^dot(r,previous))
            # Both state maps are symmetric XOR-circulants. Eliminate F
            # independently through F=lambda*A_effective*C.
            ay=matvec(A_EFFECTIVE,y)
            moment=tuple(a^MUL[FACTORED_LAMBDA][b] for a,b in zip(h,ay))
            self.assertEqual(back,tuple(MUL[a][b] for a,b in zip(matvec(C,moment),dx)))
            self.assertEqual(previous,tuple(MUL[a^b][d] for a,b,d in zip(h,ay,ds)))

    def test_cell_upper_dominates_exact_weights(self):
        g=metadata();w=WITNESSES[0][1];factor=prefactor(w,g)
        upper=cell_upper(50381,50384,g,w,factor)
        points=sum((cell_upper(v,v,g,w,factor) for v in range(50381,50385)),interval(F(0)))
        self.assertLessEqual(points.hi,upper.hi)

    def test_complete_middle_and_dense(self):
        a,b=middle(),dense()
        self.assertEqual(a['weights'],[50381,1343487])
        self.assertEqual(b['weights'],[1343488,1679360])
        self.assertGreaterEqual(a['margin_bits_lower'],40)
        self.assertGreaterEqual(b['margin_bits_lower'],40)
        self.assertEqual(b['completion_factor_log2'],3358720)

    def receipts(self):
        g=metadata()
        a=dict(**g,complete=True,rigorous=True,margin_bits_lower=70,active_range=[1,2015],
               failure_upper_dyadic_hex=['0x1',hex(1<<70)])
        b=dict(**g,complete=True,rigorous=True,margin_bits_lower=80,weights=[50381,1343487],
               failure_upper=str(F(1,1<<80)))
        c=dict(**g,complete=True,rigorous=True,margin_bits_lower=90,weights=[1343488,1679360],
               failure_upper_dyadic_hex=['0x1',hex(1<<90)],omitted_symbols=419840,
               completion_factor_log2=3358720,
               dense_inner_hypothesis='invertible full causal inner; independent nonzero input diagonals')
        return a,b,c

    def test_aggregate_adds_failure_and_covers_all_sources(self):
        result=combine(*self.receipts())
        self.assertEqual(result['margin_bits_lower'],69)
        self.assertEqual(result['first_remaining_outer_weight'],50400)
        self.assertEqual(result['minimum_distance'],125952)
        self.assertTrue(result['whole_code'])

    def test_rejects_gaps_and_wrong_metadata(self):
        for index,key,value in ((0,'active_range',[1,2014]),(1,'weights',[50381,1343486]),
                               (2,'completion_factor_log2',0),(2,'outputs',4),
                               (0,'matrix_sha256','wrong'),(1,'rigorous',False)):
            receipts=deepcopy(self.receipts());receipts[index][key]=value
            with self.assertRaises((ValueError,ArithmeticError)):
                combine(*receipts)


if __name__=='__main__':
    unittest.main()

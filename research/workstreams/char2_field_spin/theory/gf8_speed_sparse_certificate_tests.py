from fractions import Fraction as F
import math
import unittest

from flint import ctx
from gf8_punctured_screen import local_matrices
from gf8_speed_sparse_certificate import ANCHORS,certify,certify_zone,polynomial_power,witnesses


class SpeedSparseCertificateTests(unittest.TestCase):
    def test_exact_two_packet_polynomial(self):
        local=local_matrices(F(3,4),F(9),4,3,4)
        previous=ctx.prec
        try:
            ctx.prec=512
            result=polynomial_power(local,2,8)
            for h in range(9):
                for i,j in ((0,0),(0,1),(1,0),(1,1)):
                    expected=sum((F(math.comb(4,a)*math.comb(4,h-a))*local[a][i][k]*local[h-a][k][j]
                        for a in range(max(0,h-4),min(4,h)+1) for k in range(2)),F(0))
                    value=result[2*i+j][h]
                    lo,hi=value.lower().fmpq(),value.upper().fmpq()
                    self.assertLessEqual(F(int(lo.p),int(lo.q)),expected)
                    self.assertGreaterEqual(F(int(hi.p),int(hi.q)),expected)
        finally:
            ctx.prec=previous

    def test_low_range_metadata(self):
        receipt=certify(last=16,precision=1024)
        self.assertTrue(receipt['complete'] and receipt['rigorous'])
        self.assertFalse(receipt['whole_code'])
        self.assertEqual(receipt['outer'],[64,40])
        self.assertEqual(receipt['nominal_rate'],'5/6')
        self.assertEqual(receipt['rows'],26240)
        self.assertEqual(receipt['input_n'],1679360)
        self.assertEqual(receipt['output_n'],1259520)
        self.assertEqual(receipt['output_threshold'],125951)
        self.assertGreaterEqual(receipt['margin_bits_lower'],60)

    def test_shallow_sparse_dip(self):
        receipt=certify_zone(64,64,64,ANCHORS[64],precision=1024)
        self.assertGreaterEqual(receipt['certified_integer_bits'],67)
        p,q=(int(x,16) for x in receipt['failure_upper_dyadic_hex'])
        self.assertLessEqual(F(p,q),F(1,2**67))

    def test_dense_witnesses_and_invalid_geometry(self):
        anchors=witnesses(2015)
        self.assertTrue(all(a in anchors for a in range(32,193)))
        with self.assertRaises(ValueError):
            polynomial_power([[[F(0)]*2 for _ in range(2)]]*9,1,1)
        with self.assertRaises(ValueError):
            certify(last=0)
        with self.assertRaises(ValueError):
            certify_zone(2,1,1,ANCHORS[1])


if __name__=='__main__':
    unittest.main()

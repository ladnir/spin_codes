from fractions import Fraction as F
import unittest
from gf8_rate_climb_sparse import certify,certify_zone,ANCHORS88
from gf8_punctured_sparse_certificate import certify_zone as frozen_zone,ANCHORS


class RateClimbSparseTests(unittest.TestCase):
    def test_frozen84_point_agreement(self):
        frozen=frozen_zone(1,1,1,ANCHORS[1],precision=512)
        generic=certify_zone(1,1,1,ANCHORS[1],dimension=84,precision=512)
        self.assertEqual(frozen,generic)

    def test_directed_low_range88(self):
        receipt=certify(last=16,precision=1024)
        self.assertTrue(receipt['complete'] and receipt['rigorous'])
        self.assertFalse(receipt['whole_code'])
        self.assertEqual(receipt['outer'],[128,88])
        self.assertEqual(receipt['nominal_rate'],'11/14')
        self.assertEqual(receipt['output_n'],1340416)
        self.assertGreaterEqual(receipt['margin_bits_lower'],75)

    def test_directed_dip88(self):
        receipt=certify_zone(224,224,224,ANCHORS88[224],precision=1024)
        self.assertGreaterEqual(receipt['certified_integer_bits'],134)
        p,q=(int(x,16) for x in receipt['failure_upper_dyadic_hex'])
        self.assertLessEqual(F(p,q),F(1,2**134))

    def test_keep6_low_range(self):
        receipt=certify(last=16,dimension=80,outputs=6,precision=1024)
        self.assertTrue(receipt['complete'] and receipt['rigorous'])
        self.assertEqual(receipt['nominal_rate'],'5/6')
        self.assertEqual(receipt['input_n'],1679360)
        self.assertEqual(receipt['output_n'],1259520)
        self.assertEqual(receipt['output_threshold'],125951)
        self.assertGreaterEqual(receipt['margin_bits_lower'],264)

    def test_keep5_low_range(self):
        receipt=certify(last=16,dimension=67,outputs=5,precision=1024)
        self.assertTrue(receipt['complete'] and receipt['rigorous'])
        self.assertEqual(receipt['nominal_rate'],'67/80')
        self.assertEqual(receipt['input_n'],2007040)
        self.assertEqual(receipt['output_n'],1254400)
        self.assertEqual(receipt['output_threshold'],125439)
        self.assertGreaterEqual(receipt['margin_bits_lower'],368)


if __name__=='__main__':
    unittest.main()

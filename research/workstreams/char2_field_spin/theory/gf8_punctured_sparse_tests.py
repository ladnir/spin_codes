from fractions import Fraction as F
import unittest

from gf8_punctured_sparse_certificate import metadata,certify_zone,certify,ANCHORS
from gf8_punctured_screen import local_matrices


class PuncturedSparseTests(unittest.TestCase):
    def test_shared_geometry(self):
        g=metadata()
        self.assertEqual(g['rows'],12544)
        self.assertEqual(g['input_n'],1605632)
        self.assertEqual(g['output_n'],1404928)
        self.assertEqual(g['output_threshold'],140492)
        self.assertEqual(g['outer'],[128,84])
        self.assertEqual(g['padded_k']*4,g['output_n']*3)

    def test_exact_local_coefficients(self):
        matrices=local_matrices(F(9,10),F(6),8,7,8)
        for matrix in matrices:
            for row in matrix:
                self.assertTrue(all(isinstance(v,(int,F)) and v>=0 for v in row))
        self.assertEqual(matrices[0][0],[1,0])
        self.assertTrue(all(matrices[w][0][0]==0 for w in range(1,9)))
        # Eight inputs may map to zero in the seven retained outputs.
        self.assertGreater(matrices[8][0][1],0)

    def test_directed_singleton_and_short_cover(self):
        point=certify_zone(1,1,1,ANCHORS[1],precision=512)
        self.assertGreaterEqual(point['certified_integer_bits'],232)
        result=certify(last=16,precision=1024)
        self.assertTrue(result['complete'])
        self.assertTrue(result['rigorous'])
        self.assertFalse(result['whole_code'])
        self.assertEqual(result['active_range'],[1,16])
        self.assertGreaterEqual(result['margin_bits_lower'],76)
        p,q=(int(x,16) for x in result['failure_upper_dyadic_hex'])
        self.assertLessEqual(F(p,q),F(1,2**76))


if __name__=='__main__':
    unittest.main()

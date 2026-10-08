from fractions import Fraction as F
from itertools import product
import unittest

from flint import ctx
from gf8_mds_packet import fullstate_matrices
from gf8_packet_maps import determinant
from gf8_outer_rate_dense import geometry, interval_margin as old_margin, rational
from gf8_punctured_dense import certify, interval_margin
from gf8_punctured_profile import metadata


class PuncturedDenseTests(unittest.TestCase):
    def test_exact_geometry_and_feedthrough(self):
        g = metadata()
        self.assertEqual((g['rows'],g['input_n'],g['output_n'],g['padded_k']),
                         (12544,1605632,1404928,1053696))
        self.assertEqual(F(g['padded_k'],g['output_n']),F(3,4))
        self.assertNotEqual(determinant(fullstate_matrices()[0]),0)

    def test_completion_count(self):
        for q in (2,3,4):
            for radius in (0,1,2):
                visible = [v for v in product(range(q),repeat=2) if sum(x!=0 for x in v)<=radius]
                completed = [v for v in product(range(q),repeat=4)
                             if sum(x!=0 for x in v[:2])<=radius]
                self.assertEqual(len(completed),len(visible)*q**2)

    def test_dense_margin_changes_only_output_union(self):
        previous = ctx.prec
        try:
            ctx.prec = 192
            g = metadata()
            old = geometry(84)
            actual = interval_margin(g,750000,752000,F(10))
            reference = old_margin(old,750000,752000,F(10))
            loss = F(673,100)*(g['input_n']-g['output_n'])
            difference = actual-reference
            self.assertTrue((difference+rational(loss)).contains(0))
        finally:
            ctx.prec = previous

    def test_directed_complete_replay(self):
        result = certify()
        self.assertEqual(result['weights'],[722535,1605632])
        self.assertTrue(result['complete'] and result['rigorous'])
        self.assertFalse(result['whole_code'])
        self.assertGreaterEqual(result['margin_bits_lower'],40)
        self.assertEqual(result['completion_factor_log2'],1605632)

    def test_reject_float_density(self):
        with self.assertRaises(ValueError):
            certify(.45)


if __name__=='__main__':
    unittest.main()

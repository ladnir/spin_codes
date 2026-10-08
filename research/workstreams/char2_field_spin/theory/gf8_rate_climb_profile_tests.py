from fractions import Fraction as F
import unittest

from gf8_punctured_profile import metadata as frozen
from gf8_rate_climb_profile import metadata


class ClimbProfileTests(unittest.TestCase):
    def test_checkpoint_identity(self):
        self.assertEqual(metadata(),frozen())

    def test_geometry(self):
        for outputs in (6,7,8):
            for dimension in (int(12*outputs),int(13*outputs)):
                for k_log in (16,18,20):
                    g=metadata(dimension,outputs,k_log)
                    self.assertEqual(g['input_n'],128*g['rows'])
                    self.assertEqual(g['output_n'],16*outputs*g['rows'])
                    self.assertEqual(g['padded_k'],dimension*g['rows'])
                    self.assertEqual(F(g['nominal_rate']),F(g['padded_k'],g['output_n']))
                    self.assertLess(10*g['output_threshold'],g['output_n'])
                    self.assertGreaterEqual(10*(g['output_threshold']+1),g['output_n'])
                    self.assertGreaterEqual(g['padded_k'],2**k_log)
                    self.assertEqual(g['rows']%64,0)

    def test_invalid(self):
        for args in ((0,7,20),(112,7,20),(84,9,20),(84,7,19)):
            with self.assertRaises(ValueError):metadata(*args)


if __name__=='__main__':unittest.main()

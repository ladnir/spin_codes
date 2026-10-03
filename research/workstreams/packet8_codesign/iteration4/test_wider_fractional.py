import unittest
from math import comb,log
import wider_fractional as wf


class WiderFractionalTests(unittest.TestCase):
    def test_complete_geometries(self):
        for bits,groups,regions in ((16,512,32),(32,256,64)):
            shape=wf.wider.geometry(symbol_bits=bits)
            self.assertEqual(shape['outer_groups'],groups)
            self.assertEqual(shape['regions'],regions)
            self.assertEqual(shape['physical_steps_total'],2048)
            self.assertEqual(groups*shape['group_dimension'],65536)
            self.assertEqual(regions*groups*8,131072)

    def test_union_is_not_powered(self):
        shape=wf.wider.geometry(symbol_bits=32)
        got=wf.exponent(shape,8,.06,.4,256*log(2),-30)
        expected=log(comb(256,8))+.4*(8*256*log(2)+.06*13107)-30
        self.assertEqual(got,expected)

    def test_bad_point_rejected(self):
        for s in ('1:.1:.4','2:0:.4','2:.1:0','2:.1:1.2'):
            with self.assertRaises(ValueError): wf.point(s)


if __name__=='__main__': unittest.main()

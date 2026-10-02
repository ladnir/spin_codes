from fractions import Fraction
import unittest
from flint import arb,arb_mat,ctx
import attack_cache as cache


class AttackCacheTests(unittest.TestCase):
    def test_roundtrip_exact_types_and_upper_endpoints(self):
        ctx.prec=192
        record=({(1,2):[Fraction(3,7),arb(1)/3]},arb_mat([[1,arb(2)/3]]))
        restored=cache.decode(cache.encode(record))
        self.assertEqual(restored[0][1,2][0],Fraction(3,7))
        self.assertTrue(restored[0][1,2][1].is_exact())
        self.assertGreaterEqual(restored[0][1,2][1],record[0][1,2][1])
        self.assertTrue(restored[1][0,1].is_exact())

    def test_reject_unknown_node(self):
        with self.assertRaises(ValueError):cache.decode({'python':'bad'})


if __name__=='__main__':unittest.main()

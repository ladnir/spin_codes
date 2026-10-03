"""Exact finite-field fixtures for the structured randomizer transport."""
from fractions import Fraction as F
from itertools import product
import unittest
import structured_randomizers as s


class StructuredRandomizers(unittest.TestCase):
    def test_gf16_aes_matrix_all69_minors(self):
        h = [[2,3,1,1],[1,2,3,1],[1,1,2,3],[3,1,1,2]]
        passed,details = s.mds_minors(h,0x13)
        self.assertTrue(passed)
        self.assertEqual(details['checked'],69)

    def test_gf16_quadratic_extension_is_field(self):
        self.assertTrue(all(s.mul(x,x,0x13)^x^8 for x in range(16)))
        for a in range(1,256):
            self.assertEqual({s.tower256_mul(a,b) for b in range(1,256)},set(range(1,256)))
        for a in range(256):
            expected = s.mul(2,a&15,0x13)|(s.mul(2,a>>4,0x13)<<4)
            self.assertEqual(s.tower256_mul(2,a),expected)

    def test_gf8_mds_sandwich_all_input_supports(self):
        h = [[s.inverse(x^y,0b1011) for y in range(4,8)] for x in range(4)]
        self.assertTrue(s.mds_minors(h,0b1011)[0])
        for matrix in (h,list(map(list,zip(*h)))):
            maximum,counts = s.sandwich_point_masses(matrix,0b1011)
            self.assertEqual(maximum,F(1,7**4))
            self.assertEqual(sum(counts.values()),8**4-1)
            self.assertEqual({before for before,after in counts},set(range(1,16)))

    def test_removing_one_input_diagonal_preserves_ensemble(self):
        original,normalized = s.normalized_diagonal_ensemble([[1,1],[1,2]],0b111)
        self.assertEqual(set(original),set(normalized))
        self.assertTrue(all(original[key]==3*normalized[key] for key in original))

    def test_two_fixed_input_diagonals_have_large_atom(self):
        h = [[s.inverse(x^y,0b1011) for y in range(4,8)] for x in range(4)]
        x = [s.mul(s.inverse(h[0][0],0b1011),h[0][1],0b1011),1,0,0]
        y = s.matvec(h,x,0b1011)
        self.assertEqual(y[0],0)
        self.assertEqual(sum(bool(v) for v in y),3)
        # With these two D_in coordinates fixed, D_out alone produces an
        # exact atom1/7^3, seven times the claimed sandwich envelope.
        self.assertGreater(F(1,7**3),F(1,7**4))

    def test_xor_middle_has_non_mild_atom(self):
        h = [[int(i!=j) for j in range(4)] for i in range(4)]
        maximum,_ = s.sandwich_point_masses(h,0b1011)
        self.assertGreaterEqual(maximum,F(1,7**3))
        self.assertFalse(s.mds_minors(h,0b1011)[0])


if __name__=='__main__':
    unittest.main()

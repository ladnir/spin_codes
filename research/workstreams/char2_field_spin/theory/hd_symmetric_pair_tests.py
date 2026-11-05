from collections import Counter
from fractions import Fraction as F
from itertools import product
import unittest

from hd_symmetric_pair import apply_pair, support_law, atom_majorant, positive_mask_envelope, primitive_counts
from hd_outer_pair_layers import field_multiply


class SymmetricPairTests(unittest.TestCase):
    def test_exact_GF4_value_law_and_involution(self):
        for a, b in product(range(4), repeat=2):
            mask = bool(a) | (bool(b) << 1)
            counts = Counter()
            for alpha in (2, 3):
                y = apply_pair(a, b, alpha, modulus=0b111)
                counts[bool(y[0]) | (bool(y[1]) << 1)] += 1
                self.assertEqual(apply_pair(*y, alpha, modulus=0b111), (a, b))
            law = support_law(mask, equal_nonzero=(a == b), q=4)
            self.assertEqual(tuple(F(counts[i], 2) for i in range(4)), law)

    def test_pointwise_product_majorant_and_general_envelope(self):
        for mask in range(4):
            for equal in (False, True):
                law = support_law(mask, equal_nonzero=equal, q=4)
                envelope = positive_mask_envelope(mask, q=4)
                self.assertTrue(all(x <= y for x, y in zip(law, envelope)))
                for t0, t1 in product((F(0), F(1, 3), F(1)), repeat=2):
                    actual = law[0]+law[1]*t0+law[2]*t1+law[3]*t0*t1
                    self.assertLessEqual(actual, atom_majorant(mask, t0, t1, q=4))
        self.assertEqual(sum(positive_mask_envelope(3)), 1+F(2, 254))

    def test_equality_and_marker_restrictions_are_material(self):
        with self.assertRaises(ValueError):
            support_law(3)
        with self.assertRaises(ValueError):
            atom_majorant(3, 2, 2)
        # At markers2,2 the unequalGF4 law weighs2, while equal inputs weigh4.
        law = support_law(3, equal_nonzero=False, q=4)
        self.assertEqual(law[1]*2+law[2]*2+law[3]*4, 2)

    def test_old_input_gauge_family_is_not_closed(self):
        # No output diagonal and alpha' restores M_2 diag(1,2) to M_alpha'.
        def mul(a, b):
            return field_multiply(a, b, modulus=0b111)
        target_columns = (apply_pair(1, 0, 2, modulus=0b111), apply_pair(0, 2, 2, modulus=0b111))
        candidates = []
        for alpha, d0, d1 in product((2, 3), range(1, 4), range(1, 4)):
            columns = tuple(tuple(mul(d, value) for d, value in zip((d0, d1), apply_pair(*x, alpha, modulus=0b111)))
                            for x in ((1, 0), (0, 1)))
            candidates.append(columns)
        self.assertNotIn(target_columns, candidates)

    def test_cost_and_invalid_key(self):
        self.assertEqual((primitive_counts()['variable_products'], primitive_counts()['binary_xors']), (1, 3))
        with self.assertRaises(ValueError):
            apply_pair(1, 1, 1)


if __name__ == '__main__':
    unittest.main()

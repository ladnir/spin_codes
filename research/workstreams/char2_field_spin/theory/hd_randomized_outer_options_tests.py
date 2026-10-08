"""Exhaustive audits of scalar entropy, row normalization, and descent limits."""

from fractions import Fraction
import unittest

from hd_randomized_outer_options import (
    base_quiet_dense_receipt, code_image_law, evaluation_code,
    randomized_measure, multiply, power,
)


class RandomizedOuterOptionsTests(unittest.TestCase):
    def test_arbitrary_subset_and_nonuniform_laws(self):
        for field, length, dimension, law in (
            (4, 3, 2, {1:Fraction(1, 2), 2:Fraction(1, 2)}),
            (4, 3, 2, {1:Fraction(1, 3), 2:Fraction(2, 3)}),
            (8, 3, 1, {1:Fraction(1, 3), 2:Fraction(1, 3), 3:Fraction(1, 3)}),
        ):
            code = evaluation_code(field, length, dimension)
            measure = randomized_measure(code, law, field)
            cap = max(law.values())**(length-dimension)
            self.assertLessEqual(max(measure.values()), cap)
            self.assertEqual(sum(measure.values()), field**dimension-1)
            self.assertNotIn((0,)*length, measure)

    def test_one_fixed_scalar_preserves_entire_code_image_law(self):
        code = evaluation_code(4, 3, 2)
        image_law = code_image_law(code, (1, 2, 3), 4)
        self.assertEqual(image_law,
                         code_image_law(code, (1, 2, 3), 4, fix_first=True))
        self.assertEqual(len(image_law), 3**2)
        self.assertEqual(set(image_law.values()), {Fraction(1, 3**2)})

    def test_subfield_orbit_does_not_give_full_field_envelope(self):
        generator = next(x for x in range(2, 16) if power(x, 3, 16) == 1)
        subgroup = (1, generator, generator ^ 1)
        code = evaluation_code(16, 3, 1)
        law = {x:Fraction(1, 3) for x in subgroup}
        measure = randomized_measure(code, law, 16)
        self.assertEqual(measure[(1, 1, 1)], Fraction(1, 9))
        self.assertGreater(measure[(1, 1, 1)], Fraction(1, 15**2))
        self.assertEqual(code_image_law(code, subgroup, 16),
                         code_image_law(code, subgroup, 16, fix_first=True))

    def test_shared_row_scalar_does_not_randomize_the_code(self):
        code = evaluation_code(4, 3, 2)
        for scalar in (1, 2, 3):
            scaled = {tuple(multiply(scalar, x, 4) for x in word) for word in code}
            self.assertEqual(scaled, set(code))

    def test_extension_coefficient_can_destroy_base_field_witnesses(self):
        generator = next(x for x in range(2, 16) if power(x, 3, 16) == 1)
        base = {0, 1, generator, generator ^ 1}
        eta = next(x for x in range(16) if x not in base)
        line = {(a, multiply(eta, a, 16)) for a in range(16)}
        base_words = {word for word in line if all(x in base for x in word)}
        self.assertEqual(base_words, {(0, 0)})
        self.assertTrue(any(any(word) for word in line))

    def test_directed_base_quiet_dense_tail(self):
        result = base_quiet_dense_receipt()
        self.assertTrue(result['directed'])
        self.assertEqual(result['first'], 7519)
        self.assertFalse(result['whole_code_certificate'])
        with self.assertRaises(ArithmeticError):
            base_quiet_dense_receipt(7518)


if __name__ == '__main__':
    unittest.main()

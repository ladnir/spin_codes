"""Exact arithmetic checks for the GF16 RS[16,8] reference circuit."""

import unittest

import rs16_maps as rs


def dot(a, b):
    return sum((x & y).bit_count() for x, y in zip(a, b)) & 1


class Rs16MapTests(unittest.TestCase):
    def test_field_inverse_and_distributivity(self):
        for a in range(1, 16):
            self.assertEqual(rs.multiply(a, rs.inverse(a)), 1)
        for a in range(16):
            for b in range(16):
                for c in range(16):
                    self.assertEqual(rs.multiply(a, b ^ c),
                                     rs.multiply(a, b) ^ rs.multiply(a, c))

    def test_systematic_interpolation(self):
        for degree in range(8):
            def value(x):
                result = 1
                for _ in range(degree):
                    result = rs.multiply(result, x)
                return result
            self.assertEqual(rs.encode_row([value(x) for x in range(8)]),
                             tuple(value(x) for x in range(16)))

    def test_factor_and_exact_counts(self):
        check = rs.check_parity()
        self.assertEqual(rs.PARITY_CONVOLUTION, (15, 2, 12, 5, 10, 4, 3, 8))
        self.assertEqual(rs.SQUARE_ZERO_COEFFICIENTS, (1, 11, 2, 13, 5, 12, 11, 8))
        self.assertEqual(check['direct_nontrivial_products'], 64)
        self.assertEqual(check['factored_nontrivial_products'], 19)
        self.assertEqual(check['change_of_basis_xors'], 24)
        self.assertEqual(rs.superset_zeta(rs.superset_zeta(list(range(8)))),
                         list(range(8)))

    def test_full_binary_adjoint(self):
        for source_bit in range(32):
            source = [0] * 8
            source[source_bit // 4] = 1 << (source_bit % 4)
            encoded = rs.encode_row(source)
            for target_bit in range(64):
                target = [0] * 16
                target[target_bit // 4] = 1 << (target_bit % 4)
                self.assertEqual(dot(encoded, target),
                                 dot(source, rs.transpose_row(target)))

    def test_adjoint_is_not_ordinary_multiplication(self):
        self.assertNotEqual(rs.adjoint_multiply(2, 1), rs.multiply(2, 1))
        for coefficient in range(16):
            for x in range(16):
                for y in range(16):
                    self.assertEqual(dot((rs.multiply(coefficient, x),), (y,)),
                        dot((x,), (rs.adjoint_multiply(coefficient, y),)))


if __name__ == '__main__':
    unittest.main()

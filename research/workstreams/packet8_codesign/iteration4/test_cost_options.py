"""Exact algebra checks for proposed low-cost maps; no timing benchmark."""
import random
import unittest

import cost_options as cost


class CostOptionsTests(unittest.TestCase):
    def test_cubic_has_no_gf256_root(self):
        for x in range(256):
            self.assertNotEqual(cost.mul8(cost.mul8(x, x), x)^x^1, 0)

    def test_six_products_match_polynomial_arithmetic(self):
        randomizer = random.Random(2217)
        for _ in range(512):
            left, right = randomizer.randrange(1 << 24), randomizer.randrange(1 << 24)
            self.assertEqual(cost.mul24_six(left, right), cost.mul24_reference(left, right))

    def test_scalar_binary_adjoint_circuit(self):
        randomizer = random.Random(2906)
        scalars = [1 << bit for bit in range(24)]+[randomizer.randrange(1, 1 << 24) for _ in range(16)]
        for scalar in scalars:
            columns = [cost.mul24_reference(1 << bit, scalar) for bit in range(24)]
            for row in range(24):
                expected = sum(((column >> row) & 1) << bit for bit, column in enumerate(columns))
                self.assertEqual(cost.adj24_six(1 << row, scalar), expected)

    def test_nonzero_scalar_inverse_samples(self):
        for scalar in (1, 2, 0x100, 0x10000, 0x123456, 0xffffff):
            inverse = cost.pow24(scalar, (1 << 24)-2)
            self.assertEqual(cost.mul24_six(scalar, inverse), 1)

    def test_two_band_fold_matches_binary_transpose(self):
        # Check all output basis vectors for every byte scale. Linearity then
        # proves the circuit for all64-bit inputs at each of those255 scales.
        for d in range(1, 256):
            scales = (1,)*4+(d,)*4
            columns = [cost.expansion24(1 << bit, scales) for bit in range(24)]
            for row in range(64):
                expected = sum(((column >> row) & 1) << bit for bit, column in enumerate(columns))
                self.assertEqual(cost.two_band_expansion_adjoint(1 << row, d), expected)

    def test_quadratic_field_and_three_product_circuit(self):
        for x in range(256):
            self.assertNotEqual(cost.mul8(x, x)^x^0x20, 0)
        randomizer = random.Random(4096)
        for _ in range(512):
            a, b = randomizer.randrange(65536), randomizer.randrange(65536)
            self.assertEqual(cost.mul16_three(a, b), cost.mul16_reference(a, b))

    def test_quadratic_binary_adjoint_circuit(self):
        for scalar in [1 << bit for bit in range(16)]+[0x1234, 0xffff, 0x2718]:
            columns = [cost.mul16_reference(1 << bit, scalar) for bit in range(16)]
            for row in range(16):
                expected = sum(((column >> row) & 1) << bit for bit, column in enumerate(columns))
                self.assertEqual(cost.adj16_three(1 << row, scalar), expected)

    def test_two_band_CA_formula(self):
        for d in (1, 2, 15, 0x53, 0xff):
            t = cost.mul8(6, 1^d)
            for bit in range(24):
                a, b, c = cost.bytes_of(1 << bit, 3)
                expected = cost.word_of((0, cost.mul8(t, c), cost.mul8(t, b)))
                self.assertEqual(cost.feedback24(cost.expansion24(1 << bit, (1,)*4+(d,)*4)), expected)

    def test_operation_counts_include_reverse_boundary(self):
        counts = cost.gfni_counts()
        self.assertEqual(counts['inner_gfni']['baseline16_GL2'], 65524)
        self.assertEqual(counts['inner_gfni']['state24_scalar'], 90092)
        self.assertEqual(counts['small_outer_gfni'], 96256)
        self.assertEqual(counts['wide_outer_shared_gfni'], 120832)
        self.assertEqual(counts['outer_randomizer_extra_xors'], 73728)
        self.assertEqual(counts['wide_plus_scalar24_extra_gfni'], 49144)
        self.assertEqual(counts['inner_gfni']['width32_state16_scalar'], 73720)
        self.assertEqual(counts['inner_gfni']['width32_state16_GL2'], 81910)


if __name__ == '__main__':
    unittest.main()

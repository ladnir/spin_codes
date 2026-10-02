"""Exact field and binary-adjoint checks; no timing or certificate claim."""
import random
import unittest
import field_symbol_randomizer as tower


class FieldSymbolRandomizerTests(unittest.TestCase):
    def test_quadratic_irreducibility(self):
        self.assertEqual(tower.trace(tower.MU, 8, tower.mul8), 1)
        self.assertEqual(tower.trace(tower.NU, 16, tower.mul16), 1)
        self.assertTrue(all(tower.mul8(x, x) ^ x != tower.MU for x in range(256)))

    def test_gf256_adjoint_exhaustive(self):
        for coefficient in range(256):
            matrix = tower.adjoint_matrix8(coefficient)
            images = [tower.mul8(coefficient, 1 << bit) for bit in range(8)]
            for value in range(256):
                expected = sum(((image & value).bit_count() & 1) << bit
                               for bit, image in enumerate(images))
                self.assertEqual(tower.affine(value, matrix), expected)

    def test_nine_map_factor_and_binary_adjoint(self):
        rng = random.Random(0x543332)
        scalars = [1 << bit for bit in range(32)] + [rng.randrange(1, 1 << 32) for _ in range(32)]
        for scalar in scalars:
            rows = tower.adjoint_rows(scalar)
            for value in [1 << bit for bit in range(32)] + [rng.getrandbits(32) for _ in range(8)]:
                expected = sum(((image & value).bit_count() & 1) << bit
                               for bit, image in enumerate(rows))
                self.assertEqual(tower.transpose32(value, scalar), expected)

    def test_field_inverse_and_frobenius(self):
        rng = random.Random(0x544F574552)
        for value in [1 << bit for bit in range(32)] + [rng.randrange(1, 1 << 32) for _ in range(32)]:
            self.assertEqual(tower.power(value, 1 << 32), value)
            self.assertEqual(tower.mul32(value, tower.power(value, (1 << 32) - 2)), 1)
        for _ in range(128):
            a, b, c = (rng.getrandbits(32) for _ in range(3))
            self.assertEqual(tower.mul32(a, b), tower.mul32(b, a))
            self.assertEqual(tower.mul32(a, b ^ c), tower.mul32(a, b) ^ tower.mul32(a, c))
            self.assertEqual(tower.mul32(a, tower.mul32(b, c)), tower.mul32(tower.mul32(a, b), c))

    def test_nonzero_contract(self):
        for scalar in (0, -1, 1 << 32):
            with self.assertRaises(ValueError):
                tower.adjoint_coefficients(scalar)

    def test_third_quadratic_irreducibility(self):
        self.assertEqual(tower.trace(tower.THETA, 32, tower.mul32), 1)

    def test_twenty_seven_map_factor_and_binary_adjoint(self):
        rng = random.Random(0x543634)
        scalars = [1 << bit for bit in range(64)] + [rng.randrange(1, 1 << 64) for _ in range(16)]
        for scalar in scalars:
            rows = tower.adjoint_rows64(scalar)
            self.assertEqual(len(tower.adjoint_coefficients64(scalar)), 27)
            for value in [1 << bit for bit in range(64)] + [rng.getrandbits(64) for _ in range(4)]:
                expected = sum(((image & value).bit_count() & 1) << bit
                               for bit, image in enumerate(rows))
                self.assertEqual(tower.transpose64(value, scalar), expected)

    def test_field64_inverse_and_frobenius(self):
        rng = random.Random(0x544F5745523634)
        for value in [1 << bit for bit in range(64)] + [rng.randrange(1, 1 << 64) for _ in range(8)]:
            self.assertEqual(tower.power(value, 1 << 64, tower.mul64), value)
            self.assertEqual(tower.mul64(value, tower.power(value, (1 << 64) - 2, tower.mul64)), 1)
        for _ in range(32):
            a, b, c = (rng.getrandbits(64) for _ in range(3))
            self.assertEqual(tower.mul64(a, b), tower.mul64(b, a))
            self.assertEqual(tower.mul64(a, b ^ c), tower.mul64(a, b) ^ tower.mul64(a, c))
            self.assertEqual(tower.mul64(a, tower.mul64(b, c)), tower.mul64(tower.mul64(a, b), c))

    def test_field64_nonzero_contract(self):
        for scalar in (0, -1, 1 << 64):
            with self.assertRaises(ValueError):
                tower.adjoint_coefficients64(scalar)

    def test_byte_factored_multiplier_and_rows(self):
        rng = random.Random(0x425954453332)
        scalars = [1 << bit for bit in range(32)] + [rng.randrange(1, 1 << 32) for _ in range(32)]
        for scalar in scalars:
            rows = tower.multiply_rows32(scalar)
            for value in [1 << bit for bit in range(32)] + [rng.getrandbits(32) for _ in range(8)]:
                expected = tower.mul32(value, scalar)
                self.assertEqual(tower.multiply32_factored(value, scalar), expected)
                self.assertEqual(sum(((row & value).bit_count() & 1) << j for j, row in enumerate(rows)), expected)

    def test_adjoint_family_scalar_action_is_bijective(self):
        # Forfixednonzero x, c -> M_c^T*x must havefullbinaryrank. This checks
        # the replacement's exact-transitivity argument onbasis+randominputs.
        rng = random.Random(0x5452414E534954)
        values = [1 << bit for bit in range(32)] + [rng.randrange(1, 1 << 32) for _ in range(32)]
        for value in values:
            pivots = {}
            for bit in range(32):
                image = tower.transpose32(value, 1 << bit)
                while image:
                    pivot = image.bit_length() - 1
                    if pivot not in pivots:
                        pivots[pivot] = image
                        break
                    image ^= pivots[pivot]
            self.assertEqual(len(pivots), 32)

    def test_byte_factored_multiplier64_and_rows(self):
        rng = random.Random(0x425954453634)
        scalars = [1 << bit for bit in range(64)] + [rng.randrange(1, 1 << 64) for _ in range(8)]
        for scalar in scalars:
            rows = tower.multiply_rows64(scalar)
            for value in [1 << bit for bit in range(64)] + [rng.getrandbits(64) for _ in range(4)]:
                expected = tower.mul64(value, scalar)
                self.assertEqual(tower.multiply64_factored(value, scalar), expected)
                self.assertEqual(sum(((row & value).bit_count() & 1) << j for j, row in enumerate(rows)), expected)


if __name__ == "__main__":
    unittest.main()

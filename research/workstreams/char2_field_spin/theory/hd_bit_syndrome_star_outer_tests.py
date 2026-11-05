"""Small exact checks only; no regional provider or native benchmark."""
import random
import unittest

import hd_bit_syndrome_outer as old
import hd_bit_syndrome_star_outer as new
from gf8_packet_maps import INV, MUL, determinant, matmul, matvec


class StarOuterTests(unittest.TestCase):
    def test_factorization_involution_and_no_zero_entries(self):
        for variant in new.QUADRATIC_COEFFICIENTS:
            v, b, h = new.syndrome_rows(), new.coefficient_matrix(variant), new.matrix(variant)
            self.assertEqual(determinant(b), determinant(old.coefficient_matrix()))
            self.assertNotEqual(determinant(b), 0)
            self.assertTrue(all(x == 0 for row in matmul(tuple(zip(*v)), v) for x in row))
            n = matmul(matmul(v, b), tuple(zip(*v)))
            self.assertEqual(h, tuple(tuple(n[i][j] ^ int(i == j) for j in range(16)) for i in range(16)))
            self.assertEqual(matmul(h, h), tuple(tuple(int(i == j) for j in range(16)) for i in range(16)))
            self.assertTrue(all(x for row in h for x in row))
        c = new.quadratic_coefficients('coupled')
        self.assertEqual(determinant(c), 1)
        self.assertEqual(determinant(tuple(tuple(c[i][j] ^ int(i == j)
                                                 for j in range(3)) for i in range(3))), 1)
        self.assertTrue(all(x != 1 for row in old.matrix() for x in row))

    def test_shared_circuit_basis_random_and_operation_count(self):
        for variant in new.QUADRATIC_COEFFICIENTS:
            for i in range(16):
                for bit in range(8):
                    x = [0] * 16
                    x[i] = 1 << bit
                    self.assertEqual(new.circuit(x, variant), new.direct(x, variant))
            rng = random.Random(20261009)
            for _ in range(32):
                x = [rng.randrange(256) for _ in range(16)]
                self.assertEqual(new.circuit(x, variant), new.direct(x, variant))
        self.assertEqual(new.circuit_counts(), {'binary_xors': 74, 'fixed_products': 8})
        self.assertEqual(new.circuit_counts('coupled'), {'binary_xors': 76, 'fixed_products': 8})
        self.assertEqual(old.circuit_counts(), {'binary_xors': 66, 'fixed_products': 8})

    def test_rank_is_direct_sum_of_two_halves(self):
        rows = tuple(sum(v << j for j, v in enumerate(row)) for row in new.syndrome_rows())
        # All half supports, plus deterministic mixed supports, test the
        # actual eight-column matrix rather than reusing its split formula.
        masks = {sum(((m >> i) & 1) << (2 * i + p) for i in range(8))
                 for m in range(256) for p in (0, 1)}
        masks.update(random.Random(33).sample(range(1 << 16), 512))
        for mask in masks:
            self.assertEqual(new.syndrome_rank(mask), new._binary_rank(
                rows[i] for i in range(16) if mask >> i & 1))
            self.assertGreaterEqual(new.syndrome_rank(mask), old.affine_rank(mask))
        self.assertEqual(new.syndrome_rank(0xffff), 8)
        self.assertEqual(new.syndrome_rank(0x5555), 4)
        self.assertEqual(new.syndrome_rank(0xaaaa), 4)

    def test_exact_binary_quiet_enumerator_by_brute_words(self):
        rows = new.syndrome_rows()
        check = [0] * 17
        for mask in range(1 << 16):
            syndrome = 0
            for i in range(16):
                if mask >> i & 1:
                    syndrome ^= sum(x << j for j, x in enumerate(rows[i]))
            if syndrome == 0:
                check[mask.bit_count()] += 1
        self.assertEqual(tuple(check), new.quiet_weight_enumerator(2))
        self.assertEqual(sum(check), 2 ** 8)
        self.assertEqual(next(i for i in range(1, 17) if check[i]), 4)

    def test_histograms_and_quiet_subset_domination(self):
        from math import comb
        expected = {
            7: {4: 16, 5: 448, 6: 4704, 7: 6272},
            8: {4: 2, 5: 128, 6: 1764, 7: 7840, 8: 3136},
            9: {5: 16, 6: 448, 7: 4704, 8: 6272},
        }
        # Each half has 14 affine-plane four-sets and 56 rank-four
        # four-sets; combining these gives the independently derived table.
        for a in expected:
            self.assertEqual(new.rank_histogram(a), expected[a])
            self.assertEqual(sum(new.rank_histogram(a).values()), comb(16, a))
        current, previous = new.quiet_weight_enumerator(), new.original_quiet_weight_enumerator()
        self.assertEqual(sum(current), 256 ** 8)
        self.assertEqual(sum(previous), 256 ** 11)
        self.assertTrue(all(0 <= a <= b for a, b in zip(current, previous)))
        # A full half retains the old four-dimensional quiet subspace.
        self.assertEqual(new.quiet_support_count(0x5555), old.quiet_support_counts()[0x5555])

    def test_projective_pair_law_and_involution_duality(self):
        for i, j in ((0, 1), (0, 15), (1, 3), (7, 15)):
            law = new.pair_mask_counts(i, j)
            self.assertEqual(sum(count for _, count in law), 255)
            self.assertTrue(all(mask for mask, _ in law))
        h = new.matrix()
        # Any weight-two output gives an input with the same tabulated mask
        # after inversion; it is nonquiet whenever that mask has weight>2.
        y = [0] * 16
        y[1], y[3] = 1, 1
        x = matvec(h, y)
        self.assertEqual(new.direct(x), tuple(y))
        self.assertNotEqual(tuple(x), tuple(y))

    def test_complete_pair_screen_records_deterioration(self):
        result = new.pair_summary()
        self.assertEqual(result['ratio_count'], 30600)
        self.assertEqual(result['maximum_pair_holes'], 8)
        self.assertEqual(result['input_arities_for_nonquiet_weight_two_output'],
                         {8: 24, 10: 44, 12: 34, 14: 21, 15: 1286, 16: 29191})
        y = [0] * 16
        y[3], y[5] = 1, MUL[new.LABELS[3]][INV[new.LABELS[5]]]
        x = new.direct(y)
        self.assertEqual(sum(bool(v) for v in x), 8)
        self.assertEqual(sum(bool(v) << i for i, v in enumerate(x)), 5437)
        self.assertEqual(new.direct(x), tuple(y))
        coupled = new.pair_summary('coupled')
        self.assertEqual(coupled['ratio_count'], 30600)
        self.assertEqual(coupled['maximum_pair_holes'], 6)
        self.assertEqual(coupled['input_arities_for_nonquiet_weight_two_output'],
                         {10: 17, 12: 71, 14: 37, 15: 1460, 16: 29015})
        self.assertNotIn(8, coupled['input_arities_for_nonquiet_weight_two_output'])
        self.assertNotEqual(new.construction_identity(), new.construction_identity('coupled'))

    def test_guards(self):
        for mask in (-1, 1 << 16, True):
            with self.assertRaises(ValueError):
                new.syndrome_rank(mask)
        for q in (0, 3, True):
            with self.assertRaises(ValueError):
                new.quiet_weight_enumerator(q)
        with self.assertRaises(ValueError):
            new.circuit([0] * 15)
        with self.assertRaises(ValueError):
            new.pair_mask_counts(3, 3)
        with self.assertRaises(ValueError):
            new.matrix('unknown')


if __name__ == '__main__':
    unittest.main()

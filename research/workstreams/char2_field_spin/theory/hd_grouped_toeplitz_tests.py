import random
import unittest
from fractions import Fraction
from itertools import product
from hd_grouped_toeplitz import binary_field, matrix, nonsingular, circuit, toy_counts, column_cap, group_cap
from hd_grouped_outer import atom_cap


class ToeplitzColumnTests(unittest.TestCase):
    def test_all_toy_atoms_and_nonsingularity_counts(self):
        for q, g in ((2, 2), (2, 3), (2, 4), (4, 2), (4, 3), (4, 4)):
            rows, accepted = toy_counts(q, g)
            self.assertEqual(accepted, (q-1)*q**(2*g-2))
            cap = column_cap(q, g)
            for row in rows[1:]:
                self.assertEqual(row[0], 0)
                self.assertEqual(sum(row), accepted)
                self.assertLessEqual(Fraction(max(row), accepted), cap)
        rows, accepted = toy_counts(4, 4)
        self.assertEqual(accepted, 12288)
        self.assertEqual(max(map(max, rows[1:])), 64)

    def test_recursive_circuit_all_gf4_coefficients_on_basis(self):
        mul = binary_field(4)
        for g in (2, 4):
            for t in product(range(4), repeat=2*g-1):
                m = matrix(t)
                # Both circuits are field-linear, so these determine the map.
                for j in range(g):
                    x = [int(i == j) for i in range(g)]
                    self.assertEqual(circuit(t, x, mul), tuple(row[j] for row in m))

    def test_gf256_arithmetic_including_zero_coefficients(self):
        mul = binary_field(256)
        rng = random.Random(113)
        for g in (2, 4):
            for _ in range(256):
                t = [rng.randrange(256) for _ in range(2*g-1)]
                x = [rng.randrange(256) for _ in range(g)]
                expected = []
                for row in matrix(t):
                    y = 0
                    for a, b in zip(row, x):
                        y ^= mul[a][b]
                    expected.append(y)
                self.assertEqual(circuit(t, x, mul), tuple(expected))

    def test_old_whole_group_envelope_still_dominates(self):
        for q in (4, 256, 65536):
            for g in (2, 4):
                new = group_cap(q, g, 256//g, 128//g)
                old = atom_cap(q, g, 256//g, 128//g)
                self.assertLess(new, old)
                self.assertEqual(new / old, Fraction(q-1, q)**(128-128//g))

    def test_complete_two_row_rs_group(self):
        counts, den = toy_counts(4, 2)
        support = [tuple((y, c) for y, c in enumerate(row) if c) for row in counts]
        mul = binary_field(4)
        atoms = [0]*65536
        for a, b, c, d in product(range(4), repeat=4):
            if not (a or b or c or d):
                continue
            columns = [(a ^ mul[b][j]) + 4*(c ^ mul[d][j]) for j in range(4)]
            for y0, c0 in support[columns[0]]:
                for y1, c1 in support[columns[1]]:
                    low = y0 + 16*y1
                    mass = c0*c1
                    for y2, c2 in support[columns[2]]:
                        mid = low + 256*y2
                        m = mass*c2
                        for y3, c3 in support[columns[3]]:
                            atoms[mid + 4096*y3] += m*c3
        self.assertEqual(atoms[0], 0)
        self.assertEqual(sum(atoms), 255*den**4)
        self.assertLessEqual(Fraction(max(atoms), den**4), group_cap(4, 2, 4, 2))

    def test_invertibility_is_essential(self):
        mul = binary_field(4)
        self.assertFalse(nonsingular((0,)*7, mul))
        self.assertTrue(nonsingular((0, 0, 0, 1, 0, 0, 0), mul))
        # Accepting singular mixers would invalidate the shortened-information
        # proof's inference that an output zero column has a zero preimage.
        self.assertEqual(circuit((0,)*7, (1, 0, 0, 0), mul), (0,)*4)

    def test_invalid_parameters(self):
        for args in ((1, 4), (4, 0)):
            with self.assertRaises(ValueError):
                column_cap(*args)
        with self.assertRaises(ValueError):
            group_cap(256, 4, 32, 32)


if __name__ == '__main__':
    unittest.main()

"""Independent algebra checks for the input-weighted GRS[8,4] baseline."""

from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import arithmetic as a


def rank(rows, modulus):
    """Fraction-free elimination over the field, independent of inversion."""
    rows = [list(row) for row in rows]
    height = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(height, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[height], rows[pivot] = rows[pivot], rows[height]
        pivot_value = rows[height][column]
        for i in range(height + 1, len(rows)):
            scale = rows[i][column]
            for j in range(column + 1, len(rows[i])):
                rows[i][j] = (a.multiply(pivot_value, rows[i][j], modulus)
                              ^ a.multiply(scale, rows[height][j], modulus))
            rows[i][column] = 0
        height += 1
        if height == len(rows):
            break
    return height


def generator(modulus):
    rows = []
    for x in range(8):
        s1 = a.multiply(x, x, modulus) ^ x
        rows.append((1, x, s1, a.multiply(x, s1, modulus)))
    return rows


class ReturnAudit(unittest.TestCase):
    def test_every_grs_information_set(self):
        for modulus in a.MODULI.values():
            rows = generator(modulus)
            for selection in combinations(range(8), 4):
                self.assertEqual(rank([rows[j] for j in selection], modulus), 4)

    def test_every_shortened_dimension(self):
        # This directly checks all 256 coordinate supports over all four fields.
        for modulus in a.MODULI.values():
            rows = generator(modulus)
            for mask in range(256):
                zeros = [rows[j] for j in range(8) if not (mask >> j) & 1]
                self.assertEqual(4 - rank(zeros, modulus), max(0, mask.bit_count() - 4))

    def test_pair_evaluation_identity(self):
        for modulus in a.MODULI.values():
            for u in ((1, 2, 4, 8), (0, 0, 0, 1), (19, 77, 201, 18)):
                for v in (0, 2, 4, 6):
                    s1 = a.multiply(v, v, modulus) ^ v
                    low = u[0] ^ a.multiply(s1, u[2], modulus)
                    high = u[1] ^ a.multiply(s1, u[3], modulus)
                    even = low ^ a.multiply(v, high, modulus)
                    rows = generator(modulus)
                    literal = []
                    for point in (v, v ^ 1):
                        result = 0
                        for message, coefficient in zip(u, rows[point]):
                            result ^= a.multiply(message, coefficient, modulus)
                        literal.append(result)
                    self.assertEqual(literal, [even, even ^ high])

    def test_generic_return_factor_in_characteristic_two(self):
        modulus, coefficients = 0x7, (1, 2, 3)
        forms = [form for form in product(range(4), repeat=2) if form != (0, 0)]
        for ell, state in product(forms, repeat=2):
            for weight_input in (False, True):
                zero_forms = 0
                for coefficient in coefficients:
                    weighted, unweighted = (ell, state) if weight_input else (state, ell)
                    result = tuple(a.multiply(coefficient, x, modulus) ^ y
                                   for x, y in zip(weighted, unweighted))
                    zero_forms += result == (0, 0)
                self.assertLessEqual(zero_forms, 1)
                for tau in (F(1), F(3, 2), F(512)):
                    average = (zero_forms + F(3 - zero_forms, 1) / tau) / 3
                    self.assertLessEqual(average, F(1, 3) + F(2, 3) / tau)


if __name__ == "__main__":
    unittest.main()

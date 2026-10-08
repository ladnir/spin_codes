"""Exact small-field checks for HIGH_DISTANCE_GATE.md; no timing code."""

import itertools
import unittest
from fractions import Fraction


def multiply(a, b, modulus, degree):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a >> degree:
            a ^= modulus
    return result


def inverse(a, modulus, degree):
    result, power = 1, (1 << degree) - 2
    while power:
        if power & 1:
            result = multiply(result, a, modulus, degree)
        a = multiply(a, a, modulus, degree)
        power >>= 1
    return result


def rank(rows, modulus, degree):
    rows = [list(row) for row in rows]
    pivot = 0
    for column in range(len(rows[0])):
        candidate = next((i for i in range(pivot, len(rows)) if rows[i][column]), None)
        if candidate is None:
            continue
        rows[pivot], rows[candidate] = rows[candidate], rows[pivot]
        scale = inverse(rows[pivot][column], modulus, degree)
        rows[pivot] = [multiply(value, scale, modulus, degree) for value in rows[pivot]]
        for i in range(pivot + 1, len(rows)):
            scale = rows[i][column]
            rows[i] = [a ^ multiply(scale, b, modulus, degree)
                       for a, b in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def predicted_histogram(field_degree, coefficient_degree):
    if field_degree < 2 * coefficient_degree - 1:
        raise ValueError("polynomial products may wrap modulo the field polynomial")
    q = 1 << coefficient_degree
    result = {1: (1 << field_degree) - q * q // 2, q: 1}
    for height in range(1, coefficient_degree):
        result[1 << (coefficient_degree - height)] = 3 * (1 << (2 * height - 1))
    return {size: count for size, count in result.items() if count}


def projective_moment(field_degree, coefficient_degree, rows):
    q = 1 << coefficient_degree
    result = Fraction(2, q ** rows)
    for size, count in predicted_histogram(field_degree, coefficient_degree).items():
        result += count * Fraction(size, q * q) ** rows
    return result


def defect_bound(field_degree, coefficient_degree, rows, tau=None):
    field_size = 1 << field_degree
    if tau is None:
        tau = field_size
    q = 1 << coefficient_degree
    return (1 + (tau - 1) * projective_moment(field_degree, coefficient_degree, rows)
            + Fraction(tau * tau - tau, q ** (2 * rows)))


def affine_values(offset, first, second, coefficient_degree, modulus, field_degree):
    q = 1 << coefficient_degree
    return [tuple(z ^ multiply(a, x, modulus, field_degree)
                  ^ multiply(b, y, modulus, field_degree)
                  for z, x, y in zip(offset, first, second))
            for a in range(q) for b in range(q)]


class PolynomialBoxRankMomentTests(unittest.TestCase):
    def test_exact_intersections_and_moments(self):
        for field_degree, modulus in ((2, 0x7), (4, 0x13), (8, 0x11B)):
            for coefficient_degree in range(1, field_degree // 2 + 1):
                q = 1 << coefficient_degree
                histogram = {}
                for scalar in range(1, 1 << field_degree):
                    size = sum(multiply(scalar, value, modulus, field_degree) < q
                               for value in range(q))
                    histogram[size] = histogram.get(size, 0) + 1
                self.assertEqual(histogram, predicted_histogram(field_degree, coefficient_degree))
                for rows in (2, 3, 4):
                    actual = Fraction(2, q ** rows) + sum(
                        count * Fraction(size, q * q) ** rows for size, count in histogram.items())
                    self.assertEqual(actual, projective_moment(field_degree, coefficient_degree, rows))

    def test_closed_second_and_third_moments(self):
        for field_degree in (8, 32, 64, 128):
            field_size = 1 << field_degree
            for coefficient_degree in range(1, field_degree // 2 + 1):
                q = 1 << coefficient_degree
                second = Fraction(2 + 3 * coefficient_degree, 2 * q * q) + Fraction(field_size, q ** 4)
                third = Fraction(9, 2 * q ** 3) - Fraction(7, 2 * q ** 4) + Fraction(field_size, q ** 6)
                self.assertEqual(second, projective_moment(field_degree, coefficient_degree, 2))
                self.assertEqual(third, projective_moment(field_degree, coefficient_degree, 3))

    def test_affine_batches_with_offsets_and_different_planes(self):
        # All 16^r coefficient assignments are enumerated for each fixture.
        # Every pair of feedback forms is independent before the whole batch.
        field_degree, coefficient_degree, modulus = 4, 2, 0x13
        families = (
            ((0, 0, 0), (1, 0, 0), (0, 1, 0)),
            ((7, 2, 0), (1, 1, 0), (0, 2, 0)),
            ((0, 1, 0), (0, 1, 0), (0, 0, 1)),
            ((2, 0, 8), (1, 0, 1), (0, 1, 3)),
        )
        fixtures = ((0, 0), (0, 1), (0, 2), (1, 3),
                    (0, 0, 0), (0, 1, 2), (1, 2, 3))
        for indices in fixtures:
            value_sets = []
            for index in indices:
                offset, first, second = families[index]
                self.assertEqual(rank((first, second), modulus, field_degree), 2)
                value_sets.append(affine_values(offset, first, second, coefficient_degree,
                                                modulus, field_degree))
            counts = [0, 0, 0, 0]
            for rows in itertools.product(*value_sets):
                counts[rank(rows, modulus, field_degree)] += 1
            for tau in (1, 2, 4, 16):
                actual = sum(count * Fraction(tau) ** (2 - row_rank)
                             for row_rank, count in enumerate(counts)) / (16 ** len(indices))
                self.assertLessEqual(actual, defect_bound(field_degree, coefficient_degree, len(indices), tau))

    def test_native_half_width_bounds(self):
        for field_degree in (8, 32, 64, 128):
            coefficient_degree = field_degree // 2
            q = 1 << coefficient_degree
            self.assertLess(defect_bound(field_degree, coefficient_degree, 2),
                            Fraction(8 + 3 * coefficient_degree, 2))
            self.assertLess(defect_bound(field_degree, coefficient_degree, 3) - 1, Fraction(9, 2 * q))

    def test_degree_condition_is_required(self):
        with self.assertRaises(ValueError):
            predicted_histogram(4, 3)


if __name__ == "__main__":
    unittest.main()

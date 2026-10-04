"""Small exact checks for the outward replay; no state census or benchmark."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

import replay_fractional as replay


def exact(value):
    """Convert an exact Arb endpoint or FLINT rational without floating point."""
    if isinstance(value, arb):
        if not value.is_exact():
            raise ValueError('an exact endpoint is required')
        value = value.fmpq()
    return F(str(value))


def rational(value):
    value = F(value)
    return fmpq(value.numerator, value.denominator)


class OutwardReplayTests(unittest.TestCase):
    def setUp(self):
        self.previous_precision = ctx.prec
        ctx.prec = 256

    def tearDown(self):
        ctx.prec = self.previous_precision

    def test_exact_dyadic_union_retains_arbitrarily_small_terms(self):
        def decode(record):
            return F(int(record['numerator_hex'], 16)) * F(2) ** record['exponent']

        values = [F(1, 1 << 60), F(3, 1 << 10000), F(7, 1 << 301)]
        total = sum(values, F(0))
        encoded = replay.encode_exact(total)
        self.assertEqual(decode(encoded), total)
        self.assertGreater(decode(encoded), values[0])
        for value in (F(0), F(11 << 400), F(3, 1 << 10000)):
            self.assertEqual(decode(replay.encode_exact(value)), value)
        with self.assertRaises(ArithmeticError):
            replay.encode_exact(F(1, 3))

    def test_fractional_upper_by_exact_integer_comparison(self):
        # y >= x^(4/5) iff y^5 >= x^4 for nonnegative rational x,y.
        # This comparison does not use another approximate fractional power.
        exponent = replay.aq(F(4, 5))
        self.assertLessEqual(exact(exponent.lower()), F(4, 5))
        self.assertGreaterEqual(exact(exponent.upper()), F(4, 5))
        previous = F(0)
        for value in (F(1, 1 << 10000), F(1, 32), F(3, 8), F(1), F(5, 2), F(1 << 10000)):
            enclosure = replay.aq(value) ** exponent
            bound = exact(replay.upper(enclosure))
            self.assertGreater(bound, 0)
            self.assertGreaterEqual(bound ** 5, value ** 4)
            self.assertLessEqual(bound ** 5, value ** 4 * (1 + F(1, 1 << 220)))
            self.assertGreaterEqual(bound, previous)
            previous = bound

    def test_loaded_tiny_dyadics_are_exact_at_replay_precision(self):
        raw = [[[[0, 0] for _ in range(11)] for _ in range(11)] for _ in range(17)]
        raw[0][0][0] = [1, 0]
        raw[1][4][7] = [1, -10000]
        raw[2][3][6] = [(1 << 191) + 1, -12000]
        record = dict(schema='k18-certified22-outward-physical-export-1',
            source_sha256={}, map_record=dict(map_sha256=replay.MAP), precision=192,
            matrices=raw, tilt='17/200')
        # Mock only file reading: the fixture tests the numerical import,
        # independently of the real receipt's separate source authentication.
        with patch.object(Path, 'read_text', return_value=json.dumps(record)):
            matrices, tilt, _ = replay.load_local(Path('numeric-fixture.json'))
        self.assertEqual(tilt, F(17, 200))
        for k, i, j in ((1, 4, 7), (2, 3, 6)):
            self.assertGreater(matrices[k][i, j], 0)
            self.assertTrue(matrices[k][i, j].is_exact())
            self.assertEqual(exact(matrices[k][i, j]), replay.fraction(raw[k][i][j]))

    def test_alpha_one_matches_global_thinning_after_active_convolution(self):
        first = fmpq_mat([[1, fmpq(1, 4)], [0, fmpq(1, 2)]])
        second = fmpq_mat([[fmpq(1, 8), 0], [fmpq(1, 4), fmpq(3, 8)]])
        self.assertNotEqual(first * second, second * first)
        physical = [(first if j % 2 == 0 else second) / (j + 1) for j in range(17)]
        active = []
        # Independently convolve active occupancies first. Globally thin a
        # uniform potential subset afterward, reversing the producer's order
        # of probabilistic averaging while preserving the matrix chronology.
        for total in range(33):
            value = fmpq_mat(2, 2)
            for left in range(max(0, total - 16), min(16, total) + 1):
                coefficient = rational(F(comb(16, left) * comb(16, total - left), comb(32, total)))
                value += coefficient * (physical[left] * physical[total - left])
            active.append(value)
        actual = replay.fractional_blocks([arb_mat(matrix) for matrix in physical], F(1))
        for total, bound in enumerate(actual):
            truth = fmpq_mat(2, 2)
            for active_count in range(total + 1):
                probability = rational(F(comb(total, active_count) * 15 ** active_count, 16 ** total))
                truth += probability * active[active_count]
            for i in range(2):
                for j in range(2):
                    difference = exact(bound[i, j]) - exact(truth[i, j])
                    self.assertGreaterEqual(difference, 0)
                    self.assertLess(difference, F(1, 1 << 220))

    def test_fractional_blocks_preserve_tiny_positive_perfect_powers(self):
        for exponent in (-5, -5000):
            value = F(2) ** exponent
            physical = [arb_mat([[replay.aq(value)]]) for _ in range(17)]
            result = replay.fractional_blocks(physical, F(4, 5))
            truth = F(2) ** (2 * exponent * 4 // 5)
            for matrix in result:
                bound = exact(matrix[0, 0])
                self.assertGreater(bound, 0)
                self.assertGreaterEqual(bound, truth)
                self.assertLessEqual(bound, truth * (1 + F(1, 1 << 220)))

    def test_polynomial_coefficient_preserves_chronological_order(self):
        a = fmpq_mat([[1, 2], [0, 1]])
        b = fmpq_mat([[1, 0], [3, 1]])
        c = fmpq_mat([[2, 0], [1, 1]])
        d = fmpq_mat([[1, 4], [0, 3]])
        actual = replay.placement._product([a, b], [c, d], 2, fmpq_mat, lambda value: value)
        self.assertEqual(actual, [a * c, a * d + b * c, b * d])
        self.assertNotEqual(actual[1], d * a + c * b)


if __name__ == '__main__':
    unittest.main()

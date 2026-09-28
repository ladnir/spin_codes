"""Exhaustive small binary-map tests for the complete density-column bound."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement, permutations, product
from math import lcm
import unittest

from flint import arb, arb_mat, ctx

from column_density import density_caps, exact_bound, bound, refine
from occupancy_memory import M, C, U
from mature_tail import L48, L56


def linear(columns, value):
    result = 0
    for i, column in enumerate(columns):
        if value >> i & 1:
            result ^= column
    return result


def packets(shape, windows, feedback):
    """Uniform labeled distinct windows and uniform two-bit weight masks."""
    for positions in permutations(range(windows), len(shape)):
        choices = [tuple(mask for mask in (1, 2, 3) if mask.bit_count() == b) for b in shape]
        for masks in product(*choices):
            word = sum(mask << (2*position) for mask, position in zip(masks, positions))
            yield word, linear(feedback, word)


def columns(shape, windows, expansion, feedback, z):
    entries = list(packets(shape, windows, feedback))
    return [sum((z**((linear(expansion, t ^ syndrome) ^ word).bit_count())
                 for word, syndrome in entries), Q(0))/len(entries)
            for t in range(1 << len(expansion))]


def window_records(windows, expansion, feedback, z):
    result = {}
    for b in (1, 2):
        values = columns((b,), windows, expansion, feedback, z)
        nonzero, zero = max(values[1:]), values[0]
        denominator = lcm(nonzero.denominator, zero.denominator)
        # Mass is unused by this refinement; the trivial upper 1 suffices.
        result[b] = (denominator, int(nonzero*denominator), int(zero*denominator), denominator)
    return result


class ColumnDensity(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 160

    def tearDown(self):
        ctx.prec = self.precision

    def test_all_states_shapes_and_distinct_slots(self):
        cases = [(3, (0b010101, 0b001111), (1, 2, 3, 1, 3, 2)),
                 (4, (0b01010101, 0b00110011), (1, 2, 3, 1, 2, 3, 3, 1)),
                 (3, (0, 0), (1, 2, 3, 1, 3, 2)),
                 (4, (255,), (0,)*8)]
        for windows, expansion, feedback in cases:
            for z in (Q(1, 2), Q(4, 5), Q(1)):
                caps = density_caps(window_records(windows, expansion, feedback, z))
                for j in range(2, windows+1):
                    for shape in combinations_with_replacement((1, 2), j):
                        upper = exact_bound(shape, caps, z, windows=windows, packet_bits=2)
                        actual = columns(shape, windows, expansion, feedback, z)
                        self.assertTrue(all(value <= upper <= 1 for value in actual))
                        self.assertEqual(upper, exact_bound(shape[::-1], caps, z, windows=windows, packet_bits=2))

    def test_arbitrary_and_uniform_source_measures(self):
        windows, expansion = 3, (0b010101, 0b001111)
        feedback, z, alpha = (1, 2, 3, 1, 3, 2), Q(3, 4), Q(1, 4)
        caps = density_caps(window_records(windows, expansion, feedback, z))
        for shape in ((1, 1), (1, 2), (2, 2), (1, 1, 2)):
            upper = exact_bound(shape, caps, z, windows=windows, packet_bits=2)
            entries = list(packets(shape, windows, feedback))
            measures = [(Q(0), Q(1, 8), Q(3, 8), Q(1, 4))]
            for subset in ((1,), (1, 2), (2, 3), (1, 2, 3)):
                measures.append(tuple(Q(int(s in subset), len(subset)) for s in range(4)))
            for source in measures:
                for target in range(1, 4):
                    value = sum((source[target ^ syndrome]*z**((linear(expansion, target ^ syndrome)^word).bit_count())
                                 for word, syndrome in entries), Q(0))*alpha/len(entries)
                    self.assertLessEqual(value, alpha*max(source)*upper)

    def test_zero_target_is_required_in_complete_column_cap(self):
        z = Q(1, 2)
        records = window_records(4, (255,), (0,)*8, z)
        correct = exact_bound((1, 1), density_caps(records), z, windows=4, packet_bits=2)
        wrong_caps = {b: Q(record[1], record[3]) for b, record in records.items()}
        wrong = exact_bound((1, 1), wrong_caps, z, windows=4, packet_bits=2)
        zero_column = columns((1, 1), 4, (255,), (0,)*8, z)[0]
        self.assertLess(wrong, zero_column)
        self.assertLessEqual(zero_column, correct)
        self.assertGreater(density_caps(records)[1], wrong_caps[1])

    def test_outward_arithmetic(self):
        records = {1: (100, 3, 5, 100), 2: (100, 9, 2, 100)}
        shape, tilt = (1, 1, 2), Q(13, 250)
        for precision in (32, 64, 192):
            ctx.prec = precision
            upper = bound(shape, records, tilt, windows=4, packet_bits=2)
            ctx.prec = 512
            exact_scale = arb(tilt.numerator)/tilt.denominator
            candidates = [arb(1)]
            for b, cap in density_caps(records).items():
                candidates.append(arb(2)*arb(cap.numerator)/cap.denominator
                                  *(exact_scale*(sum(shape)-b)).exp())
            self.assertTrue(upper >= min(candidates))
            self.assertTrue(upper.is_exact())
        zero_tilt = bound(shape, records, 0, windows=4, packet_bits=2)
        endpoint = zero_tilt.fmpq()
        self.assertGreaterEqual(Q(int(endpoint.p), int(endpoint.q)), Q(1, 10))
        self.assertAlmostEqual(float(zero_tilt), .1)

    def test_refines_only_the_requested_entries(self):
        records = {1: (1000, 10, 12, 1000), 2: (1000, 9, 8, 1000)}
        before = arb_mat([[arb(1)/2 for _ in range(11)] for _ in range(11)])
        for source in (M, L48, L56):
            before[source, C] = 0
        spectrum = {48: 5, 56: 7, 64: 11, 72: 13, 80: 17}
        after = refine(before, (1, 2), records, '.052', spectrum=spectrum)
        allowed = {(C, C)} | {(U+i, C) for i in range(5)}
        for i, j in product(range(11), repeat=2):
            if (i, j) in allowed:
                self.assertLess(after[i, j], before[i, j])
            else:
                self.assertEqual(after[i, j], before[i, j])
        only_density = refine(before, (1, 2), records, '.052')
        self.assertTrue(all(only_density[U+i, C] == before[U+i, C] for i in range(5)))
        self.assertEqual(refine(before, (), records, '.052'), before)
        self.assertEqual(refine(before, (1,), records, '.052'), before)
        self.assertEqual(before[C, C], arb(1)/2)

    def test_invalid_inputs_and_coupled_columns(self):
        records = {1: (10, 1, 2, 10)}
        for shape, tilt in (((1,), 0), ((1, 1), -1), ((1, 1), .05), ((1, 0), 0)):
            with self.assertRaises(ValueError):
                bound(shape, records, tilt)
        with self.assertRaises(ValueError):
            exact_bound((1, 1), {1: Q(1, 2)}, Q(0))
        with self.assertRaises(ValueError):
            density_caps({1: (1, 11, 0, 10)})
        with self.assertRaises(ValueError):
            bound((1, 2), records, 0)
        for source in (M, L48, L56):
            matrix = arb_mat(11, 11)
            matrix[source, C] = 1
            with self.assertRaises(ValueError):
                refine(matrix, (1, 1), records, 0)
        for rounds in (0, -1, 2.0, True):
            with self.assertRaises(ValueError):
                refine(arb_mat(11, 11), (1, 1), records, 0, rounds=rounds)


if __name__ == '__main__':
    unittest.main()

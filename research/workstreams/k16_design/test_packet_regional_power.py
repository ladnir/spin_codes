"""Exact and outward toy checks; no state census or large proof computation."""
from itertools import combinations
from math import comb
import unittest

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat
import packet_q1 as q1
import packet_regional_power as power


def identity(value):
    return value


def rational_operators():
    return [fmpq_mat([[fmpq(2, 3), fmpq(1, 7)], [0, fmpq(3, 5)]]),
            fmpq_mat([[fmpq(1, 11), 0], [fmpq(2, 7), fmpq(1, 3)]]),
            fmpq_mat([[fmpq(1, 5), fmpq(1, 2)], [fmpq(1, 13), 0]])]


class RegionalPowerTests(unittest.TestCase):
    def test_exact_noncommuting_coefficients_match_linear_placement(self):
        ops = rational_operators()
        self.assertNotEqual(ops[0] * ops[1], ops[1] * ops[0])
        for epochs in (1, 2, 3, 5, 7, 8, 13):
            for degree in sorted({0, 1, min(5, 2*epochs), 2*epochs}):
                expected = q1.placement(ops, epochs, 2, fmpq_mat, identity,
                                        maximum_groups=degree)
                actual = power.placement_power(ops, epochs, 2, fmpq_mat, identity,
                                               maximum_groups=degree)
                self.assertEqual(actual, expected)

    def test_direct_slot_subset_average(self):
        ops = rational_operators()
        epochs, windows = 3, 2
        actual = power.placement_power(ops, epochs, windows, fmpq_mat, identity,
                                       maximum_groups=epochs*windows)
        for q in range(epochs*windows + 1):
            total = fmpq_mat(2, 2)
            for slots in combinations(range(epochs*windows), q):
                product = fmpq_mat([[1, 0], [0, 1]])
                for epoch in range(epochs):
                    count = sum(slot // windows == epoch for slot in slots)
                    product *= ops[count]
                total += product
            self.assertEqual(actual[q], total / comb(epochs*windows, q))

    def test_truncation_retains_every_requested_coefficient(self):
        ops = rational_operators()
        full = power.placement_power(ops, 7, 2, fmpq_mat, identity, maximum_groups=14)
        for degree in (0, 1, 2, 4, 9, 13):
            truncated = power.placement_power(ops, 7, 2, fmpq_mat, identity,
                                              maximum_groups=degree)
            self.assertEqual(truncated, full[:degree+1])
        # Degree one does not need the unretained local degree-two operator.
        short = power.placement_power(ops[:2], 7, 2, fmpq_mat, identity, maximum_groups=1)
        self.assertEqual(short, full[:2])

    def test_convolution_order_and_rounding_after_full_coefficient_sum(self):
        a, b, c = rational_operators()
        completed = []
        def record(value):
            completed.append(value)
            return value
        result = power._product([a, b], [c, a], 2, fmpq_mat, record)
        self.assertEqual(result, [a*c, a*a+b*c, b*a])
        self.assertEqual(completed, result)
        self.assertEqual(len(completed), 3)

    def test_arb_intervals_contain_exact_at_192_and_256_bits(self):
        exact_ops = rational_operators()
        previous = ctx.prec
        try:
            for precision in (192, 256):
                ctx.prec = precision
                intervals = [arb_mat(op) for op in exact_ops]
                for epochs in (3, 8, 13):
                    degree = min(9, 2*epochs)
                    exact = power.placement_power(exact_ops, epochs, 2, fmpq_mat,
                                                   identity, maximum_groups=degree)
                    balls = power.placement_power(intervals, epochs, 2, arb_mat,
                                                   identity, maximum_groups=degree)
                    for truth, enclosure in zip(exact, balls):
                        for i in range(2):
                            for j in range(2):
                                # contains(fmpq) first encloses the rational as an Arb
                                # ball; compare exact dyadic endpoints instead.
                                self.assertLessEqual(enclosure[i, j].lower().fmpq(), truth[i, j])
                                self.assertGreaterEqual(enclosure[i, j].upper().fmpq(), truth[i, j])
        finally:
            ctx.prec = previous

    def test_outward_endpoints_dominate_exact_at_192_and_256_bits(self):
        exact_ops = rational_operators()
        previous = ctx.prec
        try:
            for precision in (192, 256):
                ctx.prec = precision
                envelopes = [power.rounded(arb_mat(op)) for op in exact_ops]
                for epochs in (1, 3, 8, 13):
                    degree = min(9, 2*epochs)
                    exact = q1.placement(exact_ops, epochs, 2, fmpq_mat, identity,
                                          maximum_groups=degree)
                    result = power.placement_power(envelopes, epochs, 2,
                                                     maximum_groups=degree)
                    prior = q1.placement(envelopes, epochs, 2, arb_mat, power.rounded,
                                           maximum_groups=degree)
                    for truth, endpoint, old in zip(exact, result, prior):
                        for i in range(2):
                            for j in range(2):
                                value = endpoint[i, j]
                                self.assertTrue(value.is_exact())
                                self.assertGreaterEqual(value.fmpq(), truth[i, j])
                                self.assertGreaterEqual(old[i, j].fmpq(), truth[i, j])
                                # Both are valid envelopes; neither must dominate the other.
                                self.assertLess(abs(value-old[i, j]), arb(2)**(-precision//2))
        finally:
            ctx.prec = previous

    def test_local_inputs_remain_unchanged_and_default_degree_matches(self):
        ops = rational_operators()
        copies = [fmpq_mat(op) for op in ops]
        actual = power.placement_power(ops, 5, 2, fmpq_mat, identity)
        expected = q1.placement(ops, 5, 2, fmpq_mat, identity)
        self.assertEqual(actual, expected)
        self.assertEqual(ops, copies)

    def test_invalid_geometries_missing_operators_and_negative_entries(self):
        ops = rational_operators()
        for options in (dict(epochs=0), dict(epochs=True), dict(windows=0),
                        dict(maximum_groups=-1), dict(maximum_groups=129),
                        dict(maximum_groups=True)):
            with self.assertRaises(ValueError):
                settings = dict(windows=2, matrix=fmpq_mat, rounding=identity)
                settings.update(options)
                power.placement_power(ops, **settings)
        with self.assertRaises(ValueError):
            power.placement_power(ops[:2], 3, 2, fmpq_mat, identity, maximum_groups=2)
        with self.assertRaises(ValueError):
            power.placement_power([fmpq_mat([[-1]])], 3, 2, fmpq_mat, identity)
        with self.assertRaises(ValueError):
            power.placement_power([arb_mat([[arb('inf')]])])
        with self.assertRaises(ValueError):
            power.placement_power([fmpq_mat(2, 3)], matrix=fmpq_mat, rounding=identity)


if __name__ == '__main__':
    unittest.main()

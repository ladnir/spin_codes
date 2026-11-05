"""Finite-field checks for the XOR-column theory probe; no benchmarks."""
from fractions import Fraction
from itertools import product
from math import comb
import unittest

from hd_xor_column_probe import (column_support_counts, column_weight_counts, membership_mass,
                                 mds_prescribed_support, two_active_row_weight_enumerator,
                                 uniform_column_support_counts, xor_column)


def mul4(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return out


class XorColumnTests(unittest.TestCase):
    def test_involution(self):
        for values in product(range(4), repeat=4):
            self.assertEqual(xor_column(xor_column(values)), values)

    def test_exact_column_kernel(self):
        for q, group in ((2, 8), (4, 4), (8, 4)):
            for active in range(group+1):
                expected = [0]*(group+1)
                for values in product(range(1, q), repeat=active):
                    output = xor_column(values+(0,)*(group-active))
                    expected[sum(bool(x) for x in output)] += 1
                self.assertEqual(column_weight_counts(q, group, active), expected)

    def test_exact_support_kernel(self):
        for mask in range(16):
            active = [i for i in range(4) if mask >> i & 1]
            expected = [0]*16
            for nonzero in product(range(1, 4), repeat=len(active)):
                values = [0]*4
                for i, value in zip(active, nonzero):
                    values[i] = value
                output = xor_column(values)
                expected[sum(bool(value) << i for i, value in enumerate(output))] += 1
            self.assertEqual(column_support_counts(4, 4, mask), expected)

    def test_uniform_comparison_kernel(self):
        for mask in range(16):
            active = [i for i in range(4) if mask >> i & 1]
            expected = [0]*16
            for amplitudes in product(range(4), repeat=len(active)):
                values = [0]*4
                for i, value in zip(active, amplitudes):
                    values[i] = value
                output = xor_column(values)
                expected[sum(bool(value) << i for i, value in enumerate(output))] += 1
            self.assertEqual(uniform_column_support_counts(4, 4, mask), expected)
        q = 256
        counts = uniform_column_support_counts(q, 8, 3)
        self.assertEqual({i: count for i, count in enumerate(counts) if count},
                         {0: 1, 3: q-1, 253: q-1, 254: q-1, 255: (q-1)*(q-2)})
        full = uniform_column_support_counts(q, 8, 255)
        self.assertEqual(full, [(q-1)**mask.bit_count() for mask in range(256)])

    def test_mds_support_count(self):
        counts = {}
        for a, b in product(range(4), repeat=2):
            word = tuple(a ^ mul4(b, x) for x in range(4))
            support = tuple(i for i, value in enumerate(word) if value)
            counts[support] = counts.get(support, 0)+1
        for support_bits in range(16):
            support = tuple(i for i in range(4) if support_bits >> i & 1)
            self.assertEqual(counts.get(support, 0), mds_prescribed_support(4, 4, 2, len(support)))

    def test_membership_direct_setup_enumeration(self):
        # Two RS[2,1] rows over GF4, every setup and nonzero message.
        actual = [0]*256
        setups = 0
        for diagonal in product(range(1, 4), repeat=4):
            setups += 1
            for message in product(range(4), repeat=2):
                if not any(message):
                    continue
                columns = tuple(xor_column(tuple(mul4(diagonal[2*j+i], message[i])
                                                       for i in range(2))) for j in range(2))
                packed = sum(value << (2*(2*j+i)) for j, column in enumerate(columns)
                             for i, value in enumerate(column))
                actual[packed] += 1
        for packed, count in enumerate(actual):
            columns = tuple(tuple((packed >> (2*(2*j+i))) & 3 for i in range(2))
                            for j in range(2))
            self.assertEqual(Fraction(count, setups), membership_mass(4, 2, 1, columns))

    def test_four_row_membership_direct(self):
        # RS[2,1] means both pre-diagonal symbols in each row must be equal.
        # Enumerate every setup for selected targets, with a direct inverse.
        targets = (1, 85, 255, 1023, 4369, 21845, 65535, 12345)
        inverses = (0, 1, 3, 2)
        for packed in targets:
            columns = tuple(tuple((packed >> (2*(4*j+i))) & 3 for i in range(4))
                            for j in range(2))
            unmixed = tuple(xor_column(column) for column in columns)
            accepted = 0
            for diagonal in product(range(1, 4), repeat=8):
                first = tuple(mul4(inverses[diagonal[i]], unmixed[0][i]) for i in range(4))
                second = tuple(mul4(inverses[diagonal[4+i]], unmixed[1][i]) for i in range(4))
                accepted += first == second and any(first)
            self.assertEqual(Fraction(accepted, 3**8), membership_mass(4, 2, 1, columns))

    def test_two_row_total_mass(self):
        for q, n, k, group in ((4, 4, 2, 4), (256, 32, 16, 8)):
            actual = sum(two_active_row_weight_enumerator(q, n, k, group))
            expected = (q**k-1)**2
            self.assertAlmostEqual(actual/expected, 1.0, places=13)


if __name__ == '__main__':
    unittest.main()

import unittest

from gf8_packet_maps import matvec, rank
from hd_bit_syndrome_outer import matrix
from hd_bit_syndrome_singleton_cosets import (
    EXCEPTION_CAPS, EXCEPTION_WEIGHTS, GENERIC_CAPS, GENERIC_WEIGHTS,
    UNIVERSAL_CAPS, audit_all, counts_for, singleton_data, support_caps,
    weight_counts,
)


class SingletonCosetTests(unittest.TestCase):
    def test_all_1920_exact_audit(self):
        data = audit_all()
        self.assertEqual(data['exceptional_fibers'], 32)
        self.assertEqual(data['generic_fibers'], 1888)
        self.assertEqual(data['syndrome_binary_ranks'], {5: 224, 4: 16})
        self.assertEqual(data['universal_per_support_caps'], UNIVERSAL_CAPS)
        self.assertFalse(data['whole_code_certificate'])

    def test_representative_and_actual_saved_worst(self):
        for a, j in [(tuple(range(7)) + (8,), 7), (tuple(range(7, 15)), 15)]:
            counts = counts_for(a, j)
            self.assertEqual(weight_counts(counts), GENERIC_WEIGHTS)
            self.assertEqual(support_caps(counts), GENERIC_CAPS)
            self.assertEqual(sum(counts.values()), 256 ** 3)
            self.assertEqual(min(mask.bit_count() for mask in counts), 5)

    def test_rank4_exception_is_retained(self):
        a, j = (0, 1, 4, 6, 8, 10, 12, 14), 2
        counts = counts_for(a, j)
        self.assertEqual(singleton_data(a, j)['binary_rank'], 4)
        self.assertEqual(weight_counts(counts), EXCEPTION_WEIGHTS)
        self.assertEqual(support_caps(counts), EXCEPTION_CAPS)
        self.assertEqual(counts[sum(1 << i for i in (0, 1, 4, 10))], 1)

    def test_independent_full_field_matrix_constraints(self):
        h = matrix()
        for a, j in [(tuple(range(7)) + (8,), 7),
                     ((0, 1, 4, 6, 8, 10, 12, 14), 2)]:
            counts = counts_for(a, j)
            outside = tuple(i for i in range(16) if i not in a)
            # Independent constraints on x: H_outside,A x=e_j and
            # H_A\S,A x=0. No V/B/syndrome formula is used here.
            for local in range(256):
                support = sum(1 << a[k] for k in range(8) if local >> k & 1)
                rows = outside + tuple(i for i in a if not (support >> i & 1))
                coefficients = [[h[i][k] for k in a] for i in rows]
                augmented = [row + [int(i == j)] for row, i in zip(coefficients, rows)]
                r = rank(coefficients)
                expected = 256 ** (8 - r) if rank(augmented) == r else 0
                self.assertEqual(sum(n for mask, n in counts.items() if not (mask & ~support)), expected)

    def test_explicit_field_matrix_witness_and_scaling(self):
        x = (81, 76, 125, 1, 42, 1, 1, 0, 249, 0, 0, 0, 0, 0, 0, 0)
        y = (80, 77, 124, 0, 43, 0, 0, 1, 249, 0, 0, 0, 0, 0, 0, 0)
        self.assertEqual(matvec(matrix(), x), y)
        self.assertEqual(matvec(matrix(), y), x)
        self.assertEqual(sum(bool(v) for v in y), 6)
        # The support coefficients are for one outside value, not all255.
        self.assertEqual(sum(counts_for(tuple(range(7)) + (8,), 7).values()) * 255,
                         255 * 256 ** 3)

    def test_impossible_and_invalid_inputs(self):
        self.assertEqual(counts_for(tuple(range(8)), 8), {})
        self.assertEqual(counts_for(tuple(range(7)) + (8,), 9), {})
        with self.assertRaises(ValueError):
            counts_for(tuple(range(8)), 7)
        with self.assertRaises(ValueError):
            counts_for((0, 1, 2), 4)


if __name__ == '__main__':
    unittest.main()

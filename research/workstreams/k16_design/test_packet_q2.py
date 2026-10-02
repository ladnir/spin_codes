"""Exact small ordered-placement checks for the occupancy-two diagnostic."""
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from itertools import combinations
from math import comb
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, fmpq, fmpq_mat
import packet_q2 as screen


class PacketQ2Tests(unittest.TestCase):
    def test_ordered_coefficients_match_all_pairs_of_support_subsets(self):
        matrices = [fmpq_mat([[1, 1], [0, 1]]),
                    fmpq_mat([[1, 0], [1, 1]]),
                    fmpq_mat([[2, 1], [1, 1]])]
        self.assertNotEqual(matrices[0] * matrices[1], matrices[1] * matrices[0])
        for regions in range(5):
            actual = screen.pair_support_moments(matrices, regions,
                matrix=fmpq_mat, scalar=fmpq, rounding=lambda x: x)
            for v in range(regions + 1):
                for u in range(v + 1):
                    total = fmpq(0)
                    for first in combinations(range(regions), u):
                        for second in combinations(range(regions), v):
                            row = fmpq_mat([[1, 0]])
                            for r in range(regions):
                                row *= matrices[int(r in first) + int(r in second)]
                            total += row[0, 0] + row[0, 1]
                    self.assertEqual(actual[v][u], total / (comb(regions, u) * comb(regions, v)))

    def test_no_reset_at_region_boundary(self):
        # One occupied region births a state; the empty region keeps and expands
        # it. Resetting to zero between regions would instead give moment two.
        matrices = [fmpq_mat([[1, 0], [0, 3]]),
                    fmpq_mat([[0, 2], [0, 1]]), fmpq_mat([[1, 0], [0, 1]])]
        actual = screen.pair_support_moments(matrices, 2,
            matrix=fmpq_mat, scalar=fmpq, rounding=lambda x: x)
        self.assertEqual(actual[1][0], 4)

    def test_all_ones_scalar_normalization(self):
        matrices = [fmpq_mat([[1]]) for _ in range(3)]
        actual = screen.pair_support_moments(matrices, 8,
            matrix=fmpq_mat, scalar=fmpq, rounding=lambda x: x)
        self.assertTrue(all(value == 1 for column in actual for value in column))

    def test_outward_arithmetic_encloses_exact_rational_answer(self):
        matrices = [fmpq_mat([[1, 2], [3, 4]]) / 7,
                    fmpq_mat([[3, 1], [2, 4]]) / 11,
                    fmpq_mat([[5, 2], [1, 3]]) / 13]
        exact = screen.pair_support_moments(matrices, 4,
            matrix=fmpq_mat, scalar=fmpq, rounding=lambda x: x)
        outward = screen.pair_support_moments([arb_mat(m) for m in matrices], 4)
        for v, column in enumerate(exact):
            for u, value in enumerate(column):
                mantissa, exponent = outward[v][u].man_exp()
                bound = Q(int(mantissa)) * Q(2) ** int(exponent)
                self.assertGreaterEqual(bound, Q(int(value.numerator), int(value.denominator)))

    def test_regional_then_bivariate_matches_all_slot_pairs(self):
        # Two epochs of two slots: each group occupies one distinct slot when
        # present. This checks the regional collision operator's normalization.
        local = [fmpq_mat([[1, 1], [0, 1]]), fmpq_mat([[1, 0], [1, 1]]),
                 fmpq_mat([[2, 1], [1, 1]])]
        regional = screen.q1.placement(local, epochs=2, windows=2,
            matrix=fmpq_mat, rounding=lambda x: x, maximum_groups=2)
        expected = fmpq_mat(2, 2)
        for first in range(4):
            for second in range(4):
                if first == second:
                    continue
                expected += local[int(first < 2) + int(second < 2)] * local[int(first >= 2) + int(second >= 2)]
        self.assertEqual(regional[2], expected / 12)
        moment = screen.pair_support_moments(regional, 1,
            matrix=fmpq_mat, scalar=fmpq, rounding=lambda x: x)[1][1]
        self.assertEqual(moment, (expected[0, 0] + expected[0, 1]) / 12)

    def test_identity_kernel_counts_ordered_messages_in_unordered_groups(self):
        local = [arb_mat([[1]]) for _ in range(3)]
        with patch.object(screen.q1.kernel_t64, 'local_operators', return_value=local), redirect_stdout(StringIO()):
            record = screen.evaluate_q2((0, 5, 10), geometry=screen.q1.Geometry(32, 2, 4),
                tilts=['.01'], data={}, map_record={})
        mantissa, exponent = record['q2_upper']
        self.assertEqual(Q(mantissa) * Q(2) ** exponent, comb(32, 2) * 15**2)
        self.assertFalse(record['whole_code_certificate'])
        self.assertEqual(record['occupancy_covered'], [2])

    def test_cdf_and_bad_totals_rejected(self):
        for counts in ((0, 5, 15), (0, -1, 16), (1, 5, 9)):
            with self.assertRaises(ValueError):
                screen.evaluate_q2(counts, geometry=screen.q1.Geometry(32, 2, 4), tilts=['.01'])
        with self.assertRaises(ValueError):
            screen.pair_support_moments([arb_mat([[1]])] * 2, 1)


if __name__ == '__main__':
    unittest.main()

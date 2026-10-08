"""Exact small-field tests of the independently derived whole-code bounds."""

from fractions import Fraction as F
from itertools import combinations
from math import comb
import unittest

from batched_transfer import active_output_moment, transition
import subfield_bound_audit as audit


def gf4_one_step_regions(z):
    # L=t=4, so the regional matrix is the average over its h-lane masks.
    regions = []
    for h in range(5):
        matrix = audit.ZERO
        for support in combinations(range(4), h):
            matrix = audit.matadd(matrix, transition(support, 4, 4, 0x7, z))
        regions.append(audit.matscale(matrix, F(1, comb(4, h))))
    return regions


class SubfieldBoundAudit(unittest.TestCase):
    def test_mds_support_envelope(self):
        for q, length, dimension in ((4, 4, 2), (256, 16, 8)):
            self.assertEqual(audit.nonzero_word_enumerator(
                length, dimension, q, F(1)), q ** dimension - 1)
            for weight in range(length + 1):
                actual = audit.fixed_support_multiplicity(length, dimension, weight, q)
                exponent = weight - (length - dimension)
                envelope = F(q ** exponent) if exponent >= 0 else F(1, q ** -exponent)
                self.assertGreaterEqual(actual, 0)
                self.assertLessEqual(actual, envelope)

    def test_coupled_bound_every_gf4_active_row_count(self):
        for z in (F(1, 4), F(3, 4), F(1)):
            regional = gf4_one_step_regions(z)
            for active_rows in range(1, 5):
                literal = audit.literal_fixed_rows_moment(regional, active_rows, 4, 2, 4)
                envelope = audit.coupled_fixed_rows_bound(regional, active_rows, 4, 2, 4)
                self.assertLessEqual(literal, envelope)
                if z == 1:
                    self.assertEqual(literal, 15 ** active_rows)

    def test_scalar_potential_envelope(self):
        regional = gf4_one_step_regions(F(1, 4))
        eta = F(1, 2)
        for active_rows in range(1, 5):
            rho = next(F(i, 100) for i in range(1, 201)
                       if all(audit.potential_norm(regional[h], eta) <= F(i, 100) ** h
                              for h in range(1, active_rows + 1)))
            literal = audit.literal_fixed_rows_moment(regional, active_rows, 4, 2, 4)
            envelope = audit.scalar_enumerator_bound(
                regional, active_rows, 4, 2, 4, eta, rho)
            self.assertLessEqual(literal, envelope)

    def test_universal_bulk_row_sum_bound(self):
        for q, width in ((4, 4), (8, 4), (8, 8), (256, 16)):
            for z in (F(0), F(1, 1000), F(1, 4), F(3, 4), F(999, 1000), F(1)):
                f = z + (1 - z) / (q - 1)
                for h in range(width + 1):
                    g = active_output_moment(h, width, q, z)
                    self.assertLessEqual(g, f ** width)
                    self.assertLessEqual(g, f ** h)
                    self.assertLessEqual(z ** h, f ** h)

    def test_exact_tail_check_matches_fraction_expression(self):
        q, rows, first, z = 4, 4, 3, F(1, 10)
        result = audit.exact_bulk_tail_check(first=first, z=z, margin_bits=0,
                                            k=8, length=4, dimension=2, q=q,
                                            distance=F(1, 16))
        p = audit.nonzero_word_enumerator(4, 2, q, z + (1 - z) / (q - 1))
        ratio = F(rows - first, first + 1) * p
        upper = F(comb(rows, first), q - 1) * p ** first / (1 - ratio)
        self.assertEqual(result['rigorous'], upper < 1)

    def test_window_obstruction_exact_integer_witness(self):
        for state, active_rows, surviving in ((2, 32, 226), (4, 64, 452), (8, 128, 904)):
            result = audit.exact_window_lower_check(active_rows=active_rows, state=state)
            self.assertTrue(result['rigorous'])
            self.assertTrue(result['expectation_not_failure_probability'])
            self.assertEqual(result['surviving_dimension'], surviving)
            self.assertLess(result['maximum_output_weight'], 13108)


if __name__ == '__main__':
    unittest.main()

"""Exact order-statistic laws and directed current-geometry proof receipts."""

from itertools import combinations
from math import comb
import unittest

import hd_blockwise_interval_spans as proof


class IntervalSpanTests(unittest.TestCase):
    def test_exact_order_statistic_counts(self):
        for rows in range(1, 9):
            for active in range(1, rows+1):
                anchor, span = proof.order_statistic_counts(rows, active)
                observed_anchor, observed_span = [0]*(rows+1), [0]*(rows+1)
                for positions in combinations(range(rows), active):
                    observed_anchor[positions[-1]+1] += 1
                    observed_span[positions[-1]-positions[0]+1] += 1
                self.assertEqual(anchor, observed_anchor)
                self.assertEqual(span, observed_span)
                self.assertEqual(sum(anchor), comb(rows, active))
                self.assertEqual(sum(span), comb(rows, active))

    def test_new_geometry_actual_moment(self):
        receipt = proof.new_geometry_moment()
        self.assertEqual(receipt['dimension_lower'], 260)
        self.assertIn('740.753174823854', receipt['actual_projective_moment_log2_lower'])
        self.assertTrue(receipt['moment_exceeds_one'])
        self.assertTrue(receipt['substantial_failure_probability_not_established'])

    def test_new_geometry_complete_row_cap(self):
        receipt = proof.complete_row_window_cap()
        self.assertEqual(receipt['grouped_terms'], 31469)
        self.assertEqual(receipt['active_sizes'], list(range(3, 10)))
        self.assertIn('-528.39821182', receipt['capped_family_log2_upper'])
        self.assertTrue(receipt['below_2_to_minus_40'])
        self.assertFalse(receipt['whole_code_certificate'])
        self.assertTrue(receipt['partial_row_shortening_not_covered'])


if __name__ == '__main__':
    unittest.main()

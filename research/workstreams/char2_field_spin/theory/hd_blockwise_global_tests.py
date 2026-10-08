"""Exact small structural tests and directed K20 global gate receipts."""

from fractions import Fraction as F
from itertools import product
from math import comb
import unittest

import hd_blockwise_global as proof
from hd_blockwise_inner import add, allowed_pairs, encode, matvec, transpose


class BlockwiseGlobalTests(unittest.TestCase):
    def test_CA_zero_involution_and_packet_syndromes(self):
        # Evaluation of degree<2 polynomials over GF4; A^T A=0.
        A = tuple((1, x) for x in range(4))
        C = transpose(A)
        for column in transpose(A):
            self.assertEqual(matvec(C, column, 4), (0, 0))
        pairs = list(allowed_pairs(4, 2))
        schedule = [(pairs[17],), (pairs[53], pairs[71]), (pairs[121],)]
        for position in range(12):
            for scalar in range(1, 4):
                flat = [0]*12
                flat[position] = scalar
                x = tuple(tuple(flat[4*i:4*i+4]) for i in range(3))
                y = encode(x, A, C, schedule, 4)
                self.assertEqual(encode(y, A, C, schedule, 4), x)
                self.assertEqual(tuple(matvec(C, packet, 4) for packet in x),
                                 tuple(matvec(C, packet, 4) for packet in y))

    def test_feedback_zero_packets_are_fixed(self):
        A = tuple((1, x) for x in range(4))
        C = transpose(A)
        kernel = [x for x in product(range(4), repeat=4) if not any(matvec(C, x, 4))]
        self.assertEqual(len(kernel), 4**2)
        pairs = list(allowed_pairs(4, 2))
        schedule = [(pairs[17],), (pairs[53],)]
        for left, right in product(kernel, repeat=2):
            x = (left, right)
            self.assertEqual(encode(x, A, C, schedule, 4), x)

    def test_repeated_outer_rows_do_not_create_a_low_distance_claim(self):
        # GF4 [4,2] GRS has minimum3, the toy analogue of65/128.
        A = tuple((1, x) for x in range(4))
        pairs = list(allowed_pairs(4, 2))
        weights = []
        for message in product(range(4), repeat=2):
            if not any(message):
                continue
            word = matvec(A, message, 4)
            x = tuple((value,)*4 for value in word)
            y = encode(x, A, transpose(A), [(pairs[i],) for i in (2, 19, 37, 95)], 4)
            self.assertEqual(y, x)
            weights.append(sum(value != 0 for packet in y for value in packet))
        self.assertEqual(min(weights), 12)
        self.assertEqual(F(65, 128), F(1, 2)+F(1, 128))

    def test_sampling_without_replacement_containment_bound(self):
        for rows in range(2, 9):
            for length in range(rows+1):
                for active in range(rows+1):
                    actual = F(comb(length, active), comb(rows, active)) if active <= length else F(0)
                    self.assertLessEqual(actual, F(length, rows)**active)

    def test_connected_pattern_count_dominates_small_windows(self):
        # Enumerate all one-interval-per-region families, including boundary joins.
        # Count classification against the unrestricted-composition upper used here.
        from collections import Counter
        b, blocks, state = 4, 4, 1
        rows, width, distance = blocks*2, 2, F(2, 5)
        intervals = [(start, end) for start in range(blocks) for end in range(start+1, blocks+1)]
        actual = Counter()
        budget = (distance.numerator*b*rows-1)//distance.denominator//width
        for windows in product(intervals, repeat=b):
            if sum(end-start for start, end in windows) > budget:
                continue
            component_sizes, size = [], 1
            for c in range(b-1):
                if windows[c][1] == blocks and windows[c+1][0] == 0:
                    size += 1
                else:
                    component_sizes.append(size)
                    size = 1
            component_sizes.append(size)
            components = len(component_sizes)
            singletons = component_sizes.count(1)
            actual[components, singletons, 'charged'] += 1
            if windows[-1][1] == blocks:
                terminal = 'singleton_tail' if component_sizes[-1] == 1 else 'long_tail'
                actual[components, singletons, terminal] += 1
        upper = {(t['components'], t['singletons'], t['terminal']): t['window_families']
                 for t in proof.connected_terms(b, rows, width, state, distance)}
        self.assertTrue(actual)
        for key, count in actual.items():
            self.assertLessEqual(count, upper.get(key, 0), key)

    def test_source_set_geometric_tail(self):
        for rows, minimum, ratio in ((16, 3, F(1, 20)), (12, 1, F(1, 128)), (8, 4, F(1, 8))):
            step = F(rows-minimum, minimum+1)*ratio
            total = sum(comb(rows, a)*ratio**a for a in range(minimum, rows+1))
            bound = comb(rows, minimum)*ratio**minimum/(1-step)
            self.assertLessEqual(total, bound)

    def test_actual_first_moment_obstruction(self):
        receipt = proof.moment_obstruction(39)
        self.assertEqual(receipt['message_dimension_lower'], 464)
        self.assertIn('780.179816185', receipt['actual_projective_moment_log2_lower'])
        self.assertTrue(receipt['moment_exceeds_one'])
        self.assertTrue(receipt['distance_failure_not_established'])

    def test_capped_connected_family_directed_gate(self):
        receipt = proof.capped_full_row_traps(16)
        self.assertEqual(receipt['grouped_terms'], 7779)
        self.assertIn('-1511.827110669', receipt['capped_family_log2_upper'])
        self.assertTrue(receipt['below_2_to_minus_40'])
        self.assertTrue(receipt['partial_row_shortening_not_covered'])
        self.assertFalse(receipt['whole_code_certificate'])

    def test_partial_row_naive_union_does_not_close(self):
        receipt = proof.fixed_window_capacity()
        self.assertIn('-612.205309228', receipt['fixed_profile_log2_upper'])
        self.assertFalse(receipt['naive_union_closes'])


if __name__ == '__main__':
    unittest.main()

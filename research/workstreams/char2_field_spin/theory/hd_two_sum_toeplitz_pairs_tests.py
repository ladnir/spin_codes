"""Tiny exact-law and positive contraction controls; no distance sweep."""
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest

from hd_grouped_toeplitz import binary_field, toy_counts
from hd_two_sum_closed_kernel import uniform_mask_count
import hd_two_sum_toeplitz_pairs as target


class ToeplitzPairTests(unittest.TestCase):
    def test_exhaustive_exact_vector_mixture(self):
        for q in (2, 4):
            counts, accepted = toy_counts(q, 2)
            self.assertEqual(accepted, q*q*(q-1))
            mul = binary_field(q)
            for word in range(1, q*q):
                x = word % q, word // q
                line = {mul[c][x[0]]+q*mul[c][x[1]] for c in range(1, q)}
                for y in range(q*q):
                    expected = F(0) if not y else F(1, q*q)+F(int(y in line), q*q*(q-1))
                    self.assertEqual(F(counts[word][y], accepted), expected)

    def test_exact_support_and_pointwise_iid_majorant(self):
        for q in (2, 4):
            counts, accepted = toy_counts(q, 2)
            exact = target.pair_kernel(q)
            majorant = target.pair_kernel(q, 'iid-majorant')
            cap = target.pair_kernel(q, 'uniform-cap')
            for word in range(q*q):
                incoming = int(word % q != 0)+2*int(word // q != 0)
                observed = [F(0)]*4
                for y in range(q*q):
                    outgoing = int(y % q != 0)+2*int(y // q != 0)
                    actual = F(counts[word][y], accepted)
                    observed[outgoing] += actual
                    self.assertLessEqual(actual, majorant[incoming][outgoing]/(q-1)**outgoing.bit_count())
                self.assertEqual(tuple(observed), exact[incoming])
                self.assertTrue(all(a <= b <= c for a, b, c in zip(exact[incoming], majorant[incoming], cap[incoming])))

    def test_lazy_double_is_not_iid(self):
        counts, accepted = toy_counts(4, 2)
        # x=(1,1); same-support output(1,1) lies on the exceptional line.
        exact_double_mass = target.pair_kernel(4)[3][3]
        self.assertGreater(F(counts[5][5], accepted), exact_double_mass/9)

    def test_scalar_weight_contraction_toy(self):
        q = 4
        counts, accepted = toy_counts(q, 2)
        for word in range(q*q):
            incoming = int(word % q != 0)+2*int(word // q != 0)
            left, right = F(2, 3), F(5, 7)
            actual = sum(F(counts[word][y], accepted)*left**int(y % q != 0)*right**int(y // q != 0)
                         for y in range(q*q))
            self.assertEqual(actual, target.pair_mark(incoming, left, right, q=q))

    def test_positive_column_contraction_matches_labelled_masks(self):
        group, active = 8, 0b00101101
        matching = target.xor_matching(3, group)
        marks = tuple(F(i+1, i+2) for i in range(group))
        for mode in target.MODES:
            direct = sum(uniform_mask_count(active, mask, group)
                         *target.support_factor(mask, matching, marks, mode=mode)
                         for mask in range(1 << group))
            self.assertEqual(direct, target.column_contraction(active, matching, marks, mode=mode))

    def test_pair_mask_partition_and_marked_normalization(self):
        group, active = 8, 0b00101101
        matching = target.xor_matching(5, group)
        for marker in (F(1), F(3, 2)):
            counts = [target.pair_mask_count(active, matching, mask, group=group, input_marker=marker)
                      for mask in range(1 << (group//2))]
            self.assertEqual(counts[0], 1)
            self.assertEqual(sum(counts), (1+255*marker)**active.bit_count())
            self.assertEqual(target.column_contraction(active, matching, (F(1),)*group, input_marker=marker),
                             (1+255*marker)**active.bit_count())

    def test_full_arity_uniform_law_preserved(self):
        group = 8
        marks = tuple(F(i+1, i+2) for i in range(group))
        expected = F(1)
        for mark in marks:
            expected *= 1+255*mark
        self.assertEqual(target.column_contraction(255, target.xor_matching(1, group), marks), expected)

    def test_pair_types_and_schedule_incidence(self):
        # All15 directions partition every unordered edge, followed by d=1.
        for active in (0, 1, 3, 15, 255, 0x5555, 65535):
            internal = [sum((active >> i & 1) and (active >> j & 1) for i, j in matching)
                        for matching in target.xor_schedule()]
            self.assertEqual(sum(internal), comb(active.bit_count(), 2)+internal[-1])
        # Exact labelled enumeration independently checks the positive pair DP.
        active, matching = 15, target.xor_matching(3, 8)
        counts = [0]*16
        for mask in range(256):
            pair_mask = sum(bool(mask & ((1 << i) | (1 << j))) << p for p, (i, j) in enumerate(matching))
            counts[pair_mask] += uniform_mask_count(active, mask, 8)
        self.assertEqual(counts, [target.pair_mask_count(active, matching, m, group=8) for m in range(16)])

    def test_invalid_scope_and_exact_zero(self):
        for q in (True, 3, 0, 512):
            with self.assertRaises(ValueError):
                target.pair_kernel(q)
        with self.assertRaises(ValueError):
            target.column_contraction(1, ((0, 1),)*4, (F(1),)*8)
        with self.assertRaises(ValueError):
            target.pair_mark(1, -1, 1)
        for mode in target.MODES:
            self.assertEqual(target.pair_mark(0, 0, 0, mode=mode), 1)
            self.assertEqual(target.pair_mark(3, 0, 0, mode=mode), 0)


if __name__ == '__main__':
    unittest.main()

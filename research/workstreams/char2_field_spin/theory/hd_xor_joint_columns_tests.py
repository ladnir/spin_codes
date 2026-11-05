"""Independent support-array counting, MDS envelopes, and directed replay."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import log2
import unittest

import numpy as np
from flint import arb, ctx

from hd_xor_column_probe import mds_prescribed_support
from hd_xor_column_potential import PAIRS, output_support_law, multiset_column_matrices, category_costs
from hd_xor_joint_columns import (SAVED_WITNESSES, categories,
    column_polynomial_coefficients, joint_count_arb, joint_log_count,
    receipt, replay, arb_upper_max, resource_envelope_arb,
    resource_envelope_float, occupancy_envelope_caps_float,
    joint_multiset_code_batches, replay_envelope, replay_envelope_family,
    replay_common_pair, occupancy_envelope_receipt)


def source_category(pattern, groups):
    sizes = [((pattern >> (2*group)) & 3).bit_count() for group in range(groups)]
    return sizes.count(1), sizes.count(2)


class JointColumnTests(unittest.TestCase):
    def test_column_multiplicities_exhaust_all_binary_patterns(self):
        for groups in range(1, 6):
            caps = {category: Fraction(1, 3**(category[0]+2*category[1]+category[0]*category[1]))
                    for category in categories(groups)}
            expected = [Fraction(0) for degree in range(2*groups+1)]
            for pattern in range(1 << (2*groups)):
                expected[pattern.bit_count()] += caps[source_category(pattern, groups)]
            self.assertEqual(column_polynomial_coefficients(caps, groups), expected)
            ones = column_polynomial_coefficients(dict.fromkeys(categories(groups), 1), groups)
            self.assertEqual(sum(ones), 4**groups)

    def test_information_set_envelope_for_all_prescribed_supports(self):
        for q, n, k in ((4, 3, 2), (8, 7, 3), (256, 32, 16)):
            for weight in range(n-k+1, n+1):
                count = mds_prescribed_support(q, n, k, weight)
                self.assertLessEqual(count, (q-1)**(weight-(n-k)))
            self.assertEqual(mds_prescribed_support(q, n, k, n-k+1), q-1)

    def test_positive_tail_equals_exhaustive_support_array_relaxation(self):
        groups, n, k, q = 2, 3, 2, 4
        caps = {category: Fraction(1, 3**category[0]*5**category[1]*7**(category[0]*category[1]))
                for category in categories(groups)}
        expected = Fraction(0)
        admitted_with_a_small_row = False
        for columns in product(range(16), repeat=3):
            total = sum(column.bit_count() for column in columns)
            if total < 8:
                continue
            value = Fraction(3**(total-4))
            for column in columns:
                value *= caps[source_category(column, groups)]
            expected += value
            row_weights = [sum(column >> row & 1 for column in columns) for row in range(4)]
            admitted_with_a_small_row |= min(row_weights) < 2
        # The extension genuinely relaxes individual floors; it does not
        # assert an RS prescribed-support formula for these added arrays.
        self.assertTrue(admitted_with_a_small_row)
        with ctx.workprec(256):
            directed_caps = {key: arb(value.numerator)/value.denominator for key, value in caps.items()}
            actual = joint_count_arb(directed_caps, groups, n, k, q)
            exact = arb(expected.numerator)/expected.denominator
            self.assertTrue(actual.overlaps(exact))
            floating = joint_log_count({key: log2(value) for key, value in caps.items()}, groups, n, k, q)
            self.assertAlmostEqual(floating, log2(expected), places=12)

    def test_gf4_codeword_union_is_below_relaxed_tail(self):
        def multiply(a, b):
            result = 0
            while b:
                if b & 1:
                    result ^= a
                b >>= 1
                a <<= 1
                if a & 4:
                    a ^= 7
            return result
        supports = Counter()
        for a, b in product(range(4), repeat=2):
            if a or b:
                word = [a ^ multiply(b, x) for x in range(3)]
                supports[sum(bool(value) << column for column, value in enumerate(word))] += 1
        self.assertEqual(sum(supports.values()), 15)
        self.assertTrue(all(support.bit_count() >= 2 for support in supports))
        caps = {category: Fraction(1, 3**category[0]*5**category[1]*7**(category[0]*category[1]))
                for category in categories(2)}
        exact = Fraction(0)
        for rows in product(supports, repeat=4):
            contribution = Fraction(1)
            for row in rows:
                contribution *= supports[row]
            for column in range(3):
                pattern = sum(((row >> column) & 1) << index for index, row in enumerate(rows))
                contribution *= caps[source_category(pattern, 2)]
            exact += contribution
        with ctx.workprec(256):
            relaxed = joint_count_arb({key: arb(value.numerator)/value.denominator for key, value in caps.items()}, 2, 3, 2, 4)
            self.assertTrue(arb(exact.numerator)/exact.denominator < relaxed)

    def test_directed_r5_saved_receipt(self):
        with ctx.workprec(384):
            result = replay(5, 384)
            self.assertTrue(arb(result['failure_log2_upper']) < -129)
            self.assertEqual(result['checked_support_cases'], 749398)
            self.assertEqual(result['source_row_pair_choices'], 28**5)
            self.assertEqual(result['total_support_floor'], 170)
            self.assertEqual(result['polynomial_degree'], 320)
            self.assertEqual(result['projective_normalizations'], 1)
            self.assertTrue(result['individual_support_floors_relaxed'])
            self.assertFalse(result['region_permutation_reduction'])
            self.assertFalse(result['whole_code_certificate'])

    def test_invalid_empty_cap_and_false_nonempty_cap_are_rejected(self):
        z, potential, saved = SAVED_WITNESSES[5]
        caps = {**saved, (0, 0): '-1'}
        with self.assertRaisesRegex(ValueError, 'Invalid marks'):
            receipt(5, z, potential, caps)
        caps = {**saved, (0, 1): '-200'}
        with self.assertRaisesRegex(ArithmeticError, 'joint category failed'):
            receipt(5, z, potential, caps)

    def test_directed_max_uses_upper_endpoints_of_overlapping_balls(self):
        first, second = arb('1 +/- .25'), arb('1.1 +/- .2')
        result = arb_upper_max(first, second)
        self.assertTrue(result.is_exact())
        self.assertTrue(result >= first.upper())
        self.assertTrue(result >= second.upper())

    def test_resource_envelope_bounds_every_noncommuting_word(self):
        matrices = np.array([[[1., 0.], [0., .125]], [[0., .5], [.0625, .25]],
                             [[.125, .375], [.03125, .5]]])
        directed = [tuple(arb(float(value)) for value in matrix.flatten()) for matrix in matrices]
        potential = np.array([1., .25])
        for lower, upper in ((0, 2), (1, 2), (2, 2)):
            with np.errstate(divide='ignore'):
                floating = resource_envelope_float(np.log2(matrices), -2., lower, upper, 5)
            bounded = resource_envelope_arb(directed, arb(.25), lower, upper, 5)
            for counts in product(range(lower, upper+1), repeat=5):
                action = potential.copy()
                for count in reversed(counts):
                    action = matrices[count] @ action
                index = sum(counts)-5*lower
                self.assertTrue(all(arb(float(action[row])) <= bounded[index][row] for row in range(2)))
            np.testing.assert_allclose(np.exp2(floating),
                                       np.array([[float(value) for value in vector] for vector in bounded]), rtol=1e-14)

    def test_occupancy_constraints_hold_for_every_two_group_branch(self):
        sources = (0, *(1 << index for index in range(8)), *PAIRS)
        for first, second in product(sources, repeat=2):
            single = int(first.bit_count() == 1)+int(second.bit_count() == 1)
            double = int(first.bit_count() == 2)+int(second.bit_count() == 2)
            for (left, _), (right, _) in product(output_support_law(first), output_support_law(second)):
                cancelled = int(first.bit_count() == 2 and left.bit_count() == 2)
                cancelled += int(second.bit_count() == 2 and right.bit_count() == 2)
                occupancies = [((left >> index) & 1)+((right >> index) & 1) for index in range(8)]
                self.assertGreaterEqual(min(occupancies), double-cancelled)
                self.assertLessEqual(max(occupancies), single+double)
                self.assertEqual(sum(occupancies), 7*single+8*double-6*cancelled)

    def test_resource_caps_dominate_exact_three_group_caps(self):
        z, potential = .99859, -62.40625
        exact, _ = category_costs(multiset_column_matrices(z, 3), potential)
        bounded = occupancy_envelope_caps_float(3, z, potential)
        for category in exact:
            self.assertGreaterEqual(bounded[category]+1e-11, exact[category])

    def test_six_group_generator_keeps_all_repeated_edge_branches(self):
        codes = next(joint_multiset_code_batches(6, 0, 6, batch_size=1))[0]
        self.assertEqual(len(codes), 64)
        for branch, code in enumerate(codes):
            cancelled = branch.bit_count()
            occupancies = [int(code)//7**index % 7 for index in range(8)]
            self.assertEqual(occupancies, [6, 6]+[6-cancelled]*6)
        self.assertEqual(next(joint_multiset_code_batches(6, 0, 0)).shape, (1, 1))

    def test_saved_six_group_resource_receipt(self):
        with ctx.workprec(384):
            result = replay_envelope(6, 384)
            self.assertTrue(arb(result['failure_log2_upper']) < -127)
            self.assertEqual(result['source_multisets_enumerated'], 0)
            self.assertEqual(result['covered_source_multisets'], 5245786)
            self.assertEqual(result['resource_tables'], 28)
            self.assertEqual(result['resource_state_vectors'], 476)
            self.assertEqual(result['total_support_floor'], 204)
            self.assertEqual(result['polynomial_degree'], 384)
            self.assertFalse(result['whole_code_certificate'])

    def test_saved_resource_family_covers_six_through_sixteen(self):
        with ctx.workprec(384):
            result = replay_envelope_family(384)
            self.assertEqual([row['active_groups'] for row in result['classes']], list(range(6, 17)))
            self.assertTrue(arb(result['failure_log2_upper']) < -63)
            self.assertFalse(result['whole_code_certificate'])
            for row in result['classes']:
                self.assertEqual(row['source_multisets_enumerated'], 0)
                self.assertTrue(arb(row['failure_log2_upper']) < -40)

    def test_common_pair_seventeen_group_receipt_keeps_scope_explicit(self):
        with ctx.workprec(384):
            result = replay_common_pair(17, 384)
            self.assertTrue(arb(result['failure_log2_upper']) < -67)
            self.assertEqual(result['source_row_pair_choices'], 28)
            self.assertEqual(result['constituent_pair_geometry'], 'one shared pair')
            self.assertFalse(result['full_homogeneous_group_class'])
            self.assertFalse(result['whole_code_certificate'])
        with self.assertRaisesRegex(ArithmeticError, 'did not reach 40 bits'):
            occupancy_envelope_receipt(17, '996548217/1000000000', '-127/2', 256)


if __name__ == '__main__':
    unittest.main()

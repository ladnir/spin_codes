"""Tiny exact controls; no production fractional gate is run by these tests."""
from fractions import Fraction as F
from math import comb, log2
import unittest
from unittest.mock import patch

from hd_two_sum_input_weight import exact_support_count, marked_count
from hd_two_sum_pre_toeplitz_pairs import pre_pair_support_masses
import hd_two_sum_pre_toeplitz_fractional as target


def labelled_reference(active, matching, output, group, marker=F(1)):
    total = F(0)
    for middle in range(1 << group):
        mass = F(1)
        for left, right in matching:
            a = ((active >> left) & 1) | (((active >> right) & 1) << 1)
            b = ((middle >> left) & 1) | (((middle >> right) & 1) << 1)
            mass *= pre_pair_support_masses(a, input_marker=marker)[b]
        total += (mass/255**middle.bit_count()
                  *exact_support_count(middle.bit_count(), (middle & output).bit_count(), output.bit_count(), group))
    return total


class FractionalPairTests(unittest.TestCase):
    def test_type_partition_by_weight(self):
        for group, arity in ((8, 4), (16, 8)):
            for internal in range(arity//2+1):
                rows = target.target_types(arity, internal, group)
                self.assertEqual([sum(m for t, m, *_ in rows if t.bit_count() == w) for w in range(group+1)],
                                 [comb(group, w) for w in range(group+1)])
        self.assertEqual([len(target.target_types(8, u)) for u in range(5)], [165, 756, 1260, 1000, 225])

    def test_one_fixed_A_matching_geometry(self):
        for internal in range(5):
            matching = target.matching_for_internal(8, internal)
            self.assertEqual(sorted(x for pair in matching for x in pair), list(range(16)))
            self.assertEqual(sum(x < 8 and y < 8 for x, y in matching), internal)
        self.assertFalse(target.scope_metadata()['row_set_resampled'])

    def test_positive_DP_against_labelled_middle_enumeration(self):
        matching = target.matching_for_internal(4, 1, 8)
        for output in (0, 1, 15, 51, 85, 255):
            for marker in (F(1), F(3, 2)):
                self.assertEqual(target.fixed_mask_mass(15, matching, output, group=8, input_marker=marker),
                                 labelled_reference(15, matching, output, 8, marker))

    def test_diagonal_reduces_to_old_fixed_A(self):
        for internal in range(3):
            matching = target.matching_for_internal(4, internal, 8)
            for output, *_ in target.target_types(4, internal, 8):
                for marker in (F(1), F(3, 2)):
                    self.assertEqual(target.fixed_mask_mass(15, matching, output, group=8,
                                                           input_marker=marker, mode='diagonal'),
                                     marked_count(4, (15 & output).bit_count(), output.bit_count(), marker, 8))

    def test_grouped_exact_mass_normalization(self):
        for mode in ('iid-majorant', 'exact-support'):
            rows = target.exact_type_masses(4, 1, group=8, input_marker=F(3, 2), mode=mode)
            expected = F(1)
            for left, right in target.matching_for_internal(4, 1, 8):
                active = int(left < 4)+2*int(right < 4)
                expected *= sum(pre_pair_support_masses(active, input_marker=F(3, 2), mode=mode))
            self.assertEqual(sum(m*new for _, m, new, _ in rows), expected)
            self.assertEqual(sum(m*old for _, m, _, old in rows), (1+255*F(3, 2))**4)

    def test_multiplicity_is_not_powered(self):
        new, old = target.fractional_coefficients(((1, 7, F(3), F(2)),), group=8, alpha=F(1, 5), columns=16)
        self.assertAlmostEqual(new[1], log2(7)+3.2*log2(3))
        self.assertAlmostEqual(old[1], log2(7)+3.2)
        self.assertNotAlmostEqual(new[1], 3.2*log2(21))
        self.assertEqual(new[0], -float('inf'))

    def test_alpha_one_over_columns_scalar_reproduces_raw_moment(self):
        rows = target.exact_type_masses(4, 1, group=8)
        new, old = target.fractional_coefficients(rows, group=8, alpha=F(1, 16))
        for z in (F(1, 2), F(49, 50)):
            self.assertAlmostEqual(target.scalar_summary(new, z),
                                   log2(float(sum(m*mass*z**t.bit_count() for t, m, mass, _ in rows))), places=11)
            self.assertAlmostEqual(target.scalar_summary(old, z),
                                   log2(float(sum(m*mass*z**t.bit_count() for t, m, _, mass in rows))), places=11)

    def test_scope_and_mark_guards(self):
        for values in ((16, 8, -1), (16, 8, 5), (8, 7, 0), (8, True, 0)):
            with self.assertRaises(ValueError):
                target.checked_scope(*values)
        matching = target.matching_for_internal(4, 1, 8)
        with self.assertRaises(ValueError):
            target.fixed_mask_mass(15, matching, 0, group=8, input_marker=F(1, 2))
        with self.assertRaises(ValueError):
            target.fractional_coefficients([], alpha=F(6, 5))

    def test_zero_mass_is_exactly_one(self):
        for internal in range(5):
            matching = target.matching_for_internal(8, internal)
            self.assertEqual(target.fixed_mask_mass(255, matching, 0), 1)

    def test_driver_fixed_scope_and_partial_budget(self):
        marks = [dict(active_groups=None, output_marker='1/2', template='toy', log2_activity_marker='0')]
        toy = ((0, 1, F(1), F(1)), (1, 7, F(2), F(3)))
        with patch.object(target, 'frozen_scalar_marks', return_value=marks), \
                patch.object(target, 'exact_type_masses', return_value=toy) as contraction:
            result = target.run()
        self.assertEqual(contraction.call_count, 5)
        self.assertEqual([row['internal_pairs'] for row in result['results']], list(range(5)))
        self.assertTrue(result['complete_requested_cases'])
        self.assertFalse(result['whole_code_certificate'])
        with patch.object(target, 'frozen_scalar_marks', return_value=marks), \
                patch.object(target, 'exact_type_masses', side_effect=[toy, TimeoutError('test guard')]):
            partial = target.run()
        self.assertEqual(len(partial['results']), 1)
        self.assertFalse(partial['complete_requested_cases'])
        self.assertEqual(partial['status'], 'compute-budget-exceeded')
        for budget in (True, 0, 31):
            with self.assertRaises(ValueError):
                target.run(budget)


if __name__ == '__main__':
    unittest.main()

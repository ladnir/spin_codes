"""Tiny exact checks only: no state contraction or regional provider."""
from collections import Counter
from itertools import product
import unittest

import hd_lambda_scheduled_relations as scheduled
from hd_lambda_outer import column_branches as unscheduled_branches
from hd_pair_scheduled_outer import column_permutations, permute_mask
from hd_quaternary_outer_relations import (
    category_data as raw_categories, multiply_gf4)


def direct_scheduled(values, source, schedule, column, lambda_symbol):
    """Independent full GF4 matrix evaluation with explicit P and Q wiring."""
    before, after = column_permutations(schedule, column)
    x = [0]*16
    for label, value in zip(source, values):
        x[before[label]] = value
    y = [0]*16
    for r in range(16):
        value = 0
        for i in range(16):
            coefficient = lambda_symbol ^ int(r == i) ^ ((r & i).bit_count() & 1)
            value ^= multiply_gf4(coefficient, x[i])
        y[after[r]] = value
    return sum(bool(value) << r for r, value in enumerate(y))


def input_kernel(values):
    """Enumerate all coefficient vectors vanishing on one GF4 input tuple."""
    result = 0
    for packed in range(4**len(values)):
        value = 0
        for i, x in enumerate(values):
            value ^= multiply_gf4(packed >> (2*i) & 3, x)
        if not value:
            result |= 1 << packed
    return result


class ScheduledLambdaTests(unittest.TestCase):
    def test_stable_categories_all_schedules_and_columns(self):
        for source in ((15,), (12, 3), (8, 1, 14), (15, 2, 9, 0)):
            baseline = raw_categories(len(source))
            for schedule in scheduled.SCHEDULES:
                for column in (0, 1, 7, 15):
                    data = scheduled.column_data(source, schedule, column)
                    self.assertEqual(data['original_source'], source)
                    self.assertEqual(data['original_source_mask'], sum(1 << i for i in source))
                    for key in ('kernel_rrefs', 'kernel_bitsets', 'category_input_masks',
                                'category_input_weights', 'category_counts', 'category_order'):
                        self.assertEqual(data[key], baseline[key])
                    self.assertEqual(sum(branch[2] for branch in data['branches']), 256**len(source))
                    self.assertEqual(tuple(branch[1] for branch in data['branches']), data['category_output_masks'])

    def test_a1_a2_match_exact_unscheduled_counts_after_P_and_Q(self):
        sources = ((0,), (15,), (0, 1), (12, 3), (5, 14))
        for source, schedule, column, lam in product(sources, scheduled.SCHEDULES,
                                                    (0, 1, 7, 15), scheduled.LAMBDA_VALUES):
            before, after = column_permutations(schedule, column)
            mapped = tuple(before[i] for i in source)
            expected = Counter()
            for weight, mask, count in unscheduled_branches(mapped, lambda_value=lam):
                expected[weight, permute_mask(mask, after)] += count
            self.assertEqual(scheduled.joint_counts(source, schedule, column, lambda_value=lam), expected)

    def test_GF4_exhaustive_a1_to_a3_by_individual_kernel(self):
        sources = ((13,), (9, 2), (11, 0, 7))
        for source, schedule, column, lam in product(sources, scheduled.SCHEDULES,
                                                    (2, 13), scheduled.LAMBDA_VALUES):
            data = scheduled.column_data(source, schedule, column, lambda_value=lam, q=4)
            indices = {kernel: i for i, kernel in enumerate(data['kernel_bitsets'])}
            observed_counts = [0]*len(indices)
            merged = Counter()
            for values in product(range(4), repeat=len(source)):
                index = indices[input_kernel(values)]
                output = direct_scheduled(values, source, schedule, column,
                                          2 if lam == 0xbc else 3)
                original_mask = sum(bool(value) << i for i, value in enumerate(values))
                self.assertEqual(data['category_input_masks'][index], original_mask)
                self.assertEqual(data['category_output_masks'][index], output)
                observed_counts[index] += 1
                merged[original_mask.bit_count(), output] += 1
            self.assertEqual(tuple(observed_counts), data['category_counts'])
            self.assertEqual(scheduled.joint_counts(source, schedule, column, lambda_value=lam, q=4), merged)

    def test_source_reordering_relabels_activity_not_physical_map(self):
        source, reordered = (0, 2, 11), (11, 0, 2)
        values, values_reordered = (1, 2, 3), (3, 1, 2)
        schedule, column = 'input-inverse-output-xor', 9
        for chosen, word in ((source, values), (reordered, values_reordered)):
            data = scheduled.column_data(chosen, schedule, column, q=4)
            index = data['kernel_bitsets'].index(input_kernel(word))
            self.assertEqual(data['category_output_masks'][index],
                             direct_scheduled(values, source, schedule, column, 2))
        self.assertEqual(scheduled.joint_counts(source, schedule, column),
                         scheduled.joint_counts(reordered, schedule, column))

    def test_explicit_scope_and_zero_count_retention(self):
        data = scheduled.column_data((0, 3, 8, 15), 'input-xor', 10, q=4)
        self.assertEqual(len(data['branches']), 529)
        self.assertTrue(any(count == 0 for _, _, count in data['branches']))
        self.assertEqual(sum(count for _, _, count in data['branches']), 4**4)
        self.assertEqual(data['metadata']['lambda_aes_byte'], 0xbc)
        self.assertEqual(data['metadata']['label_field_polynomial'], 0x13)
        for key in ('whole_code_certificate', 'crossed_pair_layer_included',
                    'same_relation_kernel_across_columns_required', 'rs_normalization_included'):
            self.assertFalse(data['metadata'][key])

    def test_guards(self):
        for source in ((), (2, 2), (0, 1, 2, 3, 4), (16,)):
            with self.assertRaises(ValueError):
                scheduled.column_masks(source, 'static', 0)
        for lam in (2, 3, 1):
            with self.assertRaises(ValueError):
                scheduled.column_masks((0,), 'static', 0, lambda_value=lam)
        with self.assertRaises(ValueError):
            scheduled.column_branches((0,), 'static', 0, q=8)
        with self.assertRaises(ValueError):
            scheduled.column_masks((0,), 'unprescribed', 0)
        with self.assertRaises(ValueError):
            scheduled.column_masks((0,), 'static', 16)
        with self.assertRaises(ValueError):
            scheduled.category_data(0)


if __name__ == '__main__':
    unittest.main()

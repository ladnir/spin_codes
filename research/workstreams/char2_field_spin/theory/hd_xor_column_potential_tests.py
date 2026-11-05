"""Independent coverage, chronology, overlap, and directed-replay checks."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb, log2
import unittest

import numpy as np
from flint import arb, ctx

from hd_xor_column_probe import column_support_counts
from hd_xor_column_potential import (FULL, PAIRS, SAVED_WITNESSES,
    best_support_marks, case_plans, category_costs, chronological_words,
    column_matrices, group_log_moment, mask_code, output_support_law,
    overlap_counts, receipt, replay, arb_chronological_words,
    factored_column_vectors, factored_word_vector, multiset_plans,
    multiset_column_matrices, multiset_support_marks, multiset_replay,
    multiset_code_batches, streaming_multiset_screen)


class ColumnPotentialTests(unittest.TestCase):
    def test_exact_support_law(self):
        for source in (0, *(1 << i for i in range(8)), *PAIRS):
            counts = column_support_counts(256, 8, source)
            expected = tuple((target, Fraction(count, 255**source.bit_count()))
                             for target, count in enumerate(counts) if count)
            self.assertEqual(tuple(sorted(output_support_law(source))), expected)

    def test_all_784_pair_assignments_and_support_subsets_are_covered(self):
        plans = case_plans(2)
        laws = {category: {law for first, second, law in rows} for category, rows in plans.items()}
        self.assertEqual(len(plans[0, 2]), 28**2)
        checked = 0
        for first_pair in PAIRS:
            first_bits = [1 << i for i in range(8) if first_pair >> i & 1]
            for second_pair in PAIRS:
                second_bits = [1 << i for i in range(8) if second_pair >> i & 1]
                for first, second in product((0, *first_bits, first_pair), (0, *second_bits, second_pair)):
                    category = (int(first.bit_count() == 1)+int(second.bit_count() == 1),
                                int(first.bit_count() == 2)+int(second.bit_count() == 2))
                    expected = {}
                    left = column_support_counts(256, 8, first)
                    right = column_support_counts(256, 8, second)
                    for a, count_a in enumerate(left):
                        if not count_a:
                            continue
                        for b, count_b in enumerate(right):
                            if count_b:
                                code = mask_code(a)+mask_code(b)
                                expected[code] = expected.get(code, Fraction(0))+Fraction(
                                    count_a*count_b, 255**(first.bit_count()+second.bit_count()))
                    self.assertIn(tuple(sorted(expected.items())), laws[category])
                    checked += 1
        self.assertEqual(checked, 28**2*16)

    def test_masks_add_to_actual_occupancies(self):
        for first, second in product((0, 1, 128, 3, 129, 255), repeat=2):
            code = mask_code(first)+mask_code(second)
            for index in range(8):
                self.assertEqual(code//3**index % 3,
                                 ((first >> index) & 1)+((second >> index) & 1))

    def test_words_keep_noncommuting_order(self):
        regions = np.array([[[1, 0], [0, .125]], [[0, .5], [.0625, .25]],
                            [[.125, .375], [.03125, .5]]])
        with np.errstate(divide='ignore'):
            words = chronological_words(np.log2(regions))
        for counts in ((1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1),
                       (2, 1, 0, 2, 0, 1, 2, 1)):
            expected = np.eye(2)
            for count in counts:
                expected = expected @ regions[count]
            code = sum(count*3**index for index, count in enumerate(counts))
            np.testing.assert_allclose(np.exp2(words[code]), expected, rtol=1e-14)
        self.assertFalse(np.allclose(np.exp2(words[1]), np.exp2(words[3**7])))

    def test_overlap_counts_against_gf4_codewords(self):
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
        words = [tuple(a ^ multiply(b, x) for x in range(3))
                 for a, b in product(range(4), repeat=2) if a or b]
        expected = Counter()
        for left, right in product(words, repeat=2):
            single = sum(bool(a) != bool(b) for a, b in zip(left, right))
            double = sum(bool(a) and bool(b) for a, b in zip(left, right))
            expected[single, double] += 1
        self.assertEqual(dict(overlap_counts(4, 3, 2)), dict(expected))
        self.assertEqual(sum(count for degree, count in overlap_counts()), (256**16-1)**2)
        self.assertAlmostEqual(group_log_moment(0., 0.), 2*log2(256**16-1), places=12)

    def test_support_marks_satisfy_every_category(self):
        costs, _ = category_costs(column_matrices(.99875, 2), -62.375)
        first, second = best_support_marks(costs, 2)
        for (singles, doubles), cost in costs.items():
            self.assertGreaterEqual(singles*first+doubles*second+1e-12, cost)

    def test_mixture_occurs_before_potential(self):
        # A case's two branches may have different maximizing states.
        potential = -2.
        branches = np.array([[[.5, 0], [0, .125]], [[.125, 0], [0, .5]]])
        average = branches.mean(axis=0)
        with np.errstate(divide='ignore'):
            costs, _ = category_costs({(0, 1): np.log2(average[None])}, potential)
        self.assertAlmostEqual(2**costs[0, 1], .3125)
        self.assertLess(2**costs[0, 1], .5)

    def test_saved_receipts(self):
        with ctx.workprec(384):
            for groups, target, cases in ((1, -140, 37), (2, -108, 1109)):
                result = replay(groups, 384)
                self.assertTrue(arb(result['failure_log2_upper']) < target)
                self.assertEqual(result['checked_support_cases'], cases)
                self.assertEqual(result['source_row_pair_choices'], 28**groups)
                self.assertEqual(result['projective_normalizations'], 1)
                self.assertFalse(result['whole_code_certificate'])

    def test_false_marks_are_rejected(self):
        z, potential, first, second = SAVED_WITNESSES[2]
        with self.assertRaisesRegex(ArithmeticError, 'potential failed'):
            receipt(2, z, potential, '-100', '-100')

    def test_group_multisets_preserve_all_ordered_two_group_laws(self):
        reduced = multiset_plans(2)
        for (single, double), plans in case_plans(2).items():
            expected = {law for first, second, law in plans}
            actual = set()
            for codes in reduced[single, double][0]:
                law = {}
                for branch, code in enumerate(codes):
                    weight = Fraction(254**(double-branch.bit_count()), 255**double)
                    law[int(code)] = law.get(int(code), Fraction(0))+weight
                actual.add(tuple(sorted(law.items())))
            self.assertEqual(actual, expected)

    def test_multiset_category_counts_and_repeated_edge_probabilities(self):
        for groups in (3, 4):
            plans = multiset_plans(groups)
            self.assertEqual(sum(len(codes) for codes, weights in plans.values()), comb(36+groups, groups))
            for (single, double), (codes, weights) in plans.items():
                self.assertEqual(len(codes), comb(7+single, single)*comb(27+double, double))
            # The first edge multiset repeats the same pair in every group.
            # Distinct coin branches landing on one code must add, not collapse.
            codes = plans[0, groups][0][0]
            law = {}
            for branch, code in enumerate(codes):
                weight = Fraction(254**(groups-branch.bit_count()), 255**groups)
                law[int(code)] = law.get(int(code), Fraction(0))+weight
            full, pair = mask_code(FULL, groups+1), mask_code(PAIRS[0], groups+1)
            expected = {(groups-cancelled)*full+cancelled*pair:
                        Fraction(comb(groups, cancelled)*254**(groups-cancelled), 255**groups)
                        for cancelled in range(groups+1)}
            self.assertEqual(law, expected)
            self.assertEqual(sum(law.values()), 1)

    def test_factored_vectors_keep_noncommuting_eight_region_order(self):
        regions = [(arb(1), arb(0), arb(0), arb(1)/8),
                   (arb(0), arb(1)/2, arb(1)/16, arb(1)/4),
                   (arb(1)/8, arb(3)/8, arb(1)/32, arb(1)/2)]
        potential = arb(1)/4
        prefixes, suffixes = factored_column_vectors(regions, potential)
        words = arb_chronological_words(regions, 8)
        self.assertEqual(len(prefixes), 3**4)
        for code, matrix in enumerate(words):
            actual = factored_word_vector(code, prefixes, suffixes)
            self.assertEqual(actual, (matrix[0]+matrix[1]*potential, matrix[2]+matrix[3]*potential))
        self.assertNotEqual(factored_word_vector(1, prefixes, suffixes),
                            factored_word_vector(3**7, prefixes, suffixes))

    def test_multiset_screen_matches_ordered_screen(self):
        original, _ = category_costs(column_matrices(.99875, 2), -62.375)
        reduced, _ = category_costs(multiset_column_matrices(.99875, 2), -62.375)
        for category in original:
            self.assertAlmostEqual(original[category], reduced[category], places=12)
        first, second = multiset_support_marks(reduced)
        for (single, double), cost in reduced.items():
            self.assertGreaterEqual(single*first+double*second+1e-12, cost)

    def test_streaming_batches_match_exact_multiset_codes(self):
        for groups in (2, 3):
            for (single, double), (codes, weights) in multiset_plans(groups).items():
                streamed = np.concatenate(list(multiset_code_batches(groups, single, double, batch_size=13)))
                np.testing.assert_array_equal(streamed, codes)

    def test_streaming_screen_matches_materialized_three_group_screen(self):
        z, potential = .99859, -62.40625
        costs, _ = category_costs(multiset_column_matrices(z, 3), potential)
        result = streaming_multiset_screen(3, z, potential)
        for category, cost in costs.items():
            self.assertAlmostEqual(result['category_log2_costs'][str(category)], cost, places=11)
        self.assertEqual(result['checked_multiset_cases'], comb(39, 3))
        self.assertFalse(result['directed'])

    def test_saved_three_and_four_group_receipts(self):
        with ctx.workprec(384):
            for groups, target in ((3, -76), (4, -44)):
                result = multiset_replay(groups, 384)
                self.assertTrue(arb(result['failure_log2_upper']) < target)
                self.assertEqual(result['checked_support_cases'], comb(36+groups, groups))
                self.assertEqual(result['source_row_pair_choices'], 28**groups)
                self.assertEqual(result['cached_prefix_matrices'], (groups+1)**4)
                self.assertEqual(result['cached_suffix_vectors'], (groups+1)**4)
                self.assertEqual(result['projective_normalizations'], 1)
                self.assertTrue(result['group_permutation_reduction'])
                self.assertFalse(result['region_permutation_reduction'])
                self.assertFalse(result['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()

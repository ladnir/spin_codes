import unittest
from itertools import combinations

import hd_bit_syndrome_column_permutations as p
from hd_bit_syndrome_outer import affine_rank


class ColumnPermutationTests(unittest.TestCase):
    def test_field_inverse_and_permutations(self):
        for a in range(1, 16):
            self.assertEqual(p.mul(a, p.INVERSE[a]), 1)
        for j, permutation in enumerate(p.PERMUTATIONS):
            self.assertEqual(sorted(permutation), list(range(16)))
            for i, r in enumerate(permutation):
                self.assertEqual(p.INVERSE[r] ^ j, i)

    def test_all_four_set_plane_conditions_match_binary_rank(self):
        for choice in combinations(range(16), 4):
            mask = sum(1 << i for i in choice)
            exact = tuple(j for j in range(16) if affine_rank(p.map_mask(mask, j)) == 3)
            self.assertEqual(p.plane_columns(mask), exact)

    def test_count_scope_and_RS_implication(self):
        result = p.exact_summary()
        self.assertEqual(result['fixed_source_sets_checked'], 1820)
        self.assertEqual(result['quiet_column_histogram'], {0: 120, 1: 1200, 2: 480, 4: 20})
        self.assertEqual(result['maximum_quiet_columns'], 4)
        self.assertEqual(result['minimum_nonquiet_columns_for_four_active_RS_rows'], 5)
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['old_receipts_transfer'])

    def test_guards(self):
        for source in (0, 7, -1, 65536):
            with self.assertRaises(ValueError):
                p.plane_columns(source)
        with self.assertRaises(ValueError):
            p.map_mask(15, 16)

    def test_singleton_coloop_activity_scope(self):
        # A singleton is always a coloop, hence no nonzero quiet column.
        result = p.quiet_activity_incidence((1,), max_seconds=5)
        row = result['results'][0]
        self.assertEqual(row['minimum_activity_histogram'], {0: 16})
        self.assertEqual(row['universal_nonquiet_column_floor'], 9)
        self.assertFalse(result['sufficient_for_all_quiet_RS_word'])

    def test_complete_five_to_seven_incidence_regression(self):
        result = p.quiet_activity_incidence((5, 6, 7), max_seconds=10)
        self.assertTrue(result['all_requested_arities_complete'])
        expected = (
            (4368, {3: 960, 4: 2208, 5: 960, 7: 240}, 4368, 2),
            (8008, {9: 1920, 10: 4560, 11: 1200, 12: 328}, 0, 0),
            (11440, {12: 720, 13: 3520, 14: 6960, 15: 240}, 0, 0),
        )
        for row, (count, histogram, excluded, floor) in zip(result['results'], expected):
            self.assertEqual(row['labelled_sets_checked'], count)
            self.assertEqual(row['minimum_activity_histogram'], histogram)
            self.assertEqual(row['excluded_all_quiet_sets'], excluded)
            self.assertEqual(row['universal_nonquiet_column_floor'], floor)
            self.assertEqual(row['not_excluded_sets'], count-excluded)


if __name__ == '__main__':
    unittest.main()

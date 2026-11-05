import unittest

from hd_bit_syndrome_column_dimensions import column_ranks, dimension_summary


class ColumnDimensionTests(unittest.TestCase):
    def test_empty_full_and_four_row_rank_ranges(self):
        self.assertEqual(column_ranks(0), (0,)*16)
        self.assertEqual(column_ranks(65535), (5,)*16)
        ranks = column_ranks(195)
        self.assertEqual(ranks.count(3), 4)
        self.assertTrue(all(r in (3,4) for r in ranks))

    def test_exact_four_row_histogram_and_scope(self):
        result = dimension_summary((4,))
        row = result['records'][0]
        self.assertTrue(result['all_requested_arities_complete'])
        self.assertEqual(row['quiet_dimension_sum_histogram'], {0:120,1:1200,2:480,4:20})
        self.assertEqual(row['labelled_sets_checked'], 1820)
        self.assertFalse(result['global_setup_probability_computed'])
        self.assertFalse(result['whole_code_certificate'])

    def test_nine_row_rank_is_always_five_in_selected_columns(self):
        for mask in (511, 1022, 65535 ^ 127, 65535 ^ 5461):
            self.assertEqual(mask.bit_count(), 9)
            self.assertEqual(column_ranks(mask), (5,)*16)

    def test_guards(self):
        for arities in ((), (0,), (4,4), (17,)):
            with self.assertRaises(ValueError):
                dimension_summary(arities)
        with self.assertRaises(ValueError):
            dimension_summary((4,), max_seconds=11)


if __name__ == '__main__':
    unittest.main()

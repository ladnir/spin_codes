"""Tiny exact scheduled-relation tests; no state/proof contraction."""
from collections import Counter
from itertools import product
import unittest

import numpy as np

import hd_pair_relation_data as data
from hd_binary_outer_relations import joint_counts
from hd_pair_scheduled_outer import column_branches


def brute(forms, arity, q):
    result = Counter()
    for values in product(range(q), repeat=arity):
        mask = 0
        for row, form in enumerate(forms):
            value = 0
            for index, x in enumerate(values):
                if form >> index & 1:
                    value ^= x
            mask |= bool(value) << row
        result[sum(x != 0 for x in values), mask] += 1
    return result


class RelationDataTests(unittest.TestCase):
    def test_counts_and_original_input_masks(self):
        for arity, expected in ((1, 2), (2, 5), (3, 16), (4, 67)):
            for q in (2, 4, 256):
                categories = data.category_data(arity, q)
                self.assertEqual(len(categories), expected)
                self.assertEqual(sum(row[3] for row in categories), q**arity)
                for kernel, mask, weight, _ in categories:
                    self.assertEqual(mask, sum((not (kernel >> (1 << i) & 1)) << i for i in range(arity)))
                    self.assertEqual(weight, mask.bit_count())
        self.assertTrue(any(row[3] == 0 for row in data.category_data(4, 4)))

    def test_exact_gf4_selected_a3_a4_all_schedules(self):
        schedules = ('static', 'input-xor', 'input-inverse', 'output-xor', 'input-inverse-output-xor')
        for source, schedule, column in product(((0, 1, 7), (0, 2, 8, 15)), schedules, (1, 9)):
            forms = data.row_forms(source, schedule, column)
            masks = data.column_masks(source, schedule, column)
            observed = Counter()
            for (_, _, weight, count), mask in zip(data.category_data(len(source), 4), masks):
                if count:
                    observed[weight, int(mask)] += count
            self.assertEqual(observed, brute(forms, len(source), 4))
            self.assertEqual(observed, joint_counts(forms, len(source), q=4))

    def test_a1_a2_agree_frozen_exact_columns(self):
        for source, schedule, column in product(((0,), (15,), (0, 1), (3, 12)),
                ('static', 'input-xor', 'input-inverse-output-xor'), range(16)):
            observed = Counter()
            for (_, _, weight, count), mask in zip(data.category_data(len(source)),
                    data.column_masks(source, schedule, column)):
                observed[weight, int(mask)] += count
            expected = Counter()
            for weight, mask, count in column_branches(source, schedule, column):
                expected[weight, mask] += count
            self.assertEqual(observed, expected)

    def test_chunk_pair_layout_and_invariant_categories(self):
        sources = ((0, 1, 2), (4, 8, 15))
        prepared = data.branch_data(3, 'input-inverse-output-xor', sources=sources)
        self.assertEqual(prepared['first'].shape, (8, 512))
        self.assertEqual(prepared['counts'].shape, (2, 256))
        self.assertEqual(prepared['sources'], sources)
        cats = data.category_data(3)
        for source_index, source in enumerate(sources):
            for block in (0, 3, 7):
                left = data.column_masks(source, 'input-inverse-output-xor', 2*block)
                right = data.column_masks(source, 'input-inverse-output-xor', 2*block+1)
                for i, j in product(range(16), repeat=2):
                    index = 16*i+j
                    self.assertEqual(prepared['first'][block, source_index*256+index], left[i])
                    self.assertEqual(prepared['second'][block, source_index*256+index], right[j])
                    self.assertEqual(prepared['counts'][source_index, index], cats[i][3]*cats[j][3])
                    self.assertEqual(prepared['input_weights'][source_index, index], cats[i][2]+cats[j][2])
                    self.assertEqual(prepared['input_union_masks'][source_index, index], cats[i][1] | cats[j][1])

    def test_a4_counts_do_not_wrap_uint64_and_keep_zero_toy_categories(self):
        prepared = data.branch_data(4, sources=((0, 1, 2, 3),))
        self.assertEqual(prepared['branch_count'], 4489)
        self.assertEqual(sum(prepared['counts'][0]), 256**8)
        self.assertEqual(prepared['counts'].dtype, np.dtype(object))
        toy = data.branch_data(4, sources=((0, 1, 2, 3),), q=4)
        self.assertEqual(toy['branch_count'], 4489)
        self.assertEqual(sum(toy['counts'][0]), 4**8)
        self.assertTrue(np.any(np.isneginf(toy['log_counts'])))

    def test_bounded_iterator_and_guards(self):
        chunks = list(data.iter_branch_data(1, batch_size=5))
        self.assertEqual(tuple(source for chunk in chunks for source in chunk['sources']), data.source_sets(1))
        iterator = data.iter_branch_data(3, batch_size=2)
        self.assertEqual(next(iterator)['sources'], data.source_sets(3)[:2])
        self.assertEqual(next(iterator)['sources'], data.source_sets(3)[2:4])
        for arity in (3, 4):
            with self.assertRaises(ValueError):
                data.branch_data(arity)
        with self.assertRaises(ValueError):
            next(data.iter_branch_data(4, batch_size=2))
        with self.assertRaises(ValueError):
            data.branch_data(3, sources=((0, 1, 2), (0, 1, 2)))
        self.assertFalse(data.metadata(3, 'input-xor')['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()

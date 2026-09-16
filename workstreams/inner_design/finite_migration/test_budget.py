import copy
from fractions import Fraction as F
import unittest

import budget_dense as budget
import verify_budget


def example():
    root = budget.geometry.box(64, 512, F(0), F(1), {})
    left, right = budget.geometry.children(root, 'v')
    left['power'], right['power'] = -44, -45
    leaves, splits = {'0': left, '1': right}, {'': 'v'}
    return dict(status=budget.STATUS, minimum=64, instance={'outer_rows': 512},
                leaves=leaves, splits=splits,
                coverage=budget.geometry.check_partition(leaves, splits, 64, 512),
                upper=budget.model.base.encode(F(2)**-44 + F(2)**-45),
                budget_bits=42, budget_met=True)


class BudgetTests(unittest.TestCase):
    def test_valid_union_without_eighty_bits_per_leaf(self):
        record = example()
        self.assertEqual(budget.checked_union(record, 64, 512), 3 * F(2)**-45)

    def test_rejects_tampered_union_or_missing_coverage(self):
        for change in ('sum', 'coverage', 'status', 'power', 'target'):
            record = copy.deepcopy(example())
            if change == 'sum':
                record['upper'] = budget.model.base.encode(F(2)**-80)
            elif change == 'coverage':
                del record['leaves']['1']
            elif change == 'status':
                record['status'] = 'point-check'
            elif change == 'power':
                record['leaves']['0']['power'] = -44.0
            else:
                record['budget_bits'] = 80
            with self.assertRaises(Exception, msg=change):
                budget.checked_union(record, 64, 512)

    def test_full_union_independently_enforces_forty_bits(self):
        self.assertEqual(verify_budget.full_union([F(2)**-43] * 3), 3 * F(2)**-43)
        with self.assertRaises(AssertionError):
            verify_budget.full_union([F(2)**-41, F(2)**-42, F(2)**-42])
        with self.assertRaises(AssertionError):
            verify_budget.full_union([F(2)**-42, F(2)**-42, F(0)])


if __name__ == '__main__':
    unittest.main()

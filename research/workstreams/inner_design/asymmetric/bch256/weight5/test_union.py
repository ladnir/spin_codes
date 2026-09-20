"""Reject gaps, altered unions, and incomplete dense certificates."""
from copy import deepcopy
import unittest
import verify_full as verify
model = verify.model


class UnionTests(unittest.TestCase):
    def test_sparse_gap_rejected(self):
        bank = model.base.read(verify.candidate.HERE/'RANGE_M20.json')
        verify.sparse_union(bank,511)
        changed = deepcopy(bank)
        for job in changed['witnesses']:
            for row in job['rows']:
                if row['occupation'] == 2:
                    row['power'] = 0
        with self.assertRaises(AssertionError):
            verify.sparse_union(changed,511)

    def test_union_tamper_rejected(self):
        bank = model.base.read(verify.candidate.HERE/'RANGE_M20.json')
        bank['upper'] = model.base.encode(0)
        with self.assertRaises(AssertionError):
            verify.sparse_union(bank,511)

    def test_dense_gap_rejected(self):
        cover = model.base.read(verify.candidate.HERE/'COVER_M20.json')
        verify.dense_union(cover,512,8192)
        del cover['leaves'][next(iter(cover['leaves']))]
        with self.assertRaises((AssertionError,ValueError)):
            verify.dense_union(cover,512,8192)

    def test_weak_leaf_rejected(self):
        cover = model.base.read(verify.candidate.HERE/'COVER_M20.json')
        next(iter(cover['leaves'].values()))['power'] = -39
        with self.assertRaises(AssertionError):
            verify.dense_union(cover,512,8192)


if __name__ == '__main__':
    unittest.main()

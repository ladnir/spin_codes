import itertools
import unittest

import numpy as np

import entropy_type_split_v1 as split
import refine_bch_dense_v1 as dense


class TypeSplitTests(unittest.TestCase):
    def test_recursive_cover_preserves_every_integer_type(self):
        for total in (3, 8):
            pending = [(np.zeros(4, dtype=int), np.full(4, total, dtype=int))]
            points = []
            while pending:
                lo, hi = pending.pop()
                corners = dense.transfer.typed.vertices(lo, hi, total)
                children = split.split_box(lo, hi, corners)
                if not children:
                    points.append(tuple(corners[0]))
                    continue
                nonempty = [(a, b) for a, b in children if dense.lattice_count(a.tolist(), b.tolist(), total)]
                self.assertEqual(sum(dense.lattice_count(a.tolist(), b.tolist(), total) for a, b in nonempty),
                                 dense.lattice_count(lo.tolist(), hi.tolist(), total))
                pending.extend(nonempty)
            expected = {x for x in itertools.product(range(total+1), repeat=4) if sum(x) == total}
            self.assertEqual(set(points), expected)
            self.assertEqual(len(points), len(expected))

    def test_zero_face_is_isolated(self):
        lo = np.array([9000, 0, 0]); hi = np.array([10000, 512, 512])
        corners = dense.transfer.typed.vertices(lo, hi, 10000)
        children = split.split_box(lo, hi, corners)
        changed = np.flatnonzero(children[0][1] < corners.max(axis=0))
        self.assertEqual(len(changed), 1)
        coordinate = changed[0]
        self.assertIn(coordinate, (1, 2))
        self.assertEqual(children[0][1][coordinate], 0)
        self.assertEqual(children[1][0][coordinate], 1)


if __name__ == '__main__':
    unittest.main()

import unittest
from hyperplane_counts import bound
from test_constraint_counts import dimensions
from dual_moments import exact_tuple_cdfs
from joint_support import span


class HyperplaneCountTests(unittest.TestCase):
    def test_exhaustive_codes_and_inflated_cdfs(self):
        for basis, n in (([3, 5], 4), ([15, 51, 85], 7), ([0x97, 0x4b, 0x2d, 0x1e], 8)):
            words = span(basis)
            exact = exact_tuple_cdfs(words, n, 4)
            dims = dimensions(words, n)
            for multiplier in (1, 3):
                caps = [[multiplier*x for x in row] for row in exact]
                for h in range(2, 5):
                    for u in range(n+1):
                        upper, _ = bound(caps, dims, 4, h, u)
                        self.assertLessEqual(exact[h-1][u], upper)
                        self.assertLessEqual(upper, caps[h-1][u])


if __name__ == '__main__':
    unittest.main()

"""Ensure deferred tables preserve fixed-composition bounds."""
import unittest

import numpy as np

import occupation_composition_lazy_v1 as lazy
import occupation_refresh_v1 as original


class LazyComposition(unittest.TestCase):
    def test_bound_matches_eager_tables(self):
        args = ({2: 2, 4: 1}, 4, [[2], [4]], [.5, .8], 12)
        eager, deferred = original.CompositionBoxes(*args), lazy.CompositionBoxes(*args)
        epoch = original.Epochs(4, 2, {2: 2, 4: 1}, [1, 0, 2, 0, 1])
        regions = original.region_logs(epoch.at(.3, 12), 4, 16, 12)
        for q in (2, 6, 12):
            for a in range(q+1):
                lo, hi = [a, 0], [q, q-a]
                expected = eager.bound(regions, lo, hi, q, 16, 6, .3)
                actual = deferred.bound(regions, lo, hi, q, 16, 6, .3)
                self.assertEqual(expected, actual)
        for probability, table in zip(args[3], deferred.binomials):
            for count in (0, 6, 12):
                np.testing.assert_array_equal(table[count], original.boxes.binomial_logs(count, probability))


if __name__ == '__main__':
    unittest.main()

from collections import Counter, defaultdict
from math import comb
import unittest

from hd_xor_group_geometry import (matching_joint_counts, matching_max_counts, matching_max_cumulative_counts,
                                    subset_full_pair_range)


def full_pairs(mask, size, parity):
    return sum(bool(mask & (1 << i)) and bool(mask & (1 << ((i+1) % size)))
               for i in range(parity, size, 2))


class GroupGeometryTests(unittest.TestCase):
    def test_every_small_indexed_cycle_support(self):
        for size in (4, 6, 8, 10, 12):
            exact = defaultdict(Counter)
            joint = defaultdict(Counter)
            for mask in range(1 << size):
                exact[mask.bit_count()][max(full_pairs(mask, size, p) for p in (0, 1))] += 1
                joint[mask.bit_count()][tuple(full_pairs(mask, size, p) for p in (0, 1))] += 1
            for weight in range(size+1):
                self.assertEqual(matching_joint_counts(size, weight), joint[weight])
                counts = matching_max_counts(size, weight)
                self.assertEqual(counts, exact[weight])
                self.assertEqual(sum(counts.values()), comb(size, weight))
                cumulative = matching_max_cumulative_counts(size, weight)
                self.assertEqual(cumulative, {limit: sum(count for f, count in counts.items() if f <= limit)
                                             for limit in range(max(0, weight-size//2), weight//2+1)})

    def test_production_seventeen_partition_and_no_adjacent_supports(self):
        counts = matching_max_counts(8192, 17)
        self.assertEqual(set(counts), set(range(9)))
        self.assertEqual(sum(counts.values()), comb(8192, 17))
        numerator = 8192*comb(8192-17, 17)
        self.assertEqual(numerator % (8192-17), 0)
        self.assertEqual(counts[0], numerator//(8192-17))
        self.assertGreater(counts[0]*100, 96*comb(8192, 17))
        joint = matching_joint_counts(8192, 17)
        self.assertEqual(len(joint), 81)
        self.assertEqual(joint[8, 8], 8192)
        self.assertEqual(sum(joint.values()), comb(8192, 17))
        self.assertEqual(counts, {full: sum(count for (a, b), count in joint.items()
                                          if max(a, b) == full) for full in counts})

    def test_complement_and_cache_return_values_are_independent(self):
        sparse = matching_max_counts(8192, 17)
        self.assertEqual(matching_max_counts(8192, 8192-17),
                         {4096-17+f: count for f, count in sparse.items()})
        sparse[0] = -1
        self.assertGreater(matching_max_counts(8192, 17)[0], 0)
        self.assertEqual(matching_max_counts(8192, 0), {0: 1})
        self.assertEqual(matching_max_counts(8192, 8192), {4096: 1})
        sparse_joint = matching_joint_counts(8192, 17)
        self.assertEqual(matching_joint_counts(8192, 8192-17),
                         {(4096-17+a, 4096-17+b): count for (a, b), count in sparse_joint.items()})
        sparse_joint[0, 0] = -1
        self.assertGreater(matching_joint_counts(8192, 17)[0, 0], 0)

    def test_subset_range_is_exact_on_small_matching_supports(self):
        size = 6
        for source in range(1 << size):
            weight = source.bit_count()
            for parity in (0, 1):
                full = full_pairs(source, size, parity)
                observed = defaultdict(set)
                subset = source
                while True:
                    observed[subset.bit_count()].add(full_pairs(subset, size, parity))
                    if subset == 0:
                        break
                    subset = (subset-1) & source
                for m in range(weight+1):
                    low, high = subset_full_pair_range(weight, full, m)
                    self.assertEqual(observed[m], set(range(low, high+1)))
                self.assertEqual(subset_full_pair_range(weight, full, weight), (full, full))

    def test_invalid_geometry_and_explicit_work_guard(self):
        for size, weight in ((True, 0), (2, 1), (7, 2), (8, True), (8, -1), (8, 9), (8, 1.5)):
            with self.assertRaises(ValueError):
                matching_max_counts(size, weight)
        with self.assertRaisesRegex(ValueError, 'work budget'):
            matching_max_counts(8192, 33)
        with self.assertRaises(ValueError):
            matching_max_counts(8, 2, max_reduced_weight=True)
        self.assertEqual(sum(matching_max_counts(12, 6, max_reduced_weight=6).values()), comb(12, 6))
        for args in ((3, 2, 1), (3, 1, 4), (3, True, 1), (-1, 0, 0)):
            with self.assertRaises(ValueError):
                subset_full_pair_range(*args)


if __name__ == '__main__':
    unittest.main()

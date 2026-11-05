"""Frozen partition control; numerical primitives have separate exact tests."""
import unittest

from hd_xor_dense_holder_classes import frozen_classes, class_moments
from hd_xor_dense_holder_arity import arity_classes


class HolderArityTests(unittest.TestCase):
    def test_frozen_twelve_intersections_partition219(self):
        controls, _ = frozen_classes()
        total, seen = 0, []
        for _, _, _, masks, _ in controls:
            groups = arity_classes(masks)
            total += len(groups)
            for arity, members in groups:
                self.assertTrue(all(mask.bit_count() == arity for mask in members))
                seen.extend(members)
            self.assertEqual(sorted(m for _, group in groups for m in group), sorted(masks))
        self.assertEqual(total, 12)
        self.assertEqual(sorted(seen), [mask for mask in range(256) if mask.bit_count() >= 3])

    def test_partition_preserves_exact_moment_sum(self):
        masks = (3, 7, 10, 15)
        whole, denominator = class_moments(masks, 4, 4, 3, 2)
        parts = [class_moments(group, 4, 4, 3, 2) for _, group in arity_classes(masks)]
        self.assertTrue(all(d == denominator for _, d in parts))
        self.assertEqual(tuple(sum(part[i] for part, _ in parts) for i in range(16)), whole)


if __name__ == '__main__':
    unittest.main()

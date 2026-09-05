#!/usr/bin/env python3

import unittest

from streaming_ec_pair_collision_audit import collision_profile


class StreamingEcPairCollisionAuditTest(unittest.TestCase):
    def test_small_region_forces_repeated_collision_keys(self) -> None:
        profile = collision_profile(
            left_degree=10,
            right_degree=5,
            region_size=2,
            seed=1,
        )
        self.assertGreaterEqual(profile["maximum_cancellable_shared_edges"], 2)
        self.assertGreater(profile["repeated_projective_keys"], 0)

    def test_large_region_has_no_repeat_for_reference_seed(self) -> None:
        profile = collision_profile(
            left_degree=10,
            right_degree=5,
            region_size=1_000_003,
            seed=1,
        )
        self.assertEqual(profile["maximum_cancellable_shared_edges"], 1)
        self.assertEqual(profile["repeated_projective_keys"], 0)


if __name__ == "__main__":
    unittest.main()

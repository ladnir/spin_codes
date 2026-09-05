import unittest

from ec_placement_cap_experiment import (
    exhaustive_placement_cap,
    regional_pattern_distribution,
    trace_counts_by_equations,
)


class PlacementCapExperimentTest(unittest.TestCase):
    def test_regional_distribution_is_normalized_and_symmetric(self):
        law = regional_pattern_distribution(
            region_length=3, group_size=2, message_weight=2
        )
        self.assertAlmostEqual(sum(law.values()), 1.0)
        self.assertEqual(len({law[pattern] for pattern in law if sum(pattern) == 1}), 1)
        self.assertEqual(len({law[pattern] for pattern in law if sum(pattern) == 2}), 1)

    def test_zero_occupancy_has_one_zero_trace(self):
        counts = trace_counts_by_equations(
            occupancy=(0, 0, 0), memory=2, cutoff=0
        )
        self.assertEqual(counts[0], 1)
        self.assertEqual(sum(counts), 1)

    def test_conditional_cap_never_worsens_bound(self):
        result = exhaustive_placement_cap(
            prime=7,
            region_count=2,
            region_length=3,
            group_size=1,
            memory=2,
            message_weight=1,
            cutoff=2,
        )
        self.assertLessEqual(result.capped_mean, result.uncapped_mean)
        self.assertGreaterEqual(result.improvement_bits, 0.0)


if __name__ == "__main__":
    unittest.main()

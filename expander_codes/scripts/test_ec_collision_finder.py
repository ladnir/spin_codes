import itertools
import unittest

from ec_collision_finder import (
    effective_support,
    simulate_wrapping_candidates,
    weight_one_probability_five_draws,
    weight_three_probability_five_draws,
)


def naive_scatter(support, taps, memory):
    values = [0] * len(taps)
    for coordinate in support:
        values[coordinate] ^= 1
    for source, tap_mask in enumerate(taps):
        if not values[source]:
            continue
        for offset in range(1, memory + 1):
            if tap_mask & (1 << (offset - 1)):
                target = source + offset
                if target < len(values):
                    values[target] ^= 1
    return sum(values)


class ECCollisionFinderTests(unittest.TestCase):
    def test_effective_support_cancels_pairs(self):
        self.assertEqual(effective_support([7, 2, 7, 4, 3]), (2, 3, 4))
        self.assertEqual(effective_support([5, 5, 5, 5, 9]), (9,))

    def test_five_draw_probabilities_match_exhaustive_count(self):
        q = 4
        counts = {1: 0, 3: 0, 5: 0}
        for draws in itertools.product(range(q), repeat=5):
            counts[len(effective_support(draws))] += 1
        total = q**5
        self.assertAlmostEqual(counts[3] / total, weight_three_probability_five_draws(q))
        self.assertAlmostEqual(counts[1] / total, weight_one_probability_five_draws(q))

    def test_fixed_oldest_tap_zero_boundary(self):
        taps = [2] * 5  # memory two, only offset two is active
        self.assertEqual(
            simulate_wrapping_candidates([(0,)], taps, memory=2),
            [3],
        )

    def test_fixed_oldest_tap_stops_at_block_boundary(self):
        taps = [2] * 5  # no cyclic update from positions three or four
        self.assertEqual(
            simulate_wrapping_candidates([(3,)], taps, memory=2),
            [1],
        )

    def test_bit_sliced_simulation_matches_naive(self):
        taps = [0b101, 0b110, 0b100, 0b111, 0b101, 0b100, 0b110]
        supports = [(0, 4), (1, 3, 6), (2,)]
        expected = [naive_scatter(support, taps, 3) for support in supports]
        self.assertEqual(
            simulate_wrapping_candidates(supports, taps, memory=3),
            expected,
        )


if __name__ == "__main__":
    unittest.main()

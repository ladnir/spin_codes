import unittest

import numpy as np

from ec_region_attack import (
    boundary_response_labels,
    constraint_mask,
    decode_region_key,
    encode_region_keys,
    gf2_kernel_dependencies,
    region_boundaries,
    run_attack,
    sample_exact_rows,
)


def naive_output(support, taps):
    values = [0] * len(taps)
    for coordinate in support:
        values[coordinate] ^= 1
    for source, tap_mask in enumerate(taps):
        if not values[source]:
            continue
        pending = tap_mask
        while pending:
            bit = pending & -pending
            target = source + bit.bit_length()
            if target < len(values):
                values[target] ^= 1
            pending ^= bit
    return values


class ECRegionAttackTests(unittest.TestCase):
    def test_exact_rows_are_sorted_and_distinct(self):
        rows = sample_exact_rows(
            row_count=200,
            domain_size=17,
            row_weight=5,
            seed=9,
        )
        self.assertEqual(rows.shape, (200, 5))
        self.assertTrue(np.all(rows[:, 1:] > rows[:, :-1]))

    def test_region_key_round_trip(self):
        rows = np.array(
            [[0, 3, 4, 8, 11], [1, 2, 6, 7, 10]],
            dtype=np.uint32,
        )
        keys = encode_region_keys(rows, region_count=4, domain_size=12)
        self.assertEqual(decode_region_key(int(keys[0]), 4, 5), (0, 1, 1, 2, 3))
        self.assertEqual(decode_region_key(int(keys[1]), 4, 5), (0, 0, 2, 2, 3))

    def test_reverse_boundary_labels_match_forward_simulation(self):
        taps = [0b101, 0b110, 0b100, 0b111, 0b101, 0b100,
                0b110, 0b101, 0b111, 0b100, 0b110, 0b101, 0b100]
        boundaries = region_boundaries(len(taps), 3)
        labels = boundary_response_labels(taps, memory=3, boundaries=boundaries)
        for coordinate in range(len(taps)):
            output = naive_output((coordinate,), taps)
            expected = 0
            for region, end in enumerate(boundaries[1:-1]):
                for offset, value in enumerate(output[end - 3 : end]):
                    expected |= value << (region * 3 + offset)
            self.assertEqual(labels[coordinate], expected)

    def test_kernel_dependencies_xor_to_zero(self):
        columns = [0b001, 0b010, 0b011, 0b100, 0b101]
        dependencies = gf2_kernel_dependencies(columns)
        self.assertEqual(len(dependencies), 2)
        for dependency in dependencies:
            value = 0
            for index, column in enumerate(columns):
                if dependency & (1 << index):
                    value ^= column
            self.assertEqual(value, 0)

    def test_constraints_only_precede_unselected_gaps(self):
        mask = constraint_mask((1, 2, 4), region_count=6, memory=3)
        self.assertEqual(mask, (0b111 << 6) | (0b111 << 12))

    def test_small_end_to_end_witness_is_confined(self):
        result = run_attack(
            k=1 << 10,
            row_weight=5,
            memory=3,
            region_count=4,
            seed=3,
            max_buckets=2,
            kernel_samples=2,
            late_row_candidates=8,
        )
        witness = result["best_witness"]
        self.assertIsNotNone(witness)
        self.assertEqual(witness["outside_region_weight"], 0)
        self.assertEqual(len(witness["codeword_sha256"]), 64)
        self.assertEqual(
            witness["codeword_weight"],
            witness["message_weight"] + witness["parity_weight"],
        )


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3

import unittest

import numpy as np

from streaming_ec_distance_ablation import (
    COLLISION_FREE_VARIANTS,
    VARIANTS,
    convolve_generator,
    convolution_taps,
    null_vector,
    run_variant,
    structured_column_labels,
)


class StreamingEcDistanceAblationTest(unittest.TestCase):
    def test_null_vector(self) -> None:
        matrix = np.array([[1, 2, 3], [2, 4, 6]], dtype=np.int64)
        vector = null_vector(matrix, 127)
        self.assertIsNotNone(vector)
        self.assertTrue(np.all((matrix @ vector) % 127 == 0))
        self.assertIsNone(null_vector(np.eye(3, dtype=np.int64), 127))

    def test_periodic_taps_repeat(self) -> None:
        taps = convolution_taps(
            code_size=19,
            memory=3,
            period=5,
            prime=127,
            seed=4,
            mode="periodic",
        )
        self.assertTrue(np.array_equal(taps[:5], taps[5:10]))
        self.assertTrue(np.array_equal(taps[:5], taps[10:15]))

    def test_zero_taps_leave_generator_unchanged(self) -> None:
        matrix = np.arange(30, dtype=np.int64).reshape(5, 6) % 127
        taps = np.zeros((6, 4), dtype=np.int64)
        self.assertTrue(np.array_equal(convolve_generator(matrix, taps, 127), matrix))

    def test_structured_labels_are_nonzero_and_deterministic(self) -> None:
        for mode in (
            "column_splitmix",
            "column_affine",
            "column_sign",
            "column_unit",
        ):
            first = structured_column_labels(
                mode=mode, size=100, prime=127, seed=9
            )
            second = structured_column_labels(
                mode=mode, size=100, prime=127, seed=9
            )
            self.assertTrue(np.array_equal(first, second))
            self.assertTrue(np.all(first != 0))
            if mode == "column_sign":
                self.assertEqual(set(first), {1, 126})

    def test_all_variants_run(self) -> None:
        for variant in VARIANTS:
            result = run_variant(
                variant=variant,
                left_degree=6,
                right_degree=3,
                region_size=3,
                memory=2,
                period=5,
                prime=127,
                seed=2,
            )
            self.assertEqual(result.rank, 9)
            self.assertGreater(result.row_weight, 0)
            self.assertGreater(result.pair_weight, 0)

    def test_collision_free_variant_runs_when_geometry_allows_it(self) -> None:
        for variant in COLLISION_FREE_VARIANTS:
            result = run_variant(
                variant=variant,
                left_degree=4,
                right_degree=2,
                region_size=5,
                memory=2,
                period=4,
                prime=127,
                seed=2,
            )
            self.assertEqual(result.rank, 10)

    def test_edge_signs_remove_observed_low_degree_kernel(self) -> None:
        variants = {variant.name: variant for variant in COLLISION_FREE_VARIANTS}
        common = dict(
            left_degree=4,
            right_degree=2,
            region_size=10,
            memory=4,
            period=8,
            prime=3,
            seed=2,
        )
        base_signs = run_variant(
            variant=variants["collision_free_heuristic"], **common
        )
        edge_signs = run_variant(
            variant=variants["collision_free_edge_signs"], **common
        )
        self.assertEqual(base_signs.rank, 19)
        self.assertEqual(edge_signs.rank, 20)


if __name__ == "__main__":
    unittest.main()

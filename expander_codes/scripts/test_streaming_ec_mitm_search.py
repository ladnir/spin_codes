#!/usr/bin/env python3

import unittest

import numpy as np

from streaming_ec_mitm_search import (
    combinations_with_zero,
    mitm_search,
    projection_keys,
)


class StreamingEcMitmSearchTest(unittest.TestCase):
    def test_combinations_include_zero_once(self) -> None:
        messages = combinations_with_zero(3, 3)
        self.assertEqual(messages.shape, (27, 3))
        self.assertEqual(int(np.sum(np.all(messages == 0, axis=1))), 1)

    def test_projection_negation(self) -> None:
        words = np.array([[1, 2, 0], [2, 1, 1]], dtype=np.uint8)
        coordinates = np.array([0, 1])
        positive = projection_keys(words, coordinates, 3, negate=False)
        negative = projection_keys((-words.astype(int)) % 3, coordinates, 3, negate=False)
        direct = projection_keys(words, coordinates, 3, negate=True)
        self.assertTrue(np.array_equal(negative, direct))
        self.assertFalse(np.array_equal(positive, direct))

    def test_search_finds_weight_one_word(self) -> None:
        generator = np.array(
            [
                [1, 0, 1, 1, 0, 0],
                [0, 1, 1, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            dtype=np.int64,
        )
        best, _, _, _ = mitm_search(
            generator,
            prime=3,
            zero_coordinates=2,
            trials=10,
            seed=1,
        )
        self.assertEqual(best, 1)


if __name__ == "__main__":
    unittest.main()

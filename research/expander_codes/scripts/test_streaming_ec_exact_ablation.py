#!/usr/bin/env python3

import unittest

import numpy as np

from streaming_ec_exact_ablation import all_messages, exact_minimum_distance


class StreamingEcExactAblationTest(unittest.TestCase):
    def test_repetition_code_distance(self) -> None:
        generator = np.array([[1, 1, 1], [0, 1, 2]], dtype=np.int64)
        messages = all_messages(2, 3)
        self.assertEqual(exact_minimum_distance(generator, messages, 3), 2)

    def test_message_enumeration(self) -> None:
        messages = all_messages(3, 3)
        self.assertEqual(messages.shape, (26, 3))
        self.assertFalse(np.any(np.all(messages == 0, axis=1)))


if __name__ == "__main__":
    unittest.main()

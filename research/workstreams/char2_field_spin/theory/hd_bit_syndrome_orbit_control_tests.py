"""Cheap provenance and normalization tests; never invokes a provider."""
from fractions import Fraction as F
from hashlib import sha256
from math import log2
from unittest import TestCase, main

import hd_bit_syndrome_orbit_control as gate
from hd_bit_syndrome_gate_run import prefix


class OrbitControlTests(TestCase):
    def test_pinned_source_bytes(self):
        for name, expected in gate.SOURCE_PINS.items():
            with self.subTest(name=name):
                self.assertEqual(sha256((gate.DIRECTORY / name).read_bytes()).hexdigest(), expected)

    def test_original_and_alternative_prefix(self):
        self.assertAlmostEqual(gate.global_prefix(F(691, 700), 4), prefix(4), places=10)
        change = -gate.LIMIT * log2(float(F(497, 500) / F(691, 700)))
        self.assertAlmostEqual(gate.global_prefix(F(497, 500), 0), prefix(0) + change, places=9)
        self.assertLess(gate.global_prefix(F(497, 500), 0), 0)

    def test_invalid_marker(self):
        for mark, x in ((F(0), 0), (F(1), 0), (F(1, 2), -1), (F(1, 2), float('nan'))):
            with self.assertRaises(ValueError):
                gate.global_prefix(mark, x)


if __name__ == '__main__':
    main()

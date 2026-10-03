from itertools import product
from math import comb
import unittest

import numpy as np

import maps24
import screen24


class ExactTests(unittest.TestCase):
    def test_literal_candidate(self):
        data = maps24.make_maps()
        record = data['record']
        self.assertEqual(record['expansion_rank'], 24)
        self.assertEqual(record['feedback_rank'], 24)
        self.assertEqual(record['packet_restriction_rank_counts'],
                         {1: {8: 8}, 2: {16: 28}, 3: {24: 56}, 4: {24: 70}})
        self.assertEqual(record['CA_columns'], [0]*24)
        for h in range(8):
            for i in range(24):
                a, b, c = ((1 << i) >> (8*k) & 255 for k in range(3))
                expected = a ^ maps24.multiply(h, b) ^ maps24.multiply(maps24.multiply(h, h), c)
                self.assertEqual((data['rows'][i] >> (8*h)) & 255, expected)

    def test_chunked_profiles(self):
        data = maps24.make_maps(packet_bits=2, modulus=7, windows=4)
        counts = maps24.census(data, chunk_bits=3)
        for state in range(64):
            image = maps24.apply(data['rows'], state)
            character = sum(((state & col).bit_count() & 1) << i for i, col in enumerate(data['columns']))
            for word, indices in ((image, counts['expansion_indices']), (character, counts['character_indices'])):
                weights = [(word >> (2*h) & 3).bit_count() for h in range(4)]
                np.testing.assert_array_equal(data['profiles'][indices[state]], np.bincount(weights, minlength=3))

    def test_bounded_walsh(self):
        for bits in (1, 3, 5, 8):
            values = np.arange(1 << bits, dtype=float)**2+.37
            expected = screen24.prior.walsh(values)
            for block_bits in (1, 3, 9):
                actual = screen24.walsh_in_place(values.copy(), block_bits=block_bits)
                np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=1e-13)

    def test_compressed_birth_closure(self):
        data = maps24.make_maps(packet_bits=2, modulus=7, windows=4)
        counts = maps24.census(data, chunk_bits=4)
        z = .73
        weighted, emission = np.zeros((64, 5)), np.zeros((64, 5))
        images = [maps24.apply(data['rows'], a) for a in range(64)]
        for x in range(256):
            j = sum(bool(x >> (2*h) & 3) for h in range(4))
            den = comb(4, j)*3**j
            weighted[maps24.apply(data['columns'], x), j] += z**x.bit_count()/den
            for a in range(64):
                emission[a, j] += z**(x ^ images[a]).bit_count()/den
        direct = screen24.prior.operators(weighted, emission)
        actual, _ = screen24.local_operators(data, counts, z)
        np.testing.assert_allclose(actual, direct, atol=3e-15, rtol=3e-12)
        small, _ = screen24.local_operators(data, counts, z, only_q1=True)
        np.testing.assert_allclose(small, direct[:2, :3, :3], atol=3e-15, rtol=3e-12)


if __name__ == '__main__':
    unittest.main()

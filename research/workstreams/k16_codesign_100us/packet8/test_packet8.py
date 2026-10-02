"""Independent finite checks for the isolated eight-bit packet screen."""
import itertools
import unittest

import numpy as np
from flint import arb, ctx

import screen_packet8 as base
import screen_tail as tail


class PacketEightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 256
        cls.data, cls.record = base.construction()

    def test_permutation_spectrum_and_injective_packets(self):
        permutation = self.record['coordinate_permutation']
        self.assertEqual(sorted(permutation), list(range(64)))
        self.assertEqual([i for i in range(64) if permutation[i] != i],
                         [7, 15, 23, 31, 39, 47, 55, 63])
        self.assertEqual(self.data['packet_ranks'], [8]*8)
        self.assertEqual(sum(self.data['spectrum'].values()), 65536)
        self.assertEqual(self.data['spectrum'][0], 1)
        self.assertTrue(all((a & b).bit_count() % 2 == 0
                            for a in self.data['rows'] for b in self.data['rows']))

    def test_single_packet_polynomial_direct_enumeration(self):
        # Four output coordinates, two-bit packets, three state coordinates.
        small = base.prepare((0b1111, 0b1010, 0b1100), packet_bits=2)
        z = base.aq(base.Q(2, 3))
        for image in base.images_from_rows(small['rows']):
            if not image:
                continue
            profile = [0]*3
            for p in (0, 2):
                profile[((image >> p) & 3).bit_count()] += 1
            expected = sum((z**(image ^ (label << p)).bit_count()
                            for p in (0, 2) for label in (1, 2, 3)), arb(0))/6
            self.assertLess(abs(expected-base.emission(small, profile, z, 1)), arb(2)**-240)

    def test_walsh_matches_literal_characters(self):
        values = np.array([2, -1, 3, 0, 7, -2, 1, 4], dtype=np.int64)
        expected = [sum(int(value)*(-1)**((a & x).bit_count())
                        for x, value in enumerate(values)) for a in range(8)]
        self.assertEqual(tail.walsh(values).tolist(), expected)

    def test_q1_independent_birth_census_matches_full_fourier(self):
        full = tail.prepare_full(self.data)
        # Disable cap selection for this equality check by choosing zero caps only.
        original = tail.density.thresholds
        tail.density.thresholds = lambda values: []
        try:
            direct = base.local_operators(self.data, '.00512')
            fourier = tail.full_local(full, '.00512')
        finally:
            tail.density.thresholds = original
        for occupied in (0, 1):
            for row in range(direct[occupied].nrows()):
                for column in range(direct[occupied].ncols()):
                    self.assertLess(abs(direct[occupied][row, column]-fourier[occupied][row, column]), arb(2)**-230)

    def test_full_birth_mass_and_nonnegative_local_operators(self):
        full = tail.prepare_full(self.data)
        original = tail.density.thresholds
        tail.density.thresholds = lambda values: []
        try:
            local = tail.full_local(full, '.0256')
        finally:
            tail.density.thresholds = original
        z = (-base.aq(base.Q('.0256'))).exp()
        for j, matrix in enumerate(local):
            self.assertTrue(all(matrix[r, c] >= 0 for r in range(matrix.nrows()) for c in range(matrix.ncols())))
            total = sum((matrix[0, c] for c in range(matrix.ncols())), arb(0))
            expected = (((1+z)**8-1)/255)**j
            self.assertLess(abs(total-expected), arb(2)**-230)


if __name__ == '__main__':
    unittest.main()

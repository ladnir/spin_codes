import itertools
import math
import unittest
from unittest import mock
from fractions import Fraction

import syndrome_density_v1 as density
import test_occupation_refresh_v1 as previous


class DensityTransferTests(previous.RefreshTests):
    def setUp(self):
        self.model = density.Epochs(4, 2, {2: 2, 4: 1}, [1, 0, 2, 0, 1])

    def test_dense_selected_boxes_cover_all_integer_types(self):
        with mock.patch.object(previous.refresh, 'Epochs', density.Epochs):
            super().test_dense_selected_boxes_cover_all_integer_types()


class CharacterTests(unittest.TestCase):
    def test_krawtchouk_against_signed_subsets(self):
        for weight in range(9):
            row = density.krawtchouk_row(8, weight)
            for j in range(9):
                exact = sum((-1)**sum(i < weight for i in positions)
                            for positions in itertools.combinations(range(8), j))
                self.assertEqual(row[j], exact)

    def test_every_actual_syndrome_probability(self):
        generators = [255, 170, 204, 240]
        model = density.Epochs(8, 4, {4: 14, 8: 1}, [1, 0, 0, 0, 14, 0, 0, 0, 1])
        for j in range(9):
            counts = [0]*16
            for positions in itertools.combinations(range(8), j):
                u = sum(1 << p for p in positions)
                syndrome = sum(((u & g).bit_count() % 2) << bit for bit, g in enumerate(generators))
                counts[syndrome] += 1
            law = model.activation[j]
            upper = Fraction(law['numerator'], law['denominator'])
            self.assertGreaterEqual(upper, max(Fraction(n, math.comb(8, j)) for n in counts[1:]))
            if law['use_density']:
                mass = Fraction(law['mass_numerator'], law['denominator'])
                self.assertGreaterEqual(mass, 1-Fraction(counts[0], math.comb(8, j)))
                self.assertGreaterEqual(mass/14, max(Fraction(n, math.comb(8, j)) for n in counts[1:]))


if __name__ == '__main__':
    unittest.main()

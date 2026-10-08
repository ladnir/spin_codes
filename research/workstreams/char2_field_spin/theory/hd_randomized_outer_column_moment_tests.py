"""Integer checks for the diagnostic coupled-column spectra."""

from collections import Counter
from math import comb
import unittest

from hd_blockwise_shared_maps import multiplication_table, rank
from hd_randomized_outer_column_moment import (
    aggregate_enumerator, coupled_matrix, moment_log, pair_projective_weights,
    refine_columns,
)
from hd_randomized_outer_state_profile import local_profiles


class CoupledColumnTests(unittest.TestCase):
    def test_invertible_full_weight_columns(self):
        matrix = coupled_matrix()
        self.assertEqual(rank(matrix), 64)
        self.assertTrue(all(sum(x != 0 for x in column) == 64
                            for column in zip(*matrix)))

    def test_pair_histogram_exhaustive_base_field(self):
        table = multiplication_table()
        columns = tuple(zip(*coupled_matrix()))
        for i, j in ((0, 1), (0, 8), (0, 63), (9, 21)):
            first, second = columns[i], columns[j]
            exact = Counter(sum(a ^ table[t][b] != 0 for a, b in zip(first, second))
                            for t in range(256))
            exact[sum(b != 0 for b in second)] += 1
            self.assertEqual(pair_projective_weights(first, second, 256), exact)
            extension = exact.copy()
            extension[64] += 65536-256
            self.assertEqual(pair_projective_weights(first, second), extension)

    def test_integer_normalization_and_distance(self):
        for j in (0, 1, 2):
            counts = aggregate_enumerator(j)
            self.assertEqual(sum(counts), comb(64, j)*65536**j)
            self.assertEqual(counts[0], comb(64, j))
            self.assertAlmostEqual(moment_log(j, 0), 0, places=11)
        self.assertEqual(min(w for w in range(1, 65) if aggregate_enumerator(1)[w]), 64)
        self.assertEqual(min(w for w in range(1, 65) if aggregate_enumerator(2)[w]), 60)

    def test_refinement_cannot_increase_an_entry(self):
        old = local_profiles(-1/64)
        new = refine_columns(old, -1/64)
        for j in range(65):
            for a in range(3):
                for b in range(3):
                    self.assertLessEqual(new[j][a][b], old[j][a][b])
            if j not in (1, 2):
                self.assertEqual(new[j], old[j])

    def test_invalid_fields(self):
        for field in (255, 512, 65535, 65536.0):
            with self.assertRaises(ValueError):
                aggregate_enumerator(1, field)


if __name__ == '__main__':
    unittest.main()

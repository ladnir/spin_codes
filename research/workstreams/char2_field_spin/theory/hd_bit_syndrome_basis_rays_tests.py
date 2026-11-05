"""Tiny exact ray tests; no complete row-set preparation or contraction."""
from collections import Counter
from itertools import combinations
import unittest
from unittest.mock import patch

from gf8_packet_maps import MUL, matvec, rank
from hd_bit_syndrome_outer import affine_matrix, matrix
import hd_bit_syndrome_basis_rays as rays


class BasisRayTests(unittest.TestCase):
    SOURCES = ((0, 1, 2, 4, 8), (0, 1, 3, 4, 8), (3, 5, 6, 7, 15))

    def test_exact_five_rays_and_actual_matrix(self):
        h, v = matrix(), affine_matrix()
        for source in self.SOURCES:
            data = rays.basis_rays(source)
            self.assertEqual(data['ray_count'], 5)
            self.assertEqual(data['nonzero_vector_count'], 1275)
            self.assertEqual(tuple(row['isolated_source'] for row in data['rays']), source)
            for row in data['rays']:
                self.assertEqual(matvec(h, row['input_values']), row['output_values'])
                self.assertEqual(matvec(v, row['affine_coefficients']), row['correction_values'])
                self.assertEqual(tuple(row['correction_values'][i] for i in source),
                                 tuple(int(i != row['isolated_source']) for i in source))
                self.assertEqual(row['outside_weight'], 4)
                self.assertEqual(row['count'], 255)
                self.assertFalse(row['input_mask'] & ~data['source_mask'])
                self.assertEqual(row['output_mask'], row['outside_mask'] | row['inside_output_mask'])
            self.assertFalse(data['metadata']['complete_outer_kernel'])

    def test_independent_full_H_shortening_counts_all_330_supports(self):
        source = self.SOURCES[0]
        h = matrix()
        outside = tuple(i for i in range(16) if i not in source)
        actual = {}
        # No outside word of weight1..3 exists: any such zero set would
        # contain at least8 outside points, hence a hyperplane missing A.
        for support in combinations(outside, 4):
            zero_positions = tuple(i for i in outside if i not in support)
            constraints = tuple(tuple(h[i][j] for j in source) for i in zero_positions)
            count = 256**(5-rank(constraints))-1
            if count:
                actual[sum(1 << i for i in support)] = count
        self.assertEqual(actual, rays.basis_rays(source)['outside_weight4_counts'])

    def test_every_scalar_multiple_and_positive_terms(self):
        data = rays.basis_rays(self.SOURCES[0])
        h = matrix()
        counts = Counter()
        normalized_outputs = set()
        for row in data['rays']:
            for scale in range(1, 256):
                x = tuple(MUL[scale][value] for value in row['input_values'])
                y = matvec(h, x)
                self.assertEqual(y, tuple(MUL[scale][value] for value in row['output_values']))
                key = (sum(bool(value) for value in x), sum(1 << i for i, value in enumerate(y) if value))
                counts[key] += 1
                normalized_outputs.add(y)
        self.assertEqual(len(normalized_outputs), 1275)
        self.assertEqual(dict(counts), rays.exact_terms(self.SOURCES[0]))

    def test_five_stars_not_all_weight_four_supports(self):
        for source in self.SOURCES:
            data = rays.basis_rays(source)
            masks = tuple(row['outside_mask'] for row in data['rays'])
            self.assertTrue(all((left & right).bit_count() == 1 for left, right in combinations(masks, 2)))
            degrees = [sum(mask >> i & 1 for mask in masks) for i in range(16)
                       if data['outside_positions_mask'] >> i & 1]
            self.assertEqual(sorted(degrees), [0]+[2]*10)
            self.assertEqual(len(masks), 5)

    def test_exhaustiveness_is_checked_not_assumed(self):
        source = self.SOURCES[0]
        with patch.object(rays, 'outside_image_support_counts', return_value={0: 1}):
            with self.assertRaises(AssertionError):
                rays.basis_rays(source)

    def test_source_scope_guards(self):
        for source in ((), (0, 1, 2, 4), (0, 1, 2, 3, 4), (0, 1, 2, 4, 16),
                       (1, 0, 2, 4, 8), (0, 1, 2, 2, 8), (False, 1, 2, 4, 8)):
            with self.assertRaises(ValueError):
                rays.basis_rays(source)


if __name__ == '__main__':
    unittest.main()

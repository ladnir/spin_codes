"""Tiny exact rank4 fiber tests and the900-entry no-hyperplane certificate."""
from collections import Counter
from itertools import combinations
import unittest
from unittest.mock import patch

from gf8_packet_maps import INV, MUL, matvec, rank, xor_sum
from hd_bit_syndrome_cross_fibers import cross_fiber
from hd_bit_syndrome_outer import BASIS, OFFSET, matrix
import hd_bit_syndrome_rank4_rays as rays


class Rank4RayTests(unittest.TestCase):
    SOURCES = ((0, 1, 2, 3, 4), (0, 1, 2, 3, 8), (2, 3, 4, 5, 9))

    def test_universal_hyperplane_certificate_independent_inverse_formula(self):
        pairings = rays.hyperplane_pairings()
        w = (OFFSET,)+BASIS
        denominator = 1 ^ xor_sum(MUL[x][x] for x in w)
        direct = []
        for nmask in range(2, 32):
            n = tuple(nmask >> j & 1 for j in range(5))
            row = []
            for cmask in range(2, 32):
                c = tuple(cmask >> j & 1 for j in range(5))
                beta = MUL[matvec((w,), c)[0]][INV[denominator]]
                inverse = tuple(x ^ MUL[beta][wi] for x, wi in zip(c, w))
                row.append(matvec((n,), inverse)[0])
            direct.append(tuple(row))
        self.assertEqual(pairings, tuple(direct))
        self.assertEqual(sum(value != 0 for row in pairings for value in row), 900)

    def test_selected_rank4_weight4_branches_are_empty(self):
        for source in self.SOURCES:
            data = rays.rank4_rays(source)
            self.assertEqual(data['rays'], ())
            self.assertEqual(data['outside_weight4_counts'], {})
            self.assertEqual(data['nonzero_vector_count'], 0)
            self.assertEqual(rays.exact_terms(source), {})
            self.assertEqual(data['metadata']['normalized_fiber_cardinality'], 256)
            self.assertFalse(data['metadata']['complete_outer_kernel'])

    def test_independent_full_H_outside_four_supports(self):
        source = self.SOURCES[0]
        h = matrix()
        outside = tuple(i for i in range(16) if i not in source)
        for support in combinations(outside, 4):
            constraints = tuple(tuple(h[i][j] for j in source) for i in outside if i not in support)
            # The only solutions are the256 zero-outside quiet vectors.
            self.assertEqual(256**(5-rank(constraints)), 256)

    def test_breakpoint_counts_against_all256_parameters(self):
        h = matrix()
        for source in self.SOURCES:
            fiber = cross_fiber(source)
            kernel = fiber['input_kernel_basis'][0]
            direction = tuple(kernel[source.index(i)] if i in source else 0 for i in range(16))
            x = tuple(int(i == source[0]) for i in range(16))
            y = matvec(h, x)
            self.assertEqual(matvec(h, direction), direction)
            direct = Counter()
            for beta in range(256):
                moved_x = tuple(a ^ MUL[beta][k] for a, k in zip(x, direction))
                moved_y = matvec(h, moved_x)
                self.assertEqual(moved_y, tuple(a ^ MUL[beta][k] for a, k in zip(y, direction)))
                im = sum(1 << i for i, value in enumerate(moved_x) if value)
                om = sum(1 << i for i, value in enumerate(moved_y) if value)
                direct[im, om] += 1
            self.assertEqual(rays.line_support_counts(x, y, direction), dict(direct))
            self.assertEqual(sum(direct.values()), 256)

    def test_exhaustiveness_check_rejects_an_unexpected_branch(self):
        with patch.object(rays, 'outside_image_support_counts', return_value={15 << 4: 255}):
            with self.assertRaises(AssertionError):
                rays.rank4_rays(self.SOURCES[0])

    def test_guards(self):
        for source in ((), (0, 1, 2, 4, 8), (0, 1, 2, 3), (1, 0, 2, 3, 4),
                       (0, 1, 2, 3, 16), (False, 1, 2, 3, 4)):
            with self.assertRaises(ValueError):
                rays.rank4_rays(source)
        for x, y, k in (((), (), ()), ((0,)*16, (0,)*16, (0,)*16),
                        ((256,)+(0,)*15, (0,)*16, (1,)+(0,)*15)):
            with self.assertRaises(ValueError):
                rays.line_support_counts(x, y, k)


if __name__ == '__main__':
    unittest.main()

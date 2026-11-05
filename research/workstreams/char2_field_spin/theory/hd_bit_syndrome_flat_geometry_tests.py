"""Small exact checks; no full rank table, provider, or distance gate."""
from collections import Counter
from itertools import combinations
import unittest

from gf8_packet_maps import MUL, rank
from hd_bit_syndrome_cross_fibers import cross_fiber, information_sets
from hd_bit_syndrome_outer import affine_rank, matrix
import hd_bit_syndrome_flat_geometry as flat


def mask(labels):
    return sum(1 << i for i in labels)


class FlatGeometryTests(unittest.TestCase):
    def test_flat_enumeration_and_small_exhaustive_hulls(self):
        flats = flat.affine_flats()
        self.assertEqual(len(flats), 308)
        self.assertEqual(Counter(value.bit_count() for value in flats),
                         {0: 1, 1: 16, 2: 120, 4: 140, 8: 30, 16: 1})
        self.assertTrue(all(flat.flat_hull_mask(value) == value for value in flats))
        for bits in (1, 2, 3):
            all_hulls = {flat.flat_hull_mask(value, bits) for value in range(1 << (1 << bits))}
            self.assertEqual(all_hulls, set(flat.affine_flats(bits)))
        self.assertEqual(flat.flat_hull_mask(0), 0)
        self.assertEqual(flat.flat_hull_mask(0b0111), 0b1111)

    def test_bilinear_rank_including_overlap_and_empty(self):
        h = matrix()
        cases = ((0, 0xffff), (0x0007, 0x0780), (0x001f, 0xf000),
                 (0x000f, 0x003c), (0xffff, 0xffff), (0x5555, 0xaaaa),
                 (0x033f, 0x0f03), (0x0001, 0x0001))
        for left, right in cases:
            rows = tuple(i for i in range(16) if left >> i & 1)
            columns = tuple(i for i in range(16) if right >> i & 1)
            expected = rank(tuple(tuple(h[i][j] ^ int(i == j) for j in columns) for i in rows))
            self.assertEqual(flat.rank_cross_masks(left, right), expected)
            self.assertEqual(flat.rank_cross_masks(right, left), expected)
            if not left & right:
                self.assertEqual(expected, rank(tuple(tuple(h[i][j] for j in columns) for i in rows)))
        self.assertEqual(flat.rank_cross_masks(0xffff, 0xffff), 5)
        self.assertEqual(rank(h), 16)  # Overlapping bilinear rank is NOT H's rank.

    def test_hull_pair_cache_and_actual_cross_rank(self):
        flat._rank_flat_pair.cache_clear()
        first = flat.rank_cross_masks(0x0007, 0x0780)
        before = flat._rank_flat_pair.cache_info()
        self.assertEqual(flat.rank_cross_masks(0x0780, 0x000f), first)
        after = flat._rank_flat_pair.cache_info()
        self.assertEqual(after.misses, before.misses)
        self.assertEqual(after.hits, before.hits+1)
        for source in ((0,), (0, 3), (0, 1, 2), tuple(range(7)), tuple(range(6, 16))):
            a = mask(source)
            self.assertEqual(flat.rank_cross_masks(0xffff ^ a, a), cross_fiber(source)['cross_rank'])

    def test_zero_flat_counts_and_exact_containment_mass(self):
        self.assertEqual(flat.correction_zero_set_counts(0), {0xffff: 1})
        for source in (0x0001, 0x0009, 0x0007):
            counts = flat.correction_zero_set_counts(source)
            self.assertTrue(all(value > 0 for value in counts.values()))
            self.assertTrue(set(counts) <= set(flat.affine_flats()))
            self.assertEqual(sum(counts.values()), 256 ** affine_rank(source))
            self.assertEqual(counts[0xffff], 1)
            for required in (0, 1, 3, 15, 255, 0xffff):
                expected = 256 ** (affine_rank(source)-flat.rank_cross_masks(required, source))
                actual = sum(value for zeros, value in counts.items() if not required & ~zeros)
                self.assertEqual(actual, expected)

    def test_projective_pair_enumeration_matches_actual_outside_counts(self):
        h = matrix()
        for i, j in ((0, 1), (2, 11)):
            a = (1 << i) | (1 << j)
            outside = tuple(k for k in range(16) if not a >> k & 1)
            # Every nonzero vector in F256^2 belongs to one of257 lines.
            # C has rank2, so each outside word has one preimage.
            direct = Counter({0: 1})
            for ratio in range(256):
                support = mask(k for k in outside if h[k][i] ^ MUL[h[k][j]][ratio])
                direct[support] += 255
            direct[mask(k for k in outside if h[k][j])] += 255
            self.assertEqual(dict(direct), flat.outside_image_support_counts(a))
        for i in (0, 7, 15):
            self.assertEqual(flat.outside_image_support_counts(1 << i),
                             {0: 1, 0xffff ^ (1 << i): 255})

    def test_outside_supported_counts_without_support_enumeration(self):
        for source in (0x0007, 0x002b, 0x007f):
            outside = 0xffff ^ source
            counts = flat.outside_image_support_counts(source)
            dimension = affine_rank(source)
            self.assertEqual(counts[0], 1)
            self.assertEqual(sum(counts.values()), 256 ** dimension)
            allowed = (0, outside, outside & 0xaaaa, next(k for k in counts if k))
            for support in allowed:
                actual = sum(value for target, value in counts.items() if not target & ~support)
                expected = 256 ** (dimension-flat.rank_cross_masks(outside & ~support, source))
                self.assertEqual(actual, expected)

    def test_complement_cross_bases_are_actual_information_sets(self):
        # One high-arity case, not all-A preparation. This independently
        # checks the cheap dual-basis interface against the existing fiber.
        source = tuple(range(2, 16))
        a, outside = mask(source), 3
        dimension = flat.rank_cross_masks(outside, a)
        complements = {tuple(i for i in source if i not in basis)
                       for basis in combinations(source, dimension)
                       if flat.rank_cross_masks(outside, mask(basis)) == dimension}
        self.assertEqual(complements, set(information_sets(cross_fiber(source))))

    def test_guards(self):
        flat.flat_hull_mask(1, 4)
        flat.affine_flats(bits=4)
        with self.assertRaises(ValueError):
            flat.flat_hull_mask(True, 4)
        with self.assertRaises(ValueError):
            flat.affine_flats(bits=4.0)
        for value in (-1, 65536, True):
            with self.assertRaises(ValueError):
                flat.flat_hull_mask(value)
            with self.assertRaises(ValueError):
                flat.rank_cross_masks(value, 1)
        for value in (0, 255, 65535):
            with self.assertRaises(ValueError):
                flat.outside_image_support_counts(value)
        with self.assertRaises(ValueError):
            flat.affine_flats(5)


if __name__ == '__main__':
    unittest.main()

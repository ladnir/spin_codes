"""Small exact geometry tests; no provider or full sparse preparation."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
import unittest
from unittest.mock import patch

import numpy as np

from gf8_packet_maps import MUL
from hd_bit_syndrome_arity7_gate import actual_image_caps
from hd_bit_syndrome_arity9_gate import gap_geometry
from hd_bit_syndrome_cross_fibers import cross_fiber, information_sets
from hd_bit_syndrome_kernel import information_cover
from hd_bit_syndrome_outer import affine_matrix, matrix
import hd_bit_syndrome_sparse_geometry as geometry


def mask(indices):
    return sum(1 << i for i in indices)


class SparseGeometryTests(unittest.TestCase):
    def test_toy_zero_flat_restriction_and_constant_fiber(self):
        # Binary generator over GF4, so multiplication by entries is0/1.
        zero_sets, words = Counter(), []
        for x, y in product(range(4), repeat=2):
            word = (x, y, x ^ y, x)
            words.append(word)
            zero_sets[mask(i for i, value in enumerate(word) if not value)] += 1
        for outside, rank in ((0, 0), (1, 1), (9, 1), (6, 2), (15, 2)):
            distinct = {tuple(word[i] for i in range(4) if outside >> i & 1) for word in words}
            positions = tuple(i for i in range(4) if outside >> i & 1)
            direct = Counter(mask(positions[i] for i, value in enumerate(word) if value)
                             for word in distinct)
            result = geometry.aggregate_outside_supports(outside, zero_sets, 2, rank, q=4, size=4)
            self.assertEqual(result, dict(direct))
            self.assertEqual(result[0], 1)
            self.assertEqual(sum(result.values()), 4**rank)

    def test_actual_pair_image_by_projective_aes_enumeration(self):
        h = matrix()
        for source in ((0, 1), (3, 14)):
            outside = tuple(i for i in range(16) if i not in source)
            direct = Counter({0: 1})
            for ratio in range(256):
                support = mask(i for i in outside if h[i][source[0]] ^ MUL[h[i][source[1]]][ratio])
                direct[support] += 255
            direct[mask(i for i in outside if h[i][source[1]])] += 255
            result = geometry.actual_image_data(source)
            self.assertEqual(result['support_counts'], dict(direct))
            self.assertEqual(result['cross_rank'], 2)
            self.assertEqual(result['fiber_dimension'], 0)
            self.assertEqual(result['restriction_fiber_size'], 1)

    def test_rank_four_and_five_a7_against_old_actual_elimination(self):
        for source in (tuple(range(7)), (0, 1, 2, 3, 4, 5, 8)):
            old_caps, rank, old_counts = actual_image_caps(source)
            new = geometry.actual_image_data(source)
            outside = tuple(i for i in range(16) if i not in source)
            expected = {mask(outside[i] for i in range(9) if local >> i & 1): count
                        for local, count in enumerate(old_counts) if count}
            self.assertEqual(new['support_counts'], expected)
            self.assertEqual(new['full_caps'], old_caps)
            self.assertEqual(new['cross_rank'], rank)
            self.assertEqual(new['nonquiet_caps'], (0,)+old_caps[1:])
            self.assertEqual(sum(new['support_counts'].values()), 256**rank)

    def test_quiet_coefficients_toy_field_direct(self):
        v = affine_matrix()
        for source in ((0, 1, 2, 3), (0, 1, 2, 3, 4)):
            direct = Counter()
            for values in product(range(4), repeat=len(source)):
                syndrome = [0]*5
                for i, value in zip(source, values):
                    for j, coefficient in enumerate(v[i]):
                        if coefficient:
                            syndrome[j] ^= value
                if not any(syndrome):
                    direct[mask(i for i, value in zip(source, values) if value)] += 1
            result = geometry.quiet_data(source, q=4)
            self.assertEqual(result['support_counts'], dict(direct))
            self.assertEqual(result['caps'], geometry.caps_by_weight(direct, len(source)))
            self.assertEqual(result['caps'][0], 1)

    def test_partial_each_arity_cover_and_ordered_gaps(self):
        for arity in range(2, 8):
            row = tuple(range(arity))
            prepared = geometry.prepare_geometry(arity, rowsets=(row,))
            self.assertFalse(prepared['metadata']['all_labelled_rowsets'])
            self.assertEqual(prepared['rowsets'], (row,))
            cover = information_cover(row)
            self.assertEqual(prepared['covers'][0], tuple(cover[i] for i in row))
            fiber = cross_fiber(row)
            choices = information_sets(fiber)
            actual_cover = tuple(F(sum(i in choice for choice in choices), len(choices)) for i in row)
            self.assertEqual(prepared['covers'][0], actual_cover)
            self.assertEqual(sum(prepared['covers'][0]), prepared['fiber_dimensions'][0])
            direct = gap_geometry(prepared['rowsets'], prepared['covers'], 16)
            self.assertEqual(prepared['gap_values'], direct[0])
            self.assertEqual(prepared['gap_parents'], direct[1])
            self.assertEqual(prepared['gap_labels'], direct[2])
            np.testing.assert_array_equal(prepared['gap_codes'], direct[3])
            cls = int(prepared['classes'][0])
            self.assertEqual(prepared['outside_caps'][cls][0], 1)
            self.assertEqual(prepared['nonquiet_caps'][cls][0], 0)

    def test_complete_pair_scope_and_quiet_transport(self):
        prepared = geometry.prepare_geometry(2)
        self.assertEqual(len(prepared['rowsets']), 120)
        self.assertTrue(prepared['metadata']['all_labelled_rowsets'])
        self.assertEqual(sum(prepared['orbit_sizes']), 120)
        self.assertTrue(all(cover == (F(0), F(0)) for cover in prepared['covers']))
        self.assertTrue(all(row == (1, 0, 0) for row in prepared['quiet_caps']))
        for index in (0, 17, 119):
            row = prepared['rowsets'][index]
            cls = int(prepared['classes'][index])
            self.assertEqual(prepared['outside_caps'][cls], geometry.actual_image_data(row)['full_caps'])
        self.assertFalse(prepared['metadata']['whole_code_certificate'])

    def test_per_support_max_is_not_weight_sum(self):
        counts = {0: 1, 1: 3, 2: 3, 4: 7, 3: 9}
        self.assertEqual(geometry.caps_by_weight(counts, 3), (1, 7, 9, 0))

    def test_guards_and_stop_before_preparation(self):
        for source in ((), (0,), tuple(range(8)), (1, 0), (0, 0), (0, 16)):
            with self.assertRaises(ValueError):
                geometry.actual_image_data(source)
        for arity in (1, 8, True):
            with self.assertRaises(ValueError):
                geometry.prepare_geometry(arity)
        with self.assertRaises(ValueError):
            geometry.prepare_geometry(2, rowsets=((0, 1), (0, 1)))
        with self.assertRaises(ValueError):
            geometry.prepare_all((2, 2))
        with self.assertRaises(ValueError):
            geometry.aggregate_outside_supports(1, {15: 1}, 2, 1, size=4, q=4)
        def stop():
            raise TimeoutError('stop before any complete preparation')
        with patch.object(geometry, 'quiet_cover_orbits') as called:
            with self.assertRaises(TimeoutError):
                geometry.prepare_geometry(7, check_time=stop)
            called.assert_not_called()


if __name__ == '__main__':
    unittest.main()

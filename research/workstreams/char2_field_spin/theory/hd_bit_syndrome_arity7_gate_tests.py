"""Exact local counts and bounded small contractions; no production gate."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import log2
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_bit_syndrome_arity7_gate as gate
from hd_bit_syndrome_arity9_gate import contract as old_contract
from hd_bit_syndrome_arity9_gate_tests import exact_column
from hd_bit_syndrome_cross_fibers import _rref, cross_fiber
from hd_bit_syndrome_kernel import affine_image_counts, information_cover


def mm(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


class AritySevenGateTests(TestCase):
    def setUp(self):
        self.constants = tuple(map(F, ('2/5', '1/3', '2', '3/5', '3/7', '1/11', '2/13')))
        self.measure = dict(zip(('log2_probability', 'log2_atom00', 'log2_atom01',
            'log2_atom10', 'log2_scale01', 'log2_scale10', 'log2_scale11'), map(log2, self.constants)))

    def test_binary_image_counts_equal_labelled_gf4_enumeration(self):
        generator = ((1, 0, 1, 1), (0, 1, 1, 0))
        counts = gate.image_support_counts(generator, q=4)
        direct = Counter()
        for x, y in product(range(4), repeat=2):
            word = tuple((x if a else 0) ^ (y if b else 0) for a, b in zip(*generator))
            direct[sum(bool(v) << i for i, v in enumerate(word))] += 1
        self.assertEqual(counts, tuple(direct[i] for i in range(16)))
        self.assertEqual(sum(counts), 16)
        # A non-binary one-dimensional AES image has the same support for
        # all255 nonzero amplitudes; rank arithmetic must retain its zero.
        self.assertEqual(gate.image_support_counts(((1, 2, 7, 0),)),
                         tuple(1 if i == 0 else 255 if i == 7 else 0 for i in range(16)))

    def test_actual_exception_is_not_full_affine_outside_code(self):
        source = tuple(range(7))
        caps, rank, counts = gate.actual_image_caps(source)
        fiber = cross_fiber(source)
        ambient = affine_image_counts(fiber['outside'])
        self.assertEqual(rank, 4)
        self.assertEqual(sum(counts), 256**4)
        self.assertEqual(sum(ambient.values()), 256**5)
        self.assertEqual(fiber['inside_zero_fiber_dimension'], 3)
        self.assertTrue(fiber['zero_outside_equals_quiet'])
        self.assertEqual(sum(information_cover(source).values()), 3)
        self.assertEqual(caps[0], 1)
        self.assertTrue(all(value <= caps[mask.bit_count()] for mask, value in enumerate(counts)))
        # Independent rank restrictions validate the complete support
        # inversion, not merely its total cardinality.
        for mask in (0, 1, 3, 0x55, 0xaa, 0x100, 0x1fe, 0x1ff):
            selected_rows = tuple(row for i, row in enumerate(fiber['cross_matrix']) if not mask >> i & 1)
            restriction_rank = len(_rref(selected_rows, 7)[1])
            direct = 256**(rank-restriction_rank)
            self.assertEqual(sum(value for sub, value in enumerate(counts) if sub & ~mask == 0), direct)

    def test_generic_actual_image_equals_affine_counts(self):
        source = (0, 1, 2, 3, 4, 5, 8)
        caps, rank, counts = gate.actual_image_caps(source)
        outside = tuple(i for i in range(16) if i not in source)
        exact = affine_image_counts(outside)
        self.assertEqual(rank, 5)
        for local, count in enumerate(counts):
            global_mask = sum(1 << i for j, i in enumerate(outside) if local >> j & 1)
            self.assertEqual(count, exact.get(global_mask, 0))
        self.assertEqual(sum(counts), 256**5)
        self.assertEqual(caps[0], 1)

    def test_small_ordered_contraction_matches_old_and_rational_reference(self):
        rows = tuple(combinations(range(4), 2))
        prepared = dict(rowsets=rows, covers=((F(1), F(0)),)*len(rows),
            outside_caps=((1, 3, 9), (1, 0, 7)), classes=np.arange(len(rows), dtype=np.int8) % 2)
        answer = gate.contract(prepared, self.measure, 2, size=4, columns=3, q=4, batch_size=2)
        control = old_contract(prepared, self.measure, 2, size=4, columns=3, q=4, batch_size=3)
        np.testing.assert_array_equal(answer['column_log2_matrices'], control['column_log2_matrices'])
        self.assertEqual(answer['combined_word_log2'], control['combined_word_log2'])
        for index in (0, 2, 5):
            direct = exact_column(rows[index], prepared['covers'][index],
                prepared['outside_caps'][index % 2], 2, self.constants)
            np.testing.assert_allclose(2**answer['column_log2_matrices'][index], direct, rtol=3e-14)
        self.assertTrue(answer['outside_zero_included'])
        self.assertFalse(answer['young_split'])

    def test_nine_marked_constant_branch_retains_order(self):
        prepared = dict(rowsets=((4,),), covers=((F(1),),),
            outside_caps=((1, 0, 0, 0, 0, 0, 0, 0, 0, 0),), classes=np.array([0]))
        answer = gate.contract(prepared, self.measure, 2, size=10, columns=1, q=4)
        p, b00, b01, b10, c01, c10, c11 = self.constants
        empty = [[F(1), c01], [c10, c11]]
        weighted = [[(1-p)+p*(1+3*b00)**2, c01*(1+3*b01)**2],
                    [c10*(1+3*b10)**2, c11*16]]
        direct = [[F(1), F(0)], [F(0), F(1)]]
        for position in range(10):
            direct = mm(direct, weighted if position == 4 else empty)
        np.testing.assert_allclose(2**answer['column_log2_matrices'][0],
                                  np.asarray(direct, dtype=float), rtol=5e-14)
        self.assertEqual(answer['costs']['final_atom_count_classes'], 715)

    def test_production_scope_without_full_preparation(self):
        prepared = dict(rowsets=tuple(combinations(range(16), 7)),
            exceptional_indices=tuple(range(240)), outside_caps=((1,)*10,), orbit_sizes=(11440,))
        with patch.object(gate, 'contract', return_value={}) as called:
            answer = gate.full_row_bounds(self.measure, prepared=prepared)
        self.assertEqual(called.call_count, 1)
        meta = answer['metadata']
        self.assertEqual(meta['outside_marked_positions'], 9)
        self.assertEqual(meta['rs_normalization_per_group'], '255^-56')
        self.assertEqual(meta['aggregate_input_floor'], 1134)
        self.assertEqual(meta['exceptional_actual_image_count'], 240)
        self.assertEqual(meta['exceptional_fiber_dimension'], 3)
        self.assertFalse(meta['floor_used_numerically'])
        self.assertFalse(meta['global_prefix_included'])
        self.assertFalse(meta['whole_code_certificate'])

    def test_guards_and_early_stop(self):
        for generator, kwargs in ((((1, 0), (1, 0)), {}),
                                  (((1, 2),), {'q': 4}),
                                  (((1,)*10,), {})):
            with self.assertRaises(ValueError):
                gate.image_support_counts(generator, **kwargs)
        with self.assertRaises(ValueError):
            gate.actual_image_caps((0, 1))
        def stop():
            raise TimeoutError('no production work')
        with self.assertRaises(TimeoutError):
            gate.full_row_bounds(self.measure, check_time=stop)


if __name__ == '__main__':
    main()

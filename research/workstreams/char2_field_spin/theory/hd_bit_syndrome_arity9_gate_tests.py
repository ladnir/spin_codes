"""Small exact/equivalence controls only; no arity-nine production gate."""
from fractions import Fraction as F
from itertools import combinations, product
from math import log2
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_bit_syndrome_arity9_gate as gate
from hd_bit_syndrome_cross_fibers import cross_fiber
from hd_bit_syndrome_gate import weighted_cover_matrices
from hd_bit_syndrome_gate_tests import atom_words, elementary
from hd_bit_syndrome_kernel import affine_image_counts, information_cover
from hd_post_refresh_quiet_rows import log_matrix_power


def exact_column(source, cover, caps, groups, constants, *, size=4, q=4):
    outside = tuple(i for i in range(size) if i not in source)
    result = [[F(0), F(0)], [F(0), F(0)]]
    for initial, final, markers, mass in atom_words(size, constants):
        e = elementary(tuple(markers[i] for i in outside))
        value = sum(count*entry for count, entry in zip(caps, e))
        for i, c in zip(source, cover):
            if c == 1:
                value *= 1+(q-1)*markers[i]
            elif c == 0:
                value *= max(1, markers[i])
            else:
                raise ValueError('Rational oracle uses deterministic information sets')
        result[initial][final] += mass*value**groups
    return np.asarray(result, dtype=float)


class ArityNineGateTests(TestCase):
    def setUp(self):
        self.constants = tuple(map(F, ('2/5', '1/3', '2', '3/5', '3/7', '1/11', '2/13')))
        self.measure = dict(zip(('log2_probability', 'log2_atom00', 'log2_atom01',
            'log2_atom10', 'log2_scale01', 'log2_scale10', 'log2_scale11'), map(log2, self.constants)))

    def test_small_orbit_counts_and_transported_covers(self):
        prepared = gate.quiet_cover_orbits(8, 5, q=4)
        self.assertEqual(len(prepared['rowsets']), 56)
        self.assertEqual(sum(prepared['orbit_sizes']), 56)
        for source, cover, cls in zip(prepared['rowsets'], prepared['covers'], prepared['classes']):
            exact = information_cover(source, bits=3)
            self.assertEqual(cover, tuple(exact[i] for i in source))
            outside = tuple(i for i in range(8) if i not in source)
            caps = prepared['outside_caps'][int(cls)]
            for mask, count in affine_image_counts(outside, bits=3, q=4).items():
                self.assertLessEqual(count, caps[mask.bit_count()])
            self.assertEqual(caps[0], 1)

    def test_actual_exceptional_fiber_not_quiet(self):
        for source, expected_rank in ((tuple(range(7, 16)), 4),
                                     ((6, 7, 9, 10, 11, 12, 13, 14, 15), 5)):
            quiet = information_cover(source)
            proposed = tuple(quiet[i] for i in source)
            cover, rank, dimension, exception = gate.source_cover(source, proposed)
            direct, actual_rank, actual_dimension = gate.actual_fiber_cover(source)
            self.assertEqual((rank, dimension), (expected_rank, 9-expected_rank))
            self.assertEqual((actual_rank, actual_dimension), (rank, dimension))
            self.assertEqual(cover, direct)
            self.assertEqual(exception, expected_rank == 4)
            self.assertEqual(sum(cover), dimension)
            self.assertEqual(cross_fiber(source)['quiet_dimension'], 4)
            if exception:
                self.assertNotEqual(cover, proposed)

    def test_actual_a9_outside_image_is_full_affine_image(self):
        # A basis of the actual C image has the same dimension as V_O.
        # Each C column lies in im(V_O) by H-I=V B V^T; the rank check
        # therefore proves equality, including the rank4 exception.
        for source in (tuple(range(9)), (0, 1, 2, 3, 4, 5, 6, 8, 9)):
            fiber = cross_fiber(source)
            expected = affine_image_counts(fiber['outside'])
            self.assertEqual(sum(expected.values()), 256**fiber['cross_rank'])
            self.assertEqual(expected[0], 1)
            self.assertEqual(fiber['inside_zero_fiber_dimension']+fiber['cross_rank'], 9)

    def test_complete_branch_matches_exact_ordered_atom_words(self):
        rows = tuple(combinations(range(4), 2))
        covers = tuple((F(1), F(0)) if index % 2 else (F(0), F(1)) for index in range(len(rows)))
        prepared = dict(rowsets=rows, covers=covers,
            outside_caps=((1, 2, 3), (1, 0, 7)), classes=np.arange(len(rows), dtype=np.int8) % 2)
        answer = gate.contract(prepared, self.measure, 2, size=4, columns=3, q=4, batch_size=2)
        for index, (source, cover) in enumerate(zip(rows, covers)):
            expected = exact_column(source, cover, prepared['outside_caps'][index % 2],
                                    2, self.constants)
            np.testing.assert_allclose(2**answer['column_log2_matrices'][index], expected, rtol=3e-14)
        # Distinct labels/gaps have not been pooled as whole row-set terms.
        self.assertFalse(np.allclose(answer['column_log2_matrices'][0],
                                     answer['column_log2_matrices'][-1]))
        self.assertTrue(answer['outside_zero_included'])
        self.assertFalse(answer['young_split'])
        self.assertEqual(answer['costs']['dp_passes'], 1)

    def test_minkowski_root_sum_and_batch_invariance(self):
        rows = tuple(combinations(range(4), 2))
        prepared = dict(rowsets=rows, covers=((F(1), F(0)),)*len(rows),
            outside_caps=((1, 3, 9),), classes=np.zeros(len(rows), dtype=np.int8))
        a = gate.contract(prepared, self.measure, 2, size=4, columns=3, q=4, batch_size=1)
        b = gate.contract(prepared, self.measure, 2, size=4, columns=3, q=4, batch_size=6)
        np.testing.assert_array_equal(a['column_log2_matrices'], b['column_log2_matrices'])
        words = log_matrix_power(a['column_log2_matrices'], 3)
        values = np.logaddexp2.reduce(words[:, 0, :], axis=1)
        self.assertEqual(a['combined_word_log2'], float(2*np.logaddexp2.reduce(values/2)))
        self.assertFalse(a['extra_rowset_factor'])

    def test_fractional_cover_translates_and_above_one_atoms(self):
        # GF4 repetition code has two information sets. Enumerate every
        # affine translate, including a coordinate whose marker exceeds1.
        for left, right in product((F(0), F(1, 3), F(2)), repeat=2):
            bound = ((1+3*float(left)**2)*(1+3*float(right)**2))**.5
            for shift in range(4):
                exact = sum((left if x else 1)*(right if x ^ shift else 1) for x in range(4))
                self.assertLessEqual(float(exact), bound+1e-13)
        bank = weighted_cover_matrices(self.measure, 2, (F(0), F(1, 2), F(1)), q=4)
        self.assertTrue(np.all(np.isfinite(bank['cover_matrices'])))

    def test_production_scope_without_preparation_or_contraction(self):
        rows = tuple(combinations(range(16), 9))
        prepared = dict(rowsets=rows, exceptional_indices=tuple(range(240)),
                        outside_caps=((1, 0, 0, 0, 0, 0, 0, 1),), orbit_sizes=(11440,))
        with patch.object(gate, 'contract', return_value={}) as called:
            answer = gate.full_row_bounds(self.measure, prepared=prepared)
        self.assertEqual(called.call_count, 1)
        meta = answer['metadata']
        self.assertEqual(meta['arity'], 9)
        self.assertEqual(meta['rs_normalization_per_group'], '255^-72')
        self.assertEqual(meta['aggregate_input_floor'], 1458)
        self.assertEqual(meta['exceptional_actual_fiber_count'], 240)
        self.assertFalse(meta['floor_used_numerically'])
        self.assertFalse(meta['global_prefix_included'])
        self.assertFalse(meta['whole_code_certificate'])

    def test_guards_and_deadline(self):
        with self.assertRaises(ValueError):
            gate.gap_geometry(((0, 0),), ((F(1), F(0)),), 4)
        with self.assertRaises(ValueError):
            gate.quiet_cover_orbits(7, 3)
        bad = dict(rowsets=((0,),), covers=((F(1),),), outside_caps=((0, 3),), classes=np.array([0]))
        with self.assertRaises(ValueError):
            gate.contract(bad, self.measure, 2, size=2)
        def stop():
            raise TimeoutError('no production work')
        with self.assertRaises(TimeoutError):
            gate.full_row_bounds(self.measure, check_time=stop)


if __name__ == '__main__':
    main()

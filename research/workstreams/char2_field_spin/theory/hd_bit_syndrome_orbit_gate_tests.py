"""Exact coefficient and tiny chronological controls; no production gate."""
from fractions import Fraction as F
from math import log2
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_bit_syndrome_gate as frozen
import hd_bit_syndrome_orbit_gate as gate
from hd_bit_syndrome_gate_tests import atom_words, elementary, exact_nonquiet_column
from hd_bit_syndrome_kernel import affine_image_counts
from hd_bit_syndrome_outer import quiet_support_counts


def logs(row):
    return np.array([log2(value) if value else -np.inf for value in row])


class BitSyndromeOrbitGateTests(TestCase):
    def setUp(self):
        self.constants = tuple(map(F, ('2/5', '1/3', '2', '3/5', '3/7', '1/11', '2/13')))
        self.measure = dict(zip(('log2_probability', 'log2_atom00', 'log2_atom01',
            'log2_atom10', 'log2_scale01', 'log2_scale10', 'log2_scale11'), map(log2, self.constants)))

    def test_exact_four_class_partition_and_improved_caps(self):
        quiet, outside = gate.coefficient_tables()
        classes = gate.labelled_classes()
        self.assertEqual(tuple(np.bincount(classes)), gate.ORBIT_SIZES)
        self.assertEqual(tuple(row[8] for row in quiet), gate.FULL_QUIET_COUNTS)
        self.assertEqual(len(classes), 12870)
        self.assertFalse(classes.flags.writeable)
        old_quiet, old_outside = frozen.coefficient_tables()
        old_outside = dict(old_outside)
        for index in range(4):
            self.assertTrue(all(a <= b for a, b in zip(quiet[index], old_quiet)))
            self.assertTrue(all(a <= b for a, b in zip(outside[index], old_outside[4 if index == 0 else 5])))
        self.assertEqual(tuple(row[1] for row in outside), (0, 255, 0, 0))
        self.assertTrue(any(quiet[2][i] < old_quiet[i] for i in range(9)))

    def test_caps_dominate_exact_labelled_subsets(self):
        quiet, outside = gate.coefficient_tables()
        table = quiet_support_counts()
        rows = frozen.row_geometry()[0]
        classes = gate.labelled_classes()
        indices = set((0, 1, 11, 511, 6000, len(rows)-1))
        for cls in range(4):
            selected = np.flatnonzero(classes == cls)
            indices.update((int(selected[0]), int(selected[-1])))
        for index in indices:
            source, cls = rows[index], int(classes[index])
            mask = sum(1 << i for i in source)
            subset = mask
            while True:
                self.assertLessEqual(table[subset], quiet[cls][subset.bit_count()])
                if not subset:
                    break
                subset = (subset-1) & mask
            complement = tuple(i for i in range(16) if i not in source)
            for support, count in affine_image_counts(complement).items():
                if support:
                    self.assertLessEqual(count, outside[cls][support.bit_count()])

    def test_equal_tables_reproduce_frozen_two_branches(self):
        quiet, outside = logs((1, 2, 0, 3, 4)), logs((0, 2, 1, 3, 4))
        rows = frozen.row_geometry(8, 4)[0]
        classes = np.arange(len(rows), dtype=np.int8) % 2
        for cover in (False, True):
            reference = frozen.contract(quiet, {3: outside, 4: outside}, self.measure,
                2, 1., size=8, columns=3, q=4, batch_size=7, information_cover=cover)
            result = gate.contract(np.stack((quiet, quiet)), np.stack((outside, outside)),
                classes, self.measure, 2, 1., size=8, columns=3, q=4,
                batch_size=7, information_cover=cover)
            for branch in ('quiet', 'nonquiet'):
                np.testing.assert_array_equal(result['branch_column_log2_matrices'][branch],
                                              reference['branch_column_log2_matrices'][branch])
            self.assertEqual(result['combined_word_log2'], reference['combined_word_log2'])
            self.assertEqual(result['costs']['total_dp_passes'], 2)
            self.assertEqual(result['costs']['quiet']['dp_passes'], 1)
            self.assertEqual(result['costs']['nonquiet']['dp_passes'], 1)

    def test_class_terminal_selection_matches_independent_contractions(self):
        quiet = np.stack((logs((1, 2, 0, 3, 4)), logs((1, 0, 2, 1, 7))))
        outside = np.stack((logs((0, 2, 1, 3, 4)), logs((0, 1, 3, 2, 5))))
        rows = frozen.row_geometry(8, 4)[0]
        classes = np.array([int(0 in row) for row in rows], dtype=np.int8)
        result = gate.contract(quiet, outside, classes, self.measure,
            2, 1., size=8, columns=3, q=4, batch_size=5)
        for cls in range(2):
            reference = frozen.contract(quiet[cls], {3: outside[cls], 4: outside[cls]},
                self.measure, 2, 1., size=8, columns=3, q=4, batch_size=13)
            chosen = classes == cls
            for branch in ('quiet', 'nonquiet'):
                np.testing.assert_array_equal(result['branch_column_log2_matrices'][branch][chosen],
                                              reference['branch_column_log2_matrices'][branch][chosen])
        # Direct rational ordered atom words check the new class-selected
        # quiet DP independently of the frozen matrix implementation.
        index = 1
        source, cls = rows[index], int(classes[index])
        counts = ((1, 2, 0, 3, 4), (1, 0, 2, 1, 7))[cls]
        expected = [[F(0), F(0)], [F(0), F(0)]]
        for initial, end, markers, mass in atom_words(8, self.constants):
            e = elementary(tuple(markers[i] for i in source))
            scalar = sum(count * value * 2**weight for weight, (count, value) in enumerate(zip(counts, e)))
            expected[initial][end] += mass * scalar**2
        np.testing.assert_allclose(2**result['branch_column_log2_matrices']['quiet'][index],
                                  np.asarray(expected, dtype=float), rtol=8e-14)
        info = result['information_sets'][index]
        counts_out = ((0, 2, 1, 3, 4), (0, 1, 3, 2, 5))[cls]
        expected = exact_nonquiet_column(source, info, counts_out, 2, self.constants, size=8, q=4)
        np.testing.assert_allclose(2**result['branch_column_log2_matrices']['nonquiet'][index],
                                  np.asarray(expected, dtype=float), rtol=8e-14)

    def test_minkowski_sums_every_labelled_row_once(self):
        quiet = np.stack((logs((1, 1)), logs((1, 2))))
        outside = np.stack((logs((0, 2)), logs((0, 3))))
        result = gate.contract(quiet, outside, np.array([0, 1]), self.measure, 2, 1.,
            size=2, columns=2, q=4, lambdas=(F(3, 4), F(1, 4)))
        self.assertEqual(result['rowsets'], ((0,), (1,)))
        word = gate.log_matrix_power(result['column_log2_matrices'], 2)
        values = np.logaddexp2.reduce(word[:, 0, :], axis=1)
        np.testing.assert_array_equal(values, result['per_rowset_word_log2'])
        self.assertEqual(result['combined_word_log2'], float(2*np.logaddexp2.reduce(values/2)))
        self.assertFalse(result['extra_rowset_factor'])
        self.assertFalse(result['orbit_multiplicity_factor'])
        self.assertTrue(result['branch_sum_before_word_power'])

    def test_production_scope_without_contraction(self):
        with patch.object(gate, 'contract', return_value={}) as called:
            result = gate.full_row_bounds(self.measure, 18, 4., information_cover=True)
        self.assertEqual(called.call_count, 1)
        self.assertEqual(called.call_args.args[0].shape, (4, 9))
        self.assertEqual(called.call_args.args[1].shape, (4, 9))
        self.assertEqual(len(called.call_args.args[2]), 12870)
        self.assertTrue(called.call_args.kwargs['information_cover'])
        self.assertEqual(result['metadata']['intended_inner'], 'original-one-state')
        self.assertEqual(result['metadata']['total_dp_passes'], 2)
        self.assertFalse(result['metadata']['global_prefix_included'])
        self.assertFalse(result['metadata']['exact_sparse_restoration'])
        self.assertFalse(result['metadata']['directed'])
        self.assertFalse(result['metadata']['whole_code_certificate'])

    def test_validation_and_budget(self):
        with self.assertRaises(ValueError):
            gate.contract([[0., 0.]], [[0., 0.]], [0, 0], self.measure, 2, 0., size=2)
        with self.assertRaises(ValueError):
            gate.contract([[0., 0.]], [[-np.inf, 0.]], [0, 1], self.measure, 2, 0., size=2)
        def stop():
            raise TimeoutError('toy only')
        with self.assertRaises(TimeoutError):
            gate.full_row_bounds(self.measure, 18, 4., check_time=stop)


if __name__ == '__main__':
    main()

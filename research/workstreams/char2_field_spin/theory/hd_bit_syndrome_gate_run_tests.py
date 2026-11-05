"""Small loader/normalization tests; no production contraction."""
from copy import deepcopy
from math import log2
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_bit_syndrome_gate_run as gate


class GateRunTests(TestCase):
    def test_saved_distribution_and_moment_cover(self):
        data, curves = gate.checked_bank()
        measure = gate.entry_measure(curves)
        m = np.arange(19)
        bound = np.array([
            measure['log2_probability'] + measure['log2_atom00'] * m,
            measure['log2_scale01'] + measure['log2_atom01'] * m,
            measure['log2_scale10'] + measure['log2_atom10'] * m,
            np.full(19, measure['log2_scale11'])])
        bound[0, 0] = 0
        self.assertTrue(np.all(bound >= curves - 1e-9))
        self.assertEqual(data['source_count'], str(gate.SOURCE_COUNT))
        self.assertFalse(data['whole_code_certificate'])

    def test_reject_changed_scope(self):
        data, curves = gate.checked_bank()
        for key, value in (('source_full_pair_max', 2), ('rowsets_resampled', True),
                           ('terminal_vector', [1, 0]), ('source_count', '1'),
                           ('aggregate_input_weight_floor', 1295)):
            changed = deepcopy(data)
            changed[key] = value
            with self.subTest(key=key), patch.object(gate, 'load_saved', return_value=(changed, curves)):
                with self.assertRaises(ValueError):
                    gate.checked_bank()

    def test_reject_bad_curves(self):
        data, curves = gate.checked_bank()
        for bad in (curves[:, :-1], np.full((4, 19), np.nan), np.full((4, 19), np.inf)):
            with patch.object(gate, 'load_saved', return_value=(data, bad)):
                with self.assertRaises(ValueError):
                    gate.checked_bank()

    def test_prefix_has_exactly_one_projective_and_source_factor(self):
        expected = (log2(gate.SOURCE_COUNT) - gate.LIMIT * log2(691 / 700)
                    - 1152 * log2(255) - log2(255))
        self.assertAlmostEqual(gate.prefix(0), expected, places=10)
        self.assertAlmostEqual(gate.prefix(4), expected - 1296 * 4, places=10)
        for invalid in (-1, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                gate.prefix(invalid)

    def test_whole_word_and_fixed_rowset_minkowski(self):
        matrices = np.full((3, 2, 2), -np.inf)
        matrices[:, 0, 0] = matrices[:, 1, 1] = 1
        bound, values = gate.word_bound(matrices)
        np.testing.assert_array_equal(values, [16, 16, 16])
        self.assertAlmostEqual(bound, 16 + 18 * log2(3), places=10)


if __name__ == '__main__':
    main()

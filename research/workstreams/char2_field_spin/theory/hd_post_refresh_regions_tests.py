from itertools import product
from math import comb
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_post_refresh_regions as regions
import hd_coupled_state_regions as generic


def make_packet(terms, maximum, full_max=None, marker=1.):
    full_max = regions.checked_scope(maximum, full_max)
    result = np.zeros((8, 8, maximum+1, full_max+1))
    for s, f, matrix in terms:
        if f <= full_max:
            result[:, :, s, f] += marker**(s+f)*matrix
    return result


class PostRefreshRegionsTests(TestCase):
    def setUp(self):
        self.empty, self.single, self.full = (np.eye(8, dtype=int) for _ in range(3))
        self.empty[0, 4], self.empty[4, 5] = 1, 2
        self.single[5, 6], self.single[6, 7] = 3, 1
        self.full[7, 0], self.full[1, 2] = 2, 1
        self.terms = [(0, 0, self.empty), (1, 0, self.single), (0, 1, self.full)]

    def test_exact_noncommuting_packet_words(self):
        maximum, packets, marker = 4, 4, .125
        self.assertFalse(np.array_equal(self.empty@self.single, self.single@self.empty))
        expected = {}
        for word in product(self.terms, repeat=packets):
            single, full = sum(term[0] for term in word), sum(term[1] for term in word)
            if single+2*full > maximum:
                continue
            coefficient = np.eye(8, dtype=object)
            for _, _, matrix in word:
                coefficient = coefficient@matrix.astype(object)
            key = single, full
            expected[key] = expected.get(key, np.zeros((8, 8), dtype=object))+coefficient
        actual = regions.regions_from_packet(make_packet(self.terms, maximum, marker=marker),
            maximum, packets=packets, pairs_per_packet=1, marker=marker)['regions']
        for (single, full), coefficient in expected.items():
            denominator = comb(packets, single)*comb(packets-single, full)
            np.testing.assert_allclose(np.exp2(actual[single+2*full, full]),
                np.asarray(coefficient, dtype=float)/denominator, rtol=3e-14, atol=1e-14)

    def test_reuses_frozen_generic_power_and_exact_full_cutoff(self):
        self.assertIs(regions.power_coefficients, generic.power_coefficients)
        self.assertIs(regions.product_coefficients, generic.product_coefficients)
        full = regions.power_coefficients(make_packet(self.terms, 4), 4, 4)
        narrow = regions.power_coefficients(make_packet(self.terms, 4, 1), 4, 4, full_max=1)
        np.testing.assert_array_equal(narrow, full[:, :, :, :2])

    def test_orientation_and_identity_normalization_once(self):
        geometry = [(0, 0, 1, 'e'), (1, 0, 1, 's'), (1, 0, 1, 'f'), (0, 1, 1, 's')]
        table = dict(e=self.empty, s=self.single, f=self.full)
        base = regions.packet_from_table(table, 2, marker=.5, geometry=geometry)
        np.testing.assert_array_equal(base[:, :, 1, 0], (self.single+self.full)/4)
        np.testing.assert_array_equal(base[:, :, 0, 1], self.single/2)
        slots, maximum = 2, 4
        terms = [(s, f, comb(slots, s)*comb(slots-s, f)*np.eye(8))
                 for f in range(3) for s in range(maximum-2*f+1) if s+f <= slots]
        for marker in (1., .125):
            result = regions.regions_from_packet(make_packet(terms, maximum, marker=marker),
                maximum, packets=4, pairs_per_packet=slots, marker=marker)['regions']
            for matrix in result.values():
                np.testing.assert_allclose(np.exp2(matrix), np.eye(8), atol=2e-14)

    def test_mocked_public_bank_has_complete36_histograms_and_new_identity(self):
        # Only the all-zero packet coefficient is present. This exercises128
        # packet powers and normalization scope without any physical cap table.
        with patch.object(regions, '_packet_table', return_value={'empty': np.eye(8)}) as build, \
                patch.object(regions, 'sparse_geometry', return_value=[(0, 0, 1, 'empty')]):
            answer = regions.raw_region_bank(.9, 18, full_max=1)
        build.assert_called_once_with(.9)
        self.assertEqual(set(answer['regions']), {(m, f) for m in range(19)
            for f in range(min(1, m//2)+1)})
        self.assertEqual(len(answer['regions']), 36)
        self.assertEqual(answer['metadata']['construction'], 'proposed-g16-post-refresh-coupled-inner')
        self.assertEqual(answer['metadata']['state_order'], ['Z', 'L', 'R', 'B', 'D', 'P', 'F', 'G'])
        self.assertEqual(answer['metadata']['global_normalization_count'], 1)
        self.assertEqual(answer['metadata']['total_pairs'], 4096)
        self.assertTrue(answer['metadata']['cone_F_conditional_second_flat_density'])
        self.assertFalse(answer['metadata']['tensor_factorization'])
        self.assertTrue(np.isneginf(answer['regions'][18, 1]).all())

    def test_dimension_and_coefficient_guards(self):
        with self.assertRaises(ValueError):
            regions.regions_from_packet(np.zeros((5, 5, 3, 2)), 2)
        with self.assertRaises(ValueError):
            regions.raw_region_bank(.9, 18, full_max=10)
        with self.assertRaises(ValueError):
            regions.raw_region_bank(1., 18)
        bad = make_packet(self.terms, 2)
        bad[0, 0, 2, 1] = 1.
        with self.assertRaises(ValueError):
            regions.regions_from_packet(bad, 2)

    def test_saved_bank_serialization_roundtrip_preserves_zero(self):
        with np.errstate(divide='ignore'):
            matrix = np.log2(np.eye(8))
        bank = {(m, f): matrix.copy() for m in range(19) for f in range(min(1, m//2)+1)}
        records = regions.serialize_bank(bank)
        self.assertEqual(len(records), 36)
        self.assertEqual(records[0]['log2_matrix'][0], [0., None, None, None, None, None, None, None])
        for row in records:
            restored = np.array([[-np.inf if value is None else value for value in values]
                                 for values in row['log2_matrix']])
            np.testing.assert_array_equal(restored, bank[row['active_weight'], row['full_pairs']])
        bank.pop((18, 1))
        with self.assertRaises(ValueError):
            regions.serialize_bank(bank)

    def test_frozen_saved_builder_with_mocked_bank_no_provider(self):
        with np.errstate(divide='ignore'):
            matrix = np.log2(np.eye(8))
        bank = {(m, f): matrix for m in range(19) for f in range(min(1, m//2)+1)}
        metadata = dict(construction=regions.CONSTRUCTION, state_order=list(regions.STATE_ORDER),
            global_normalization_count=1, regional_alpha=1, directed=False, whole_code_certificate=False)
        with patch.object(regions, 'raw_region_bank', return_value=dict(regions=bank, metadata=metadata)) as build:
            result = regions.build_saved_bank()
        build.assert_called_once_with(float(regions.FROZEN_MARK), 18, full_max=1)
        self.assertTrue(result['complete_requested_bank'])
        self.assertEqual(result['mark'], '691/700')
        self.assertEqual(result['regional_histogram_count'], 36)
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['reviewed_source_sha256'], regions.SOURCE_PINS)
        self.assertEqual(len(result['source_sha256']), 64)
        with self.assertRaises(ValueError):
            regions.build_saved_bank(max_seconds=31)


if __name__ == '__main__':
    main()

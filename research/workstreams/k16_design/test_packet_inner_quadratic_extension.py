"""Light rank and mocked preparation tests; no 19..22-state census runs."""
from copy import deepcopy
import hashlib
from itertools import combinations
from pathlib import Path
import unittest
from unittest.mock import patch

import packet_inner_quadratic_extension as inner


def fixture(bits):
    rows, columns, record = inner.construction(bits)
    spectrum = {'0': 1, '16': (1 << bits)-1}
    dual = ['0']*65
    dual[0], dual[6] = '1', str((1 << (64-bits))-1)
    record.update(full_state_census=True, state_count=1 << bits, expansion_spectrum=spectrum,
        feedback_transpose_spectrum=dict(spectrum), feedback_kernel_spectrum=dual,
        minimum_expansion_weight=16, minimum_feedback_kernel_weight=6,
        birth_density='capped', exact_profiles_regenerated=True)
    physical = dict(bits=bits, columns=list(columns), map_images={1 << i: row for i, row in enumerate(rows)})
    data = dict(bits=bits, physical_step_bits=64, macro_step_bits=128, macro_windows=32,
                birth_density='capped', map_sha256=record['map_sha256'])
    return rows, columns, data, physical, record


class QuadraticExtensionTests(unittest.TestCase):
    def test_each_prefix_rank_transpose_CA_and_packets_without_enumeration(self):
        base = inner.extension.construction()[0][:16]
        previous = base
        for bits in inner.SUPPORTED_BITS:
            with patch.object(inner.extension, 'enumerate_census', side_effect=AssertionError('no census')):
                rows, columns, record = inner.construction(bits)
            self.assertEqual(rows[:len(previous)], previous)
            self.assertEqual(rows[:16], base)
            self.assertEqual(inner.s16_maps.binary_rank(rows), bits)
            self.assertEqual(inner.s16_maps.binary_rank(columns), bits)
            self.assertEqual(columns, inner.extension.columns_from_rows(rows, 64))
            self.assertEqual(len(set(columns)), 64)
            self.assertNotIn(0, columns)
            self.assertTrue(all((a & b).bit_count() % 2 == 0 for a in rows for b in rows))
            self.assertEqual(record['packet_ranks'], [4]*16)
            self.assertEqual(record['return_denominator'], (1 << bits)-1)
            self.assertEqual(record['full_RM_2_6'], bits == 22)
            self.assertFalse(record['full_state_census'])
            previous = rows

    def test_rank22_is_the_full_quadratic_space_and_s20_map_is_unchanged(self):
        monomials = [(1 << 64)-1]
        monomials.extend(sum(1 << x for x in range(64) if x >> i & 1) for i in range(6))
        monomials.extend(inner.extension.quadratic_monomial(i, j) for i, j in combinations(range(6), 2))
        rows = inner.construction(22)[0]
        self.assertEqual(len(monomials), 22)
        self.assertEqual(inner.s16_maps.binary_rank(monomials), 22)
        self.assertEqual(inner.s16_maps.binary_rank([*rows, *monomials]), 22)
        old_rows, old_columns, old_record = inner.extension.construction()
        new_rows, new_columns, new_record = inner.construction(20)
        self.assertEqual((new_rows, new_columns), (old_rows, old_columns))
        self.assertEqual(new_record['map_sha256'], old_record['map_sha256'])
        self.assertNotEqual(new_record['schema'], old_record['schema'])
        self.assertEqual(inner.extension.STATE_BITS, 20)
        self.assertEqual(inner.metadata.SUPPORTED_BITS, (17, 18))

    def test_all_construction_and_metadata_helper_pins_are_current(self):
        record = inner.construction(22)[2]
        expected = {str(inner.BASE_MAP), str(Path(inner.__file__).resolve()),
            str(Path(inner.extension.__file__).resolve()), str(Path(inner.s16_maps.__file__).resolve()),
            str(Path(inner.metadata.__file__).resolve())}
        self.assertEqual(set(record['source_sha256']), expected)
        for path, digest in record['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(path).read_bytes()).hexdigest(), digest)

    def test_rejects_invalid_dimensions_and_density_before_census(self):
        for bits in (None, True, 18, 23, 19.0, '19'):
            with self.subTest(bits=bits), patch.object(inner.extension, 'enumerate_census') as census:
                with self.assertRaises(ValueError):
                    inner.prepare(bits)
                census.assert_not_called()
        with patch.object(inner.extension, 'enumerate_census') as census, self.assertRaises(ValueError):
            inner.prepare(19, 'invalid')
        census.assert_not_called()

    def test_prepare_forwards_actual_bits_and_regenerates_profiles(self):
        for bits in inner.SUPPORTED_BITS:
            rows, columns, data, physical, record = fixture(bits)
            images = ('mocked-no-full-census',)
            spectrum = {int(w): n for w, n in record['expansion_spectrum'].items()}
            dual = tuple(map(int, record['feedback_kernel_spectrum']))
            with (patch.object(inner.extension, 'enumerate_census', return_value=(images, spectrum, dual)) as census,
                  patch.object(inner.kernel_maps, 'prepare_maps', return_value=physical) as prepare,
                  patch.object(inner.kernel_t64, 'wrap', return_value=data),
                  patch.object(inner.kernel_t64, 'authenticate', return_value=physical)):
                _, result = inner.prepare(bits)
            census.assert_called_once_with(rows, 64)
            prepare.assert_called_once_with(images, columns, bits=bits, distribution='uniform_gl', birth_density='capped')
            self.assertEqual(result['state_count'], 1 << bits)
            self.assertTrue(result['full_state_census'])
            self.assertFalse(result['whole_code_certificate'])

    def test_authentication_rejects_substituted_maps_metadata_or_missing_pins(self):
        _, _, data, physical, record = fixture(19)
        with patch.object(inner.kernel_t64, 'authenticate', return_value=physical):
            self.assertIs(inner.authenticate(data, record), physical)
        for key, value in (('s', 20), ('return_denominator', 65535), ('full_state_census', False),
                           ('state_count', 1 << 18), ('sampling', 'uniform GL16'), ('source_sha256', {})):
            changed = dict(record, **{key: value})
            with self.subTest(key=key), patch.object(inner.kernel_t64, 'authenticate', return_value=physical), self.assertRaises(ValueError):
                inner.authenticate(data, changed)
        wrong_physical = deepcopy(physical)
        wrong_physical['map_images'][1] ^= 1
        with patch.object(inner.kernel_t64, 'authenticate', return_value=wrong_physical), self.assertRaises(ValueError):
            inner.authenticate(data, record)


if __name__ == '__main__':
    unittest.main()

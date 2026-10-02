"""Structural and mocked-control-flow tests; no17/18-bit census is run."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import packet_inner_small_extension as inner


def prepared_fixture(bits):
    """Synthetic spectra exercise metadata validation, not a numerical claim."""
    rows, columns, record = inner.construction(bits)
    spectrum = {'0': 1, '16': (1 << bits)-1}
    dual = ['0']*65
    dual[0], dual[6] = '1', str((1 << (64-bits))-1)
    record.update(full_state_census=True, state_count=1 << bits,
        expansion_spectrum=spectrum, feedback_transpose_spectrum=dict(spectrum),
        feedback_kernel_spectrum=dual, minimum_expansion_weight=16,
        minimum_feedback_kernel_weight=6, birth_density='capped', exact_profiles_regenerated=True)
    physical = dict(bits=bits, columns=list(columns),
                    map_images={1 << j: row for j, row in enumerate(rows)})
    data = dict(bits=bits, physical_step_bits=64, macro_step_bits=128, macro_windows=32,
                birth_density='capped', map_sha256=record['map_sha256'])
    return rows, columns, data, physical, record


class SmallExtensionTests(unittest.TestCase):
    def test_prefixes_preserve_base_and_are_not_s20_or_each_other(self):
        base_bytes = inner.BASE_MAP.read_bytes()
        base = json.loads(base_bytes)
        full_rows, full_columns, full_record = inner.extension.construction()
        identities = {full_record['map_sha256']}
        for bits in (17, 18):
            with patch.object(inner.extension, 'enumerate_census', side_effect=AssertionError('no census')):
                rows, columns, record = inner.construction(bits)
            self.assertEqual(rows, full_rows[:bits])
            self.assertEqual(rows[:16], tuple(int(row, 16) for row in base['generator_rows_hex']))
            self.assertEqual(columns, tuple(column & ((1 << bits)-1) for column in full_columns))
            self.assertEqual(tuple(column & 65535 for column in columns), tuple(base['columns']))
            self.assertEqual(record['appended_monomials'], [[0, j] for j in range(1, bits-15)])
            self.assertEqual(record['appended_rows_hex'], list(map(hex, rows[16:])))
            self.assertEqual(record['return_denominator'], (1 << bits)-1)
            self.assertFalse(record['full_state_census'])
            self.assertFalse(record['whole_code_certificate'])
            self.assertNotIn(record['map_sha256'], identities)
            identities.add(record['map_sha256'])
        self.assertEqual(inner.BASE_MAP.read_bytes(), base_bytes)
        self.assertEqual(inner.extension.STATE_BITS, 20)

    def test_ranks_transpose_packet_ranks_and_CA_zero(self):
        for bits in (17, 18):
            rows, columns, record = inner.construction(bits)
            self.assertEqual(inner.s16_maps.binary_rank(rows), bits)
            self.assertEqual(inner.s16_maps.binary_rank(columns), bits)
            self.assertEqual(columns, inner.extension.columns_from_rows(rows, 64))
            self.assertEqual(len(set(columns)), 64)
            self.assertNotIn(0, columns)
            self.assertTrue(all((a & b).bit_count() % 2 == 0 for a in rows for b in rows))
            self.assertEqual(record['packet_ranks'], [4]*16)
            self.assertTrue(all(inner.extension._quadratic(row) for row in rows))
            self.assertEqual(record['sampling'], f'independent uniform GL{bits} for every physical t64 step')
            self.assertTrue(record['zero_initial_state'])
            self.assertFalse(record['final_flush'])
            self.assertEqual(record['macro']['state_continuity'], 'retained_between_halves')

    def test_exact_map_identity_and_current_construction_pins(self):
        expected_paths = {str(inner.BASE_MAP), str(Path(inner.__file__).resolve()),
                          str(Path(inner.extension.__file__).resolve()),
                          str(Path(inner.s16_maps.__file__).resolve())}
        for bits in (17, 18):
            rows, columns, record = inner.construction(bits)
            self.assertEqual(set(record['source_sha256']), expected_paths)
            for path, digest in record['source_sha256'].items():
                self.assertEqual(hashlib.sha256(Path(path).read_bytes()).hexdigest(), digest)
            self.assertEqual(record['source']['sha256'], inner.BASE_SHA256)
            payload = dict(bits=bits, width=64, expansion_rows=list(map(hex, rows)),
                           feedback_columns=list(columns))
            expected = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            self.assertEqual(record['map_sha256'], expected)

    def test_invalid_options_fail_before_enumeration(self):
        for bits in (None, True, 16, 19, 20, 17.0, '17'):
            with self.subTest(bits=bits):
                with self.assertRaises(ValueError):
                    inner.construction(bits)
                with patch.object(inner.extension, 'enumerate_census') as census, self.assertRaises(ValueError):
                    inner.prepare(bits)
                census.assert_not_called()
        with patch.object(inner.extension, 'enumerate_census') as census, self.assertRaises(ValueError):
            inner.prepare(17, birth_density='unsupported')
        census.assert_not_called()

    def test_prepare_requests_exact_state_width_and_fresh_profiles(self):
        for bits in (17, 18):
            for density in ('classes', 'capped'):
                rows, columns, data, physical, record = prepared_fixture(bits)
                data['birth_density'] = density
                images = ('mocked-only-no-state-census',)
                spectrum = {int(k): v for k, v in record['expansion_spectrum'].items()}
                dual = tuple(map(int, record['feedback_kernel_spectrum']))
                with (patch.object(inner.extension, 'enumerate_census', return_value=(images, spectrum, dual)) as census,
                      patch.object(inner.kernel_maps, 'prepare_maps', return_value=physical) as prepare,
                      patch.object(inner.kernel_t64, 'wrap', return_value=data),
                      patch.object(inner.kernel_t64, 'authenticate', return_value=physical)):
                    actual_data, actual_record = inner.prepare(bits, density)
                census.assert_called_once_with(rows, 64)
                prepare.assert_called_once_with(images, columns, bits=bits,
                                               distribution='uniform_gl', birth_density=density)
                self.assertIs(actual_data, data)
                self.assertEqual(actual_record['s'], bits)
                self.assertEqual(actual_record['state_count'], 1 << bits)
                self.assertEqual(actual_record['birth_density'], density)
                self.assertEqual(actual_record['return_denominator'], (1 << bits)-1)
                self.assertTrue(actual_record['full_state_census'])
                self.assertTrue(actual_record['exact_profiles_regenerated'])
                self.assertFalse(actual_record['whole_code_certificate'])

    def test_authenticate_accepts_only_matching_current_preparation(self):
        for bits in (17, 18):
            _, _, data, physical, record = prepared_fixture(bits)
            with patch.object(inner.kernel_t64, 'authenticate', return_value=physical):
                self.assertIs(inner.authenticate(data, record), physical)
            for key, value in (('bits', 16), ('physical_step_bits', 128), ('macro_windows', 16),
                               ('macro_step_bits', 64), ('birth_density', 'classes'),
                               ('map_sha256', 'wrong')):
                changed = dict(data, **{key: value})
                with self.subTest(bits=bits, key=key), patch.object(
                        inner.kernel_t64, 'authenticate', return_value=physical), self.assertRaises(ValueError):
                    inner.authenticate(changed, record)
            for field in ('columns', 'map_images'):
                bad = deepcopy(physical)
                if field == 'columns':
                    bad[field][0] ^= 1
                else:
                    bad[field][1] ^= 1
                with self.subTest(bits=bits, field=field), patch.object(
                        inner.kernel_t64, 'authenticate', return_value=bad), self.assertRaises(ValueError):
                    inner.authenticate(data, record)

    def test_saved_record_substitution_flags_spectra_and_pins_fail_closed(self):
        _, _, data, physical, record = prepared_fixture(17)
        cases = [('s', 16), ('s', 18), ('schema', inner.extension.SCHEMA),
            ('return_denominator', 65535), ('sampling', 'independent uniform GL16 for every physical t64 step'),
            ('full_state_census', False), ('full_state_census', 1), ('state_count', 1 << 16),
            ('exact_profiles_regenerated', False), ('whole_code_certificate', True),
            ('final_flush', True), ('zero_initial_state', 1), ('birth_density', 'invalid'),
            ('minimum_expansion_weight', 15), ('minimum_feedback_kernel_weight', 5),
            ('expansion_spectrum', None), ('feedback_transpose_spectrum', {}),
            ('feedback_kernel_spectrum', ['1']*65)]
        for key, value in cases:
            changed = deepcopy(record)
            changed[key] = value
            with self.subTest(key=key, value=value), patch.object(
                    inner.kernel_t64, 'authenticate', return_value=physical), self.assertRaises(ValueError):
                inner.authenticate(data, changed)
        for path in record['source_sha256']:
            for change in ('missing', 'changed'):
                changed = deepcopy(record)
                if change == 'missing':
                    del changed['source_sha256'][path]
                else:
                    changed['source_sha256'][path] = '0'*64
                with self.subTest(path=path, change=change), patch.object(
                        inner.kernel_t64, 'authenticate', return_value=physical), self.assertRaises(ValueError):
                    inner.authenticate(data, changed)
        with self.assertRaises(ValueError):
            inner.authenticate(data, None)

    def test_small_real_generic_kernel_rejects_substitution(self):
        # Eight images only: exercise the actual generic engine without an
        # expensive17/18-bit preparation or any retained map mutation.
        rows = (3, 5, 9)
        columns = inner.extension.columns_from_rows(rows, 4)
        images, _, _ = inner.extension.enumerate_census(rows, 4)
        physical = inner.kernel_maps.prepare_maps(images, columns, bits=3,
            distribution='uniform_gl', birth_density='classes')
        wrapper = inner.kernel_t64.wrap(physical)
        _, _, _, _, record = prepared_fixture(17)
        with self.assertRaises(ValueError):
            inner.authenticate(wrapper, record)
        self.assertEqual(inner.extension.STATE_BITS, 20)


if __name__ == '__main__':
    unittest.main()

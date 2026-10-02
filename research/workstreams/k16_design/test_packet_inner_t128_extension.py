"""Light structural and dispatch tests; no s20/s22 state census is run."""
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
from pathlib import Path
import unittest
from unittest.mock import patch

import packet_inner_t128_extension as inner


def fixture(bits, density='capped'):
    rows, columns, record = inner.construction(bits)
    spectrum = {'0': 1, '32': (1 << bits)-1}
    dual = ['0']*129
    dual[0], dual[6] = '1', str((1 << (128-bits))-1)
    record.update(full_state_census=True, state_count=1 << bits, expansion_spectrum=spectrum,
        feedback_transpose_spectrum=dict(spectrum), feedback_kernel_spectrum=dual,
        minimum_expansion_weight=32, minimum_feedback_kernel_weight=6,
        birth_density=density, exact_profiles_regenerated=True)
    data = dict(bits=bits, windows=32, columns=list(columns), distribution='uniform_gl',
        map_images={1 << i: row for i, row in enumerate(rows)}, birth_density=density,
        map_sha256=record['map_sha256'])
    return rows, columns, data, record


class T128ExtensionTests(unittest.TestCase):
    def test_map_ranks_transpose_CA_packets_and_preserved_prefix(self):
        previous = ()
        for bits in inner.SUPPORTED_BITS:
            with patch.object(inner, 'enumerate_census', side_effect=AssertionError('no full census')):
                rows, columns, record = inner.construction(bits)
            self.assertEqual(rows[:len(previous)], previous)
            self.assertEqual(inner.s16_maps.binary_rank(rows), bits)
            self.assertEqual(inner.s16_maps.binary_rank(columns), bits)
            self.assertEqual(columns, inner.columns_from_rows(rows))
            self.assertTrue(all(inner._quadratic(row) for row in rows))
            self.assertTrue(all((a & b).bit_count() % 2 == 0 for a in rows for b in rows))
            self.assertEqual(len(set(columns)), 128)
            self.assertNotIn(0, columns)
            self.assertEqual(record['packet_ranks'], [4]*32)
            self.assertEqual(record['geometry'], dict(physical_steps=1, physical_step_bits=128,
                physical_windows=32, proof_step_bits=128, proof_windows=32))
            self.assertEqual(record['return_denominator'], (1 << bits)-1)
            self.assertFalse(record['full_state_census'])
            previous = rows

    def test_construction_pins(self):
        record = inner.construction(20)[2]
        self.assertEqual(set(record['source_sha256']), {str(inner.BASE_MAP),
            str(Path(inner.__file__).resolve()), str(Path(inner.s16_maps.__file__).resolve())})
        for path, digest in record['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(path).read_bytes()).hexdigest(), digest)

    def test_actual_small_census_preserves_high_bits(self):
        rows = (1 << 100, 1 << 7, (1 << 80) | (1 << 3))
        images, spectrum, dual = inner.enumerate_census(rows)
        self.assertEqual(images[1], 1 << 100)
        self.assertEqual(images[4], rows[2])
        self.assertEqual(images[7], rows[0] ^ rows[1] ^ rows[2])
        self.assertEqual(spectrum, {0: 1, 1: 2, 2: 2, 3: 2, 4: 1})
        self.assertEqual(len(dual), 129)
        self.assertEqual(sum(dual), 1 << 125)
        data = inner.kernel_maps.prepare_maps(images, inner.columns_from_rows(rows), bits=3,
            distribution='uniform_gl', birth_density='capped')
        self.assertEqual(data['windows'], 32)
        self.assertEqual(data['map_images'][1], 1 << 100)
        self.assertEqual(list(map(int, data['birth_state_weights'])), [0, 1, 1, 2, 2, 3, 3, 4])
        self.assertEqual(sum(map(int, data['histogram_multiplicities'])), 7)

    def test_invalid_dimension_or_density_precedes_census(self):
        for bits in (None, True, 16, 19, 21, 23, 20.0, '20'):
            with self.subTest(bits=bits), patch.object(inner, 'enumerate_census') as census:
                with self.assertRaises(ValueError):
                    inner.prepare(bits)
                census.assert_not_called()
        with patch.object(inner, 'enumerate_census') as census, self.assertRaises(ValueError):
            inner.prepare(20, 'unknown')
        census.assert_not_called()

    def test_prepare_dispatches_once_without_wrapping(self):
        for bits in inner.SUPPORTED_BITS:
            rows, columns, data, record = fixture(bits)
            images = ('mock-full-census',)
            spectrum = {int(w): n for w, n in record['expansion_spectrum'].items()}
            dual = tuple(map(int, record['feedback_kernel_spectrum']))
            with (patch.object(inner, 'enumerate_census', return_value=(images, spectrum, dual)) as census,
                  patch.object(inner.kernel_maps, 'prepare_maps', return_value=data) as prepare,
                  patch.object(inner.kernel_maps, 'authenticate'),
                  patch.object(inner.q1.kernel_t64, 'wrap', side_effect=AssertionError('no macro wrap'))):
                got, result = inner.prepare(bits)
            self.assertIs(got, data)
            census.assert_called_once_with(rows)
            prepare.assert_called_once_with(images, columns, bits=bits,
                distribution='uniform_gl', birth_density='capped')
            self.assertEqual(result['state_count'], 1 << bits)
            self.assertEqual(result['geometry']['physical_steps'], 1)
            self.assertFalse(result['whole_code_certificate'])

    def test_authentication_rejects_changed_geometry_rows_and_record(self):
        _, _, data, record = fixture(20)
        with patch.object(inner.kernel_maps, 'authenticate'):
            self.assertIs(inner.authenticate(data, record), data)
            for key, value in (('windows', 16), ('distribution', 'two_independent_uniform_gl_physical_steps'),
                               ('map_sha256', 'wrong'), ('profile_partition', True)):
                with self.subTest(key=key), self.assertRaises(ValueError):
                    inner.authenticate(dict(data, **{key: value}), record)
            wrong = deepcopy(data)
            wrong['map_images'][1] ^= 1 << 100
            with self.assertRaises(ValueError):
                inner.authenticate(wrong, record)
            for key, value in (('full_state_census', False), ('state_count', 1 << 19),
                               ('source_sha256', {}), ('minimum_expansion_weight', 16),
                               ('geometry', {'physical_steps': 2})):
                with self.subTest(key=key), self.assertRaises(ValueError):
                    inner.authenticate(data, dict(record, **{key: value}))

    def test_local_dispatch_is_one_physical_operator_family(self):
        raw, capped = [object() for _ in range(33)], [object() for _ in range(33)]
        for density in ('classes', 'capped'):
            _, _, data, _ = fixture(20, density)
            with (patch.object(inner.kernel_maps, 'authenticate'),
                  patch.object(inner.sparse_kernel, 'outward_at_z', return_value=raw) as sparse,
                  patch.object(inner.kernel_birth_density, 'refine_local', return_value=capped) as refine,
                  patch.object(inner.q1.kernel_t64, 'convolve', side_effect=AssertionError('no convolution'))):
                result = inner.local_operators(data, Q(1, 16), Q(1, 3))
            sparse.assert_called_once()
            if density == 'capped':
                refine.assert_called_once()
                self.assertEqual(refine.call_args.args[3], Q(1, 3))
                self.assertIs(result, capped)
            else:
                refine.assert_not_called()
                self.assertIs(result, raw)


if __name__ == '__main__':
    unittest.main()

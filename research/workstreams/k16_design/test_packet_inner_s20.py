"""Light adapter tests; the million-state prepare() census is deliberately not run."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import packet_inner_s20 as inner


class S20AdapterTests(unittest.TestCase):
    def test_selected_source_bytes_and_base_rows_preserved(self):
        before = inner.BASE_MAP.read_bytes()
        rows, columns, record = inner.construction()
        self.assertEqual(hashlib.sha256(before).hexdigest(), inner.BASE_SHA256)
        source = json.loads(before)
        self.assertEqual(rows[:16], tuple(int(word, 16) for word in source['generator_rows_hex']))
        self.assertEqual(tuple(column & 0xffff for column in columns), tuple(source['columns']))
        self.assertEqual(inner.BASE_MAP.read_bytes(), before)
        self.assertFalse(record['full_state_census'])
        self.assertFalse(record['whole_code_certificate'])

    def test_declared_extension_structure(self):
        rows, columns, record = inner.construction()
        self.assertEqual((record['t'], record['s']), (64, 20))
        self.assertEqual(record['appended_monomials'], [[0, 1], [0, 2], [0, 3], [0, 4]])
        self.assertEqual(rows[16:], tuple(inner.quadratic_monomial(0, j) for j in range(1, 5)))
        self.assertEqual(inner.s16_maps.binary_rank(rows), 20)
        self.assertEqual(inner.s16_maps.binary_rank(columns), 20)
        self.assertEqual(len(set(columns)), 64)
        self.assertNotIn(0, columns)
        self.assertTrue(all((a & b).bit_count() % 2 == 0 for a in rows for b in rows))
        self.assertEqual(record['packet_ranks'], [4] * 16)
        self.assertIn('GL20', record['sampling'])
        self.assertEqual(record['macro']['physical_steps'], 2)
        self.assertFalse(record['macro']['flush'])

    def test_constructor_and_base_pins_are_current(self):
        _, _, record = inner.construction()
        self.assertEqual(set(record['source_sha256']), {str(inner.BASE_MAP), str(Path(inner.__file__).resolve())})
        for path, digest in record['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(path).read_bytes()).hexdigest(), digest)
        self.assertEqual(record['source']['sha256'], inner.BASE_SHA256)
        self.assertNotEqual(record['source']['path'], record['constructor_source']['path'])

    def test_greedy_rank_rule_skips_dependent_candidates(self):
        linear = tuple(sum(1 << x for x in range(8) if (x >> bit) & 1) for bit in range(3))
        first = inner.quadratic_monomial(0, 1, variables=3)
        rows, pairs = inner.greedy_extension((255, *linear, first), 7, variables=3)
        self.assertEqual(pairs, ((0, 2), (1, 2)))
        self.assertEqual(inner.s16_maps.binary_rank(rows), 7)
        with self.assertRaises(ValueError):
            inner.greedy_extension((255, *linear), 8, variables=3)

    def test_small_reference_enumeration_and_dual_spectrum(self):
        rows = (0b00110011, 0b01010101, 0b11110000)
        images, spectrum, dual = inner.enumerate_census(rows, 8)
        reference = []
        for state in range(8):
            image = 0
            for j, row in enumerate(rows):
                if state >> j & 1:
                    image ^= row
            reference.append(image)
        self.assertEqual(images, tuple(reference))
        self.assertEqual(sum(spectrum.values()), 8)
        actual_dual = [0] * 9
        for word in range(256):
            if all((word & row).bit_count() % 2 == 0 for row in rows):
                actual_dual[word.bit_count()] += 1
        self.assertEqual(dual, tuple(actual_dual))
        self.assertEqual(sum(dual), 2**5)
        with self.assertRaises(ValueError):
            inner.enumerate_census((1, 1), 8)

    def test_identity_matches_retained_identity_on_small_maps(self):
        rows = (3, 5, 9)
        columns = inner.columns_from_rows(rows, 4)
        images, _, _ = inner.enumerate_census(rows, 4)
        self.assertEqual(inner._identity(rows, columns), inner.kernel_maps._identity(images, columns, 3))

    def test_generic_two_step_wrapper_accepts_non16_state_parameter(self):
        # A small real preparation checks reuse without a million-state census.
        rows = (3, 5, 9)
        images, _, _ = inner.enumerate_census(rows, 4)
        physical = inner.kernel_maps.prepare_maps(images, inner.columns_from_rows(rows, 4),
            bits=3, distribution='uniform_gl', birth_density='classes')
        wrapped = inner.kernel_t64.wrap(physical)
        self.assertEqual(wrapped['bits'], 3)
        self.assertEqual(wrapped['macro_windows'], 2)
        self.assertEqual(len(inner.kernel_t64.local_operators(wrapped, Q(1, 10))), 3)
        with self.assertRaises(ValueError):
            inner.authenticate(wrapped)

    def test_prepare_forwards_twenty_bits_without_running_large_census(self):
        rows, columns, record = inner.construction()
        fake_images = ('not-enumerated-in-this-control-flow-test',)
        fake_physical = dict(marker='physical')
        fake_wrapper = dict(bits=20, physical_step_bits=64, map_sha256=record['map_sha256'])
        fake_dual = [0] * 65
        fake_dual[0] = fake_dual[6] = 1
        for density in ('classes', 'capped'):
            with (patch.object(inner, 'enumerate_census', return_value=(fake_images, {0: 1, 16: 2**20-1}, fake_dual)),
                  patch.object(inner.kernel_maps, 'prepare_maps', return_value=fake_physical) as prepare,
                  patch.object(inner.kernel_t64, 'wrap', return_value=fake_wrapper)):
                wrapped, result = inner.prepare(density)
            prepare.assert_called_once_with(fake_images, columns, bits=20,
                distribution='uniform_gl', birth_density=density)
            self.assertEqual(wrapped['bits'], 20)
            self.assertEqual(result['s'], 20)
            self.assertEqual(result['state_count'], 2**20)
            self.assertTrue(result['full_state_census'])
            self.assertEqual(result['minimum_feedback_kernel_weight'], 6)
        with self.assertRaises(ValueError):
            inner.prepare('unsupported')


if __name__ == '__main__':
    unittest.main()

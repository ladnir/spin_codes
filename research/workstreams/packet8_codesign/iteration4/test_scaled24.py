"""Small exact algebra/profile checks; no large-state gate or benchmark."""
import unittest
import numpy as np
import scaled24 as s


class ScaledMaps(unittest.TestCase):
    def test_all_one_literal_identity(self):
        old = s.maps24.make_maps()
        new = s.make_scaled((1,)*8)
        self.assertEqual(old['record']['map_sha256'], new['record']['map_sha256'])
        self.assertEqual(old['rows'], new['rows'])
        for key in ('tables', 'transpose_tables'):
            np.testing.assert_array_equal(old[key], new[key])

    def test_only_expansion_changes(self):
        old = s.maps24.make_maps()
        for name in ('bands17', 'bands2e', 'power3'):
            new = s.make_scaled(s.candidate_scales(name))
            self.assertEqual(old['columns'], new['columns'])
            np.testing.assert_array_equal(old['transpose_tables'], new['transpose_tables'])
            self.assertEqual(new['record']['expansion_rank'], 24)
            self.assertEqual(new['record']['CA_columns'], [s.maps24.apply(new['columns'], r) for r in new['rows']])

    def test_literal_field_formula_and_binary_adjoints(self):
        data = s.make_scaled(s.candidate_scales('bands17'))
        scales = data['record']['output_byte_scales']
        for state in [0, 1, 0x123456, 0xffffff]+[1 << i for i in range(24)]:
            a, b, c = state & 255, (state >> 8) & 255, state >> 16
            expected = sum(s.maps24.multiply(scales[h], a ^ s.maps24.multiply(h, b)
                ^ s.maps24.multiply(s.maps24.multiply(h, h), c)) << (8*h) for h in range(8))
            self.assertEqual(s.maps24.apply(data['rows'], state), expected)
            for word in (0, 1, 0x123456789abcdef0, 0xffffffffffffffff):
                atrans = 0
                ctrans = 0
                for h in range(8):
                    value = (word >> (8*h)) & 255
                    for k in range(3):
                        atrans ^= int(data['expansion_transpose_tables'][h, k, value]) << (8*k)
                    byte = 0
                    for k in range(3):
                        byte ^= int(data['transpose_tables'][h, k, (state >> (8*k)) & 255])
                    ctrans |= byte << (8*h)
                self.assertEqual((expected & word).bit_count() & 1, (state & atrans).bit_count() & 1)
                self.assertEqual((s.maps24.apply(data['columns'], word) & state).bit_count() & 1,
                                 (word & ctrans).bit_count() & 1)

    def test_small_field_profiles_and_local_regression(self):
        old = s.maps24.make_maps(packet_bits=2, modulus=7, windows=4)
        for scales in ((1,)*4, (1, 1, 2, 2)):
            data = s.make_scaled(scales, packet_bits=2, modulus=7, windows=4)
            census = s.maps24.census(data, chunk_bits=3)
            for state in range(64):
                image = s.maps24.apply(data['rows'], state)
                weights = [(image >> (2*h) & 3).bit_count() for h in range(4)]
                expected = np.bincount(weights, minlength=3)
                np.testing.assert_array_equal(data['profiles'][census['expansion_indices'][state]], expected)
                charimage = sum(sum(((int(col) & state).bit_count() & 1) << i
                    for i, col in enumerate(data['columns'][2*h:2*h+2])) << (2*h) for h in range(4))
                expectedchar = np.bincount([(charimage >> (2*h) & 3).bit_count() for h in range(4)], minlength=3)
                np.testing.assert_array_equal(data['profiles'][census['character_indices'][state]], expectedchar)
            if scales == (1,)*4:
                oldlocal, _ = s.screen24.local_operators(old, s.maps24.census(old), .9)
                newlocal, _ = s.screen24.local_operators(data, census, .9)
                np.testing.assert_array_equal(oldlocal, newlocal)

    def test_proxy_and_invalid_scalings(self):
        for d in (0x17, 0x2e):
            self.assertEqual(s.scale_score(d)['minimum'], 4)
            self.assertEqual(sum(s.scale_score(d)['histogram']), 255)
        for scales in ((0,)*8, (1,)*7, (256,)*8, (True,)*8):
            with self.assertRaises(ValueError):
                s.make_scaled(scales)


if __name__ == '__main__':
    unittest.main()

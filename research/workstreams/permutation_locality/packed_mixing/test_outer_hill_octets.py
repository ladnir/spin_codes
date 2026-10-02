from itertools import combinations
import unittest

import outer_hill_octets as octets


def span(rows):
    words = [0]
    for row in rows:
        words += [word ^ row for word in words]
    return set(words)


class CanonicalOctetTests(unittest.TestCase):
    def test_quotient_kernel_matches_small_codes_exhaustively(self):
        for rows, length in (([0x3F, 0xFC0], 12), ([15, 51, 85], 8), ([3, 12], 4)):
            columns, rank = octets.syndrome_columns(rows, length)
            words = span(rows)
            self.assertEqual(len(words), 1 << rank)
            for word in range(1 << length):
                syndrome = 0
                for i, column in enumerate(columns):
                    if word >> i & 1:
                        syndrome ^= column
                self.assertEqual(syndrome == 0, word in words)

    def test_near_full_scan_matches_all_words_for_small_even_code(self):
        rows = [0x3F, 0xFC0]
        columns, _ = octets.syndrome_columns(rows, 12)
        result = octets.near_full_words(columns, 4, 2)
        exact = [word for word in span(rows) if sum(bool(word & (15 << (4*b))) for b in range(3)) == 2]
        self.assertEqual(set(result), set(exact))
        self.assertEqual(len(result), 2)

    def test_scan_includes_no_deletion_and_checks_pair_is_inside_blocks(self):
        rows = [0xFF, 0xFFF000]
        columns, _ = octets.syndrome_columns(rows, 24)
        result = octets.near_full_words(columns, 4, 2)
        self.assertEqual(result, [0xFF])

    def test_weight_two_or_four_and_bad_geometry_fail(self):
        for rows in ([3], [15]):
            columns, _ = octets.syndrome_columns(rows, 8)
            with self.assertRaises(ValueError):
                octets.near_full_words(columns, 4, 2)
        columns, _ = octets.syndrome_columns([0xFF], 8)
        for width, blocks in ((2, 2), (3, 1), (4, 0), (4, 3), (True, 1)):
            with self.assertRaises(ValueError):
                octets.near_full_words(columns, width, blocks)

    def test_binary_polynomial_remainder(self):
        self.assertEqual(octets.remainder(0b10101, 0b111), 0)
        self.assertEqual(octets.remainder(0b10100, 0b111), 1)

    def test_distance_check_rejects_wrong_dimension(self):
        with self.assertRaises(ArithmeticError):
            octets.verify_distance38([0xFF])


if __name__ == '__main__':
    unittest.main()

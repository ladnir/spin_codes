from fractions import Fraction as F
import unittest

import search_scaled_split as search


class ScaledSearchTests(unittest.TestCase):
    def test_rle_exact_case_coverage(self):
        powers = [[-80, -80, -79], [-80, -80, -60, -80]]
        table = search.fast.unpack(powers, 2, 3)
        encoded = search.encode_table(table, 2, 3)
        decoded = search.decode_table(encoded, 2, 3, -80, 100)
        self.assertEqual(search.fast.pack(decoded, 2, 3), powers)
        encoded[0][0][1] -= 1
        with self.assertRaises(ValueError):
            search.decode_table(encoded, 2, 3, -80, 100)

    def test_negative_and_invalid_owner_rejected(self):
        for value in (-1, 20):
            with self.assertRaises(ValueError):
                search.decode_table([[[value, 3]]], 2, 2, 0, 19)

    def test_passing_equals_exact_case_sum(self):
        powers = [[-80, -80, -79], [-80, -80, -50, -80]]
        result = search.passing(search.fast.unpack(powers, 2, 3), 2, 3)
        self.assertEqual(result, {2: F(4, 1 << 80)})

    def test_size_scaled_proposal(self):
        powers = [[-80]*965]
        powers[0][0] = 10
        proposal = search.proposal(search.fast.unpack(powers, 964, 964), 964, 964, [], 4)
        self.assertEqual(proposal, (F(-5, 2), 964, 0))


if __name__ == '__main__':
    unittest.main()

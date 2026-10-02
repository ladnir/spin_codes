"""Small exact structural tests for the isolated S16 candidate maps."""
import unittest
import s16_maps


class MapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.images, cls.columns, cls.record = s16_maps.candidate('weight5', 0)

    def test_expansion_is_injective_with_expected_spectrum(self):
        self.assertEqual(len(set(self.images)), 1 << 16)
        self.assertEqual(self.record['expansion_spectrum'],
            {0: 1, 48: 622, 56: 13840, 64: 36623, 72: 13808, 80: 642})
        self.assertEqual(self.images[0], 0)

    def test_feedback_has_no_short_dependency(self):
        self.assertEqual(s16_maps.binary_rank(self.columns), 16)
        self.assertEqual(len(set(self.columns)), 128)
        self.assertTrue(all(column.bit_count() == 5 for column in self.columns))
        self.assertEqual(self.record['feedback_kernel_spectrum'][:4], ['1', '0', '0', '0'])
        self.assertEqual(self.record['packet_ranks'], [4]*32)
        self.assertEqual(sum(map(int, self.record['feedback_kernel_spectrum'])), 1 << 112)

    def test_sampling_is_reproducible(self):
        _, columns, record = s16_maps.candidate('weight5', 0)
        self.assertEqual(columns, self.columns)
        self.assertEqual(record, self.record)

    def test_declared_bch_feedback_is_full_rank(self):
        _, columns, record = s16_maps.candidate('bch16')
        self.assertEqual(s16_maps.binary_rank(columns), 16)
        self.assertEqual(record['feedback_kernel_spectrum'][:6], ['1', '0', '0', '0', '0', '0'])
        self.assertEqual(record['minimum_feedback_kernel_weight'], 6)


if __name__ == '__main__':
    unittest.main()

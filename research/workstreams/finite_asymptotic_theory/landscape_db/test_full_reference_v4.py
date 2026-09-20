"""Occupation ledger rejection checks for full reference replay."""
import unittest

from verify_bch_full_reference_v4 import check_intervals


class IntervalCoverage(unittest.TestCase):
    def test_small_and_dense_complete_ranges(self):
        check_intervals([(5, 64)], 64, 257)
        check_intervals([(5, 128), (129, 256)], 256, 257)
        check_intervals([(5, 256), (257, 1024), (1025, 32768)], 32768, 1025)

    def test_rejects_missing_overlapping_and_wrong_dense_boundary(self):
        for intervals, length, minimum in [([], 64, 257), ([(6, 64)], 64, 257),
                ([(5, 63)], 64, 257), ([(5, 32), (32, 64)], 64, 257),
                ([(5, 32), (34, 64)], 64, 257), ([(5, 65)], 64, 257),
                ([(5, 4), (5, 64)], 64, 257),
                ([(5, 255), (256, 1024)], 1024, 257)]:
            with self.subTest(intervals=intervals):
                with self.assertRaises(ValueError):
                    check_intervals(intervals, length, minimum)


if __name__ == '__main__':
    unittest.main()

import unittest
from fractions import Fraction as F
from math import log

from hd_blockwise_maps import matmul, power
from hd_blockwise_shared_maps import identity, inverse, rank
from hd_randomized_outer_fft16_probe import evaluation, matrices, moments, directed_endpoint


class FFT16ProbeTests(unittest.TestCase):
    def test_algebraic_identities(self):
        feedthrough, injection, feedback, translation = matrices(16)
        self.assertEqual(matmul(feedthrough, feedthrough), identity(64))
        self.assertEqual(matmul(feedthrough, injection), matmul(injection, translation))
        self.assertEqual(matmul(feedback, feedthrough), matmul(translation, feedback))
        self.assertFalse(any(any(row) for row in matmul(feedback, injection)))
        self.assertEqual(matmul(translation, translation), identity(4))

    def test_global_mds_representation(self):
        _, injection, feedback, _ = matrices(16)
        vandermonde = tuple(tuple(power(x, j) for x in range(64)) for j in range(4))
        change = matmul(tuple(tuple(row[:4]) for row in feedback),
                       inverse(tuple(tuple(row[:4]) for row in vandermonde)))
        self.assertEqual(rank(change), 4)
        self.assertEqual(matmul(change, vandermonde), feedback)
        self.assertEqual(injection, tuple(row[:4] for offset in range(64, 128, 16)
                                        for row in evaluation(offset, 16)))

    def test_tile_residual_ranks(self):
        feedthrough, _, feedback, _ = matrices(16)
        checked = 0
        for tile in range(4):
            for a in range(17):
                columns = tuple(16*tile+(i+3*tile) % 16 for i in range(a))
                for v in range(17):
                    rows = tuple(16*tile+(i+5*tile+3) % 16 for i in range(v))
                    output = tuple(tuple(feedthrough[i][j] for j in columns) for i in rows)
                    combined = output+tuple(tuple(row[j] for j in columns) for row in feedback)
                    self.assertEqual(rank(output), min(a, v))
                    self.assertEqual(rank(combined)-rank(output), min(4, max(a-v, 0)))
                    checked += 1
        self.assertEqual(checked, 1156)

    def test_distinct_residual_lines(self):
        feedthrough, _, feedback, _ = matrices(16)
        for left in range(4):
            for right in range(left+1, 4):
                for a in range(1, 17):
                    columns = tuple(16*t+(i+3*t) % 16 for t in (left, right) for i in range(a))
                    rows = tuple(16*t+(i+5*t+3) % 16 for t in (left, right) for i in range(a-1))
                    output = tuple(tuple(feedthrough[i][j] for j in columns) for i in rows)
                    combined = output+tuple(tuple(row[j] for j in columns) for row in feedback)
                    self.assertEqual(rank(combined)-rank(output), 2)

    def test_width8_regression_and_directed_arithmetic(self):
        from flint import ctx
        from hd_randomized_outer_joint_quiet import joint_moment_logs
        ctx.prec = 160
        for mark in (F(1), F(1, 2), F(1, 1024)):
            self.assertEqual(moments(mark, 8), joint_moment_logs(mark))
            logs = moments(mark, 16)
            exact = moments(mark, 16, directed=True)
            for log_row, exact_row in zip(logs, exact):
                for approximate, enclosure in zip(log_row, exact_row):
                    self.assertAlmostEqual(approximate, float(enclosure.log())/log(2), places=10)
        receipt = directed_endpoint(16, 1024, F(1, 2), F(2547111, 16777216), F(1, 16))
        self.assertFalse(receipt['below_2_to_minus_40'])
        self.assertTrue(receipt['failure_log2_enclosure'].startswith('[3855.3021455243778'))


if __name__ == '__main__':
    unittest.main()

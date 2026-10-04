"""Small arithmetic and geometry regressions; no census or large placement."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx

import outward_screen as screen


class OutwardScreenTests(unittest.TestCase):
    def setUp(self):
        self.old_precision = ctx.prec
        ctx.prec = 256
        self.regional = arb_mat([[screen.replay.aq(F(1, 4)), screen.replay.aq(F(1, 16))],
                                 [screen.replay.aq(F(1, 32)), screen.replay.aq(F(1, 2))]])

    def tearDown(self):
        ctx.prec = self.old_precision

    def test_k18_scalar_matches_frozen_replay_exactly(self):
        for q in (1, 89, 2048):
            for alpha in (F(4, 5), F(1)):
                with self.subTest(q=q, alpha=alpha):
                    expected = screen.replay.scalar_bound(self.regional, q, F(17, 200), alpha)
                    actual = screen.scalar(self.regional, q, F(17, 200), alpha, 1 << 18)
                    self.assertEqual(screen.replay.endpoint(actual), screen.replay.endpoint(expected))

    def test_k20_changes_binomial_and_floor_cutoff_outside_power(self):
        # With the same regional matrix the beta and moment cancel in this ratio.
        # The integer group-subset factor is not raised to the fractional power.
        cutoff_delta = 209715 - 52428
        for q in (1, 128, 512):
            for alpha in (F(4, 5), F(1)):
                with self.subTest(q=q, alpha=alpha):
                    tilt = F(17, 200)
                    old = screen.replay.scalar_bound(self.regional, q, tilt, alpha)
                    scale = (screen.replay.aq(F(comb(8192, q), comb(2048, q))) *
                             screen.replay.aq(alpha * tilt * cutoff_delta).exp())
                    expected = old * scale
                    actual = screen.scalar(self.regional, q, tilt, alpha, 1 << 20)
                    relative_error = abs(actual / expected - 1)
                    self.assertTrue(relative_error < screen.replay.aq(F(1, 1 << 230)))

    def test_invalid_geometry_rejected_before_loading_local_receipt(self):
        # Nonintegral groups, incomplete 32-slot macros, nonpositive length,
        # and out-of-range occupancy must not reach the numerical engine.
        invalid = ((1048577, [1]), (128, [1]), (0, [1]), (-4096, [1]),
                   (1 << 20, []), (1 << 20, [0]), (1 << 20, [8193]))
        with patch.object(Path, 'exists', return_value=False), \
                patch.object(screen.replay, 'load_local') as load:
            for k, qs in invalid:
                with self.subTest(k=k, qs=qs), self.assertRaises(ValueError):
                    screen.run(Path('unused-local.json'), Path('unused-output.json'), ['.8'], qs, k=k)
            load.assert_not_called()


if __name__ == '__main__':
    unittest.main()

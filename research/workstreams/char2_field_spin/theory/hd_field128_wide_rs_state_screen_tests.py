"""Small algebra and normalization tests; no selected-weight distance gate."""

from fractions import Fraction as F
from itertools import combinations
from math import comb
import unittest

import numpy as np

import hd_field128_wide_rs_state_screen as wide
import hd_field128_full_rs_state_screen as frozen
from hd_field128_full_rs_screen_tests import toy_moments
from hd_field128_full_rs_state_screen_tests import rank
from hd_blockwise_maps import matmul


class WideStateTests(unittest.TestCase):
    def test_actual_wide_maps_ca_zero_and_mds_minor_samples(self):
        a, c, h = wide.fixed_maps()
        self.assertEqual((len(a), len(a[0]), len(h)), (128, 64, 64))
        self.assertFalse(any(value for row in matmul(c, a, 256) for value in row))
        self.assertEqual(rank(a), 64)
        self.assertEqual(rank(h), 64)
        for size in (1, 2, 3, 8, 32, 64):
            selected = tuple(tuple(h[i][j] for j in range(64-size, 64)) for i in range(size))
            self.assertEqual(rank(selected), size)

    def test_toy_mds_ranks_and_all_refresh_minors(self):
        for state in (2, 3, 4):
            a = frozen.lagrange(range(8), range(8, 8+state), 16)
            self.assertFalse(any(value for row in matmul(tuple(zip(*a)), a, 16) for value in row))
            for size in range(1, 9):
                for selected in combinations(range(8), size):
                    self.assertEqual(rank(tuple(a[i] for i in selected), 16), min(size, state))
            h = frozen.lagrange(range(state, 2*state), range(state), 16)
            for size in range(1, state+1):
                for rows in combinations(range(state), size):
                    for columns in combinations(range(state), size):
                        self.assertEqual(rank(tuple(tuple(h[i][j] for j in columns) for i in rows), 16), size)

    def test_principal_free_return_count_bound_by_exact_toy_ranks(self):
        q, n, state = 16, 8, 3
        a = frozen.lagrange(range(n), range(n, n+state), q)
        c = tuple(zip(*a))
        ac = matmul(a, c, q)
        for positions in ((0,), (0, 3), (0, 2, 5, 7)):
            j, rho = len(positions), min(state, len(positions))
            inactive = tuple(i for i in range(n) if i not in positions)
            syndrome_rows = tuple(tuple(row[i] for i in positions) for row in c)
            for cancelled in range(j+1):
                for inside in combinations(positions, cancelled):
                    for e in range(min(state, len(inactive))+1):
                        for outside in combinations(inactive, e):
                            equations = tuple(tuple(ac[i][k]^(i == k) for k in positions)
                                              for i in inside+outside)
                            # Subtract the syndrome-zero subspace from the zero-output subspace.
                            actual = q**(j-rank(equations, q))-q**(j-rank(equations+syndrome_rows, q))
                            upper = min(q**j-q**(j-rho),
                                        (q**min(rho, max(0, state-e))-1)*q**max(0, j-cancelled-state))
                            self.assertLessEqual(actual, upper)

    def test_exact_toy_moments_and_positive_float_expansions(self):
        z, q, n = F(7, 8), 16, 8
        for state in (2, 3, 4):
            a = frozen.lagrange(range(n), range(n, n+state), q)
            for positions in ((), (0,), (0, 1), (1, 6)):
                actual_quiet, actual_return = toy_moments(a, positions, z=z, q=q)
                quiet, birth, cap, live = wide.exact_moments(len(positions), z, q, n, state)
                self.assertEqual(actual_quiet, quiet)
                self.assertLessEqual(actual_return, cap)
                self.assertLessEqual(actual_return, birth)
                self.assertLessEqual(actual_return, (q**state-1)*live)
            logs = wide.local_logs(float(z), q, n, state).reshape((n+1, 4))
            for j in range(n+1):
                quiet, birth, cap, live = wide.exact_moments(j, z, q, n, state)
                self.assertEqual(birth, (z+(1-z)/q)**j-quiet)
                inactive = n-j
                flat = sum((F(comb(inactive, e), q**min(state, e))*z**(inactive-e)*(1-z)**e
                            for e in range(inactive+1)), F(0))
                self.assertEqual(live, (z+(1-z)/q)**j*(q**state*flat-1)/(q**state-1))
                self.assertLessEqual(cap, wide.exact_moments(j, z, q, n, state, False)[2])
                np.testing.assert_allclose(logs[j], [wide.old.log_fraction(value)
                                                    for value in wide.exact_packet(j, z, q, n, state)],
                                           rtol=0, atol=1e-11)

    def test_width64_replays_frozen_oracles_and_regional_prefix(self):
        for state in frozen.STATES:
            self.assertEqual(wide.fixed_maps(64, state), frozen.fixed_maps(state))
            self.assertEqual(wide.construction(64, state), frozen.construction(state))
            for joint in (False, True):
                np.testing.assert_array_equal(wide.local_logs(.5, 256, 64, state, joint),
                                              frozen.local_logs(.5, 256, 64, state, joint))
                for j in (0, 1, 16, 33, 64):
                    self.assertEqual(wide.exact_packet(j, F(1, 2), 256, 64, state, joint),
                                     frozen.exact_packet(j, F(1, 2), 256, 64, state, joint))
            self.assertEqual(wide.regional_screen(.5, .1, 3, q=256, width=64, state=state),
                             frozen.regional_screen(.5, .1, 3, q=256, state=state))
        self.assertEqual(wide.geometry(64, width=64), frozen.fixed16.geometry(64))

    def test_pair_geometry_and_exact_regional_normalizer(self):
        maximum = 6
        for width in wide.WIDTHS:
            pairs = width//2
            geometry = wide.geometry(maximum, width=width)
            for s, f, count in geometry:
                # Two orientations per singleton are averaged by fair independent swaps.
                self.assertEqual(F(comb(pairs, f)*comb(pairs-f, s)*2**s, 2**s), count)
            # Scalar convolution independently checks the whole-region multinomial.
            counts = {(0, 0): 1}
            for unused in range(wide.routing(width)['packets_per_region']):
                next_counts = {}
                for (s, f), count in counts.items():
                    for t, g, other in geometry:
                        if s+t+2*(f+g) <= maximum:
                            key = s+t, f+g
                            next_counts[key] = next_counts.get(key, 0)+count*other
                counts = next_counts
            self.assertEqual(counts, {(s, f): comb(wide.old.PAIRS, f)*comb(wide.old.PAIRS-f, s)
                                      for f in range(maximum//2+1) for s in range(maximum-2*f+1)})
        self.assertEqual(len(wide.geometry(128, good=True)), 129)

    def test_alpha1_good_pair_composition_and_regional_replay(self):
        maximum, z, q = 3, .5, 256
        budget = wide.old.Budget(60.)
        single = wide.single_coefficients(z, 1., maximum, q, 128, 64, budget)
        single[:, 2:, :] = wide.NEG
        composed = wide.good_pair_coefficients(z, 1., maximum, q, 128, 64, budget)
        product = wide.old.polynomial_product(single, single, maximum, budget)
        np.testing.assert_allclose(composed, product, rtol=0, atol=2e-12)
        self.assertEqual(wide.regional_screen(z, 1., maximum, q=q, paired=True),
                         wide.regional_screen(z, 1., maximum, q=q, paired=False))

    def test_q128_tiny_mark_retains_live_and_return_logs(self):
        z, q = F(1, 1 << 200), 1 << 128
        for joint in (False, True):
            logs = wide.local_logs(float(z), q, joint=joint)[1].ravel()
            exact = wide.exact_packet(1, z, q, joint=joint)
            np.testing.assert_allclose(logs, [wide.old.log_fraction(value) for value in exact],
                                       rtol=0, atol=2e-9)
            self.assertTrue(np.all(np.isfinite(logs)))
            self.assertLess(logs[2], -64*128)
            self.assertGreater(exact[2], 0)
            self.assertEqual(float(exact[2]), 0.)

    def test_metadata_validation_and_source_pins(self):
        details = wide.metadata()
        self.assertEqual(details['construction'],
                         'toeplitz-fullfield-full-rs128-64-alternating-swapped-pairs-k20-rate-half-v1')
        for key, value in dict(inner_width=128, state_dimension=64, emission_minimum=65,
                               pair_slots_per_packet=64, pair_slots_per_region=4096,
                               symbols_per_region=8192, packets_per_region=64,
                               composed_pairs_per_region=32, regions=256, outer_randomizers_per_group=64,
                               refresh_diagonal_fullfield_products_per_codeword=1 << 21,
                               refresh_diagonal_key_field_elements_per_codeword=1 << 21).items():
            self.assertEqual(details[key], value)
        self.assertEqual(details['refresh_diagonal_fullfield_products_per_codeword'],
                         wide.metadata(64, 32)['refresh_diagonal_fullfield_products_per_codeword'])
        for width, state in ((64, 33), (128, 65), (256, 64), (32, 16)):
            with self.assertRaises(ValueError):
                wide.construction(width, state)
        for maximum in (0, 65):
            with self.assertRaises(ValueError):
                wide.regional_screen(.5, .1, maximum)
        with self.assertRaises(ValueError):
            wide.exact_packet(1, .5)
        pins = wide.source_pins()
        self.assertIn('hd_field128_wide_rs_state_screen.py', pins)
        self.assertIn('hd_field128_full_rs_state_screen.py', pins)
        self.assertIn('hd_field128_lift_screen.py', pins)


if __name__ == '__main__':
    unittest.main()

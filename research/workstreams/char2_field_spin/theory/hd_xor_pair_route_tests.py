"""Exact geometry, fixed-support conditioning, and directed cap controls."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import comb, log2
import unittest

import numpy as np
from flint import arb, ctx

from hd_bucket_feedback import profiles
from hd_shifted_swap_sparse import encode_arb, power_arb, sparse_geometry
from hd_xor_multigroup import nonzero_arb_coefficients, nonzero_coefficients, region_operators
from hd_xor_pair_route import (arb_pair_coefficients, arb_pair_region_matrices,
    arb_pair_regions, feasible_full_pairs, nonzero_packet_caps_arb,
    pair_coefficients_screen, pair_region_operators, pair_regions_screen,
    route_metadata)


class PairRouteTests(unittest.TestCase):
    def test_packet_geometry_against_labeled_subsets(self):
        actual = Counter()
        # Pair j occupies the corresponding slot in two neighboring RS blocks.
        pairs = tuple((j, j+16) for j in range(16))+tuple((j, j+16) for j in range(32, 48))
        pair_of = {position: i for i, pair in enumerate(pairs) for position in pair}
        for m in range(4):
            for support in combinations(range(64), m):
                occupied = Counter(pair_of[position] for position in support)
                full = sum(value == 2 for value in occupied.values())
                occupancy = tuple(sorted(sum(position//16 == block for position in support)
                                         for block in range(4)))
                actual[m-2*full, full, occupancy] += 1
        expected = Counter({(s, f, occ): count for s, f, count, occ in sparse_geometry(3)})
        self.assertEqual(actual, expected)

    def test_fair_orientation_and_unsigned_normalizer(self):
        counts = Counter()
        for s, f, count, _ in sparse_geometry(17):
            counts[s, f] += count
        for (s, f), count in counts.items():
            self.assertEqual(F(count, 2**s), comb(32, s)*comb(32-s, f))
        for m in range(18):
            total = sum(comb(4096, m-2*f)*comb(4096-m+2*f, f)*2**(m-2*f)
                        for f in range(m//2+1))
            self.assertEqual(total, comb(8192, m))

    def test_subset_full_pair_range_is_exact(self):
        # Exhaust every source and every subset on a four-pair matching.
        for source_mask in range(256):
            source = tuple(i for i in range(8) if source_mask >> i & 1)
            r = len(source)
            full = sum((source_mask >> (2*j)) & 3 == 3 for j in range(4))
            for m in range(r+1):
                attainable = set()
                for subset in combinations(source, m):
                    mask = sum(1 << i for i in subset)
                    attainable.add(sum((mask >> (2*j)) & 3 == 3 for j in range(4)))
                self.assertEqual(attainable, set(feasible_full_pairs(m, r, full)))

    def test_profile_caps_reproduce_existing_directed_scalar_caps(self):
        previous = ctx.prec
        try:
            ctx.prec = 256
            z, maximum = F(799, 800), 8
            caps = nonzero_packet_caps_arb(z, maximum)
            actual = [[arb(0) for _ in range(maximum+1)] for _ in range(4)]
            for profile, count in profiles(16, 64, 4):
                if sum(profile) <= maximum:
                    for entry in range(4):
                        actual[entry][sum(profile)] += count*caps[profile][entry]
            expected = nonzero_arb_coefficients(z, maximum)
            for entry in range(4):
                for m in range(maximum+1):
                    difference = actual[entry][m]-expected[entry][m]
                    # Both are directed implementations of the same selected
                    # cap. Endpoint staging differs by rounding, not the law.
                    self.assertTrue(abs(difference) < arb(2)**-190*max(arb(1), abs(expected[entry][m])))
        finally:
            ctx.prec = previous

    def test_oriented_pair_mixture_reproduces_scalar_packet(self):
        z, maximum = .99875, 8
        paired = pair_coefficients_screen(z, maximum)
        expected = np.exp2(nonzero_coefficients(z, maximum=maximum))
        for m in range(maximum+1):
            actual = sum((2**(m-2*f)*paired[:, m-2*f, f] for f in range(m//2+1)), np.zeros(4))
            np.testing.assert_allclose(actual.reshape(2, 2), expected[m], rtol=3e-13, atol=0)

    def test_ordered_positive_coefficient_power(self):
        previous = ctx.prec
        try:
            ctx.prec = 256
            maximum = 4
            raw = {(0, 0): (1, 2, 3, 4), (1, 0): (3, 0, 2, 1), (0, 1): (0, 4, 1, 2)}
            table = {key: np.array(value, dtype=object).reshape(2, 2) for key, value in raw.items()}
            actual = {(0, 0): np.eye(2, dtype=object)}
            for _ in range(5):
                following = {}
                for (s, f), left in actual.items():
                    for (t, g), right in table.items():
                        if s+t+2*(f+g) <= maximum:
                            key = s+t, f+g
                            following[key] = following.get(key, np.zeros((2, 2), dtype=object))+left@right
                actual = following
            encoded = encode_arb({key: tuple(arb(x) for x in values) for key, values in raw.items()}, maximum)
            powered = power_arb(encoded, 5, maximum)
            for (s, f), exact in actual.items():
                for entry, value in enumerate(exact.ravel()):
                    self.assertTrue(powered[entry][s+(2*maximum+1)*f] == arb(value))
        finally:
            ctx.prec = previous

    def test_uniform_scalar_route_is_histogram_mixture(self):
        z, maximum = .99875, 6
        paired = pair_regions_screen(z, maximum)
        scalar = region_operators(z, maximum)
        envelope = pair_region_operators(z, maximum)
        for m in range(maximum+1):
            terms = []
            for f in range(m//2+1):
                s = m-2*f
                log_probability = (log2(comb(4096, s)*comb(4096-s, f))+s-log2(comb(8192, m)))
                terms.append(paired[m, f]+log_probability)
            mixture = np.logaddexp2.reduce(terms)
            np.testing.assert_allclose(mixture, scalar[m], rtol=0, atol=2e-11)
            self.assertTrue(np.all(mixture <= envelope[m]+2e-11))

    def test_full_directed_table_two_precisions_and_float_control(self):
        previous = ctx.prec
        try:
            z, maximum = F(799, 800), 17
            ceilings = []
            for precision in (256, 384):
                ctx.prec = precision
                regions = arb_pair_regions(z, maximum)
                self.assertEqual(len(regions), 90)
                row = []
                for key in sorted(regions):
                    for value in regions[key]:
                        self.assertTrue(value >= 0)
                        row.append(None if value.is_zero() else int((value.log()/arb(2).log()*10**10).upper().ceil().unique_fmpz()))
                ceilings.append(row)
            self.assertEqual(ceilings[0], ceilings[1])
            screens = pair_regions_screen(float(z), maximum)
            for key, entries in regions.items():
                expected = np.array([float(value.log()/arb(2).log()) if not value.is_zero() else -np.inf
                                     for value in entries]).reshape(2, 2)
                np.testing.assert_allclose(screens[key], expected, rtol=0, atol=2e-10)
        finally:
            ctx.prec = previous

    def test_restricted_envelope_and_metadata(self):
        previous = ctx.prec
        try:
            ctx.prec = 256
            z, maximum, source_full = F(799, 800), 4, 1
            regions = arb_pair_regions(z, maximum)
            caps = arb_pair_region_matrices(z, maximum, source_full_pairs=source_full)
            for m, row in enumerate(caps):
                for f in feasible_full_pairs(m, maximum, source_full):
                    self.assertTrue(all(upper >= lower for upper, lower in zip(row, regions[m, f])))
            self.assertEqual(caps[maximum], regions[maximum, source_full])
            metadata = route_metadata(maximum, directed=True)
            self.assertTrue(metadata['fixed_support_conditioning'])
            self.assertFalse(metadata['source_geometry_resampled'])
            self.assertFalse(metadata['scalar_route_receipts_transferred'])
            self.assertFalse(metadata['whole_code_certificate'])
        finally:
            ctx.prec = previous

    def test_guards_and_zero_occupancy(self):
        for maximum in (-1, 18, True, 1.5):
            with self.assertRaises(ValueError):
                pair_region_operators(.99, maximum)
        for z in (0, 1.1, float('nan')):
            with self.assertRaises(ValueError):
                pair_region_operators(z, 1)
        for arguments in ((3, 2, 0), (2, 4, 3), (2, None, 1)):
            with self.assertRaises(ValueError):
                feasible_full_pairs(*arguments)
        self.assertEqual(pair_region_operators(.99, 0).shape, (1, 2, 2))
        self.assertEqual(len(arb_pair_region_matrices(F(99, 100), 0)), 1)


if __name__ == '__main__':
    unittest.main()

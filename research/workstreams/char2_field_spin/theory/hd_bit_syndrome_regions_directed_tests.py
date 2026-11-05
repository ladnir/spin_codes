"""Tiny exact tests for the new directed original-one-state adapter."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import comb
from unittest import TestCase, main
from unittest.mock import patch

from flint import arb, arb_poly, ctx

import hd_bit_syndrome_regions_directed as regions
import hd_xor_pair_route as old
from hd_shifted_swap_sparse import sparse_geometry


def mm(left, right):
    return tuple(left[2*i]*right[j]+left[2*i+1]*right[2+j] for i in range(2) for j in range(2))


def exact_power(base, exponent, maximum, full_max):
    result = {(0, 0): (1, 0, 0, 1)}
    for _ in range(exponent):
        following = {}
        for (s, f), left in result.items():
            for (t, g), right in base.items():
                if s+t+2*(f+g) <= maximum and f+g <= full_max:
                    key = s+t, f+g
                    product = mm(left, right)
                    current = following.get(key, (0, 0, 0, 0))
                    following[key] = tuple(a+b for a, b in zip(current, product))
        result = following
    return result


class BitSyndromeRegionsDirectedTests(TestCase):
    def test_packet_caps_match_reviewed_interface_and_extend_without_mutation(self):
        with ctx.workprec(256):
            expected = old.nonzero_packet_caps_arb(regions.MARK, 6)
            actual = regions.packet_caps_arb(regions.MARK, 6)
            self.assertEqual(set(actual), set(expected))
            for profile in actual:
                for left, right in zip(actual[profile], expected[profile]):
                    self.assertEqual(left, right)
            selected = regions.packet_caps_arb(regions.MARK, 18,
                selected_profiles=((0, 0, 2, 16),))[(0, 0, 2, 16)]
            self.assertTrue(selected[0].is_zero())
            birth = regions.rational(F(256**16-1, 255**16)*regions.MARK**18)
            self.assertTrue(abs(selected[1]-birth) < arb(2)**-200)
            self.assertTrue(all(value >= 0 for value in selected))
        self.assertEqual(old.MAXIMUM, 17)
        with self.assertRaises(ValueError):
            old.nonzero_packet_caps_arb(regions.MARK, 18)

    def test_geometry_counts_swaps_once(self):
        actual = Counter()
        pairs = tuple((j, j+16) for j in range(16))+tuple((j, j+16) for j in range(32, 48))
        pair_of = {position: index for index, pair in enumerate(pairs) for position in pair}
        for weight in range(4):
            for support in combinations(range(64), weight):
                occupied = Counter(pair_of[i] for i in support)
                full = sum(count == 2 for count in occupied.values())
                profile = tuple(sorted(sum(i//16 == block for i in support) for block in range(4)))
                actual[weight-2*full, full, profile] += 1
        expected = Counter({(s, f, p): count for s, f, count, p in sparse_geometry(3)})
        self.assertEqual(actual, expected)
        table = {p: (arb(1), arb(2), arb(3), arb(4)) for _, _, _, p in sparse_geometry(3)}
        with ctx.workprec(256):
            caps = regions.packet_from_table(table, 3, 1)
            for (s, f), values in caps.items():
                multiplier = comb(32, s)*comb(32-s, f)
                self.assertEqual(values, tuple(arb(multiplier*x) for x in (1, 2, 3, 4)))

    def test_ordered_positive_power_matches_exact_dictionary(self):
        base = {(0, 0): (1, 2, 0, 1), (1, 0): (2, 1, 3, 0),
                (0, 1): (0, 1, 1, 2), (4, 0): (3, 0, 0, 4), (0, 2): (9, 1, 0, 8)}
        with ctx.workprec(256):
            wide = regions.power_arb(regions.encode_arb(base, 4, 2), 5, 4, 2)
            narrow_base = {key: value for key, value in base.items() if key[1] <= 1}
            narrow = regions.power_arb(regions.encode_arb(narrow_base, 4, 1), 5, 4, 1)
            expected = exact_power(base, 5, 4, 1)
            for s, f in regions._scope(4, 1):
                for entry in range(4):
                    self.assertEqual(narrow[entry][s+9*f], arb(expected.get((s, f), (0,)*4)[entry]))
                    self.assertEqual(narrow[entry][s+9*f], wide[entry][s+9*f])
            identity = regions.power_arb(regions.encode_arb(narrow_base, 4, 1), 0, 4, 1)
            self.assertEqual(tuple(value[0] for value in identity), (arb(1), arb(0), arb(0), arb(1)))
        self.assertNotEqual(mm(base[0, 0], base[1, 0]), mm(base[1, 0], base[0, 0]))

    def test_global_normalizer_and_marker_cancel_exactly(self):
        base = {(0, 0): (1, 1, 2, 1), (1, 0): (3, 1, 0, 2), (0, 1): (1, 0, 1, 1)}
        expected = exact_power(base, 3, 3, 1)
        with ctx.workprec(256):
            outputs = []
            for marker in (F(1), F(1, 8)):
                powered = regions.power_arb(regions.encode_arb(base, 3, 1, marker), 3, 3, 1)
                actual, normalizers = regions.normalize_arb(powered, 3, 1, pairs=6, marker=marker)
                for (m, f), values in actual.items():
                    s = m-2*f
                    count = comb(6, s)*comb(6-s, f)
                    self.assertEqual(normalizers[m, f], count)
                    for entry, value in enumerate(values):
                        exact = regions.rational(F(expected.get((s, f), (0,)*4)[entry], count))
                        self.assertTrue(abs(value-exact) < arb(2)**-200)
                outputs.append(actual)
            self.assertEqual(outputs[0], outputs[1])

    def test_exact_source_counts_and_production_histograms(self):
        counts = regions.source_counts()
        self.assertEqual(counts['F0'], 4076644197931783287105561490932293125919092974583472128)
        self.assertEqual(counts['F1'], 154355220469939768130077702840442839490778956275433472)
        self.assertEqual(counts['F_le1'], counts['F0']+counts['F1'])
        self.assertEqual(counts['F_le1']+counts['F_ge2'], comb(8192, 18))
        self.assertEqual(sum(counts['all_strata'].values()), counts['all_sources'])
        self.assertEqual(len(regions._scope(18, 1)), 36)
        self.assertIn((18, 0), regions._scope(18, 1))
        self.assertIn((16, 1), regions._scope(18, 1))

    def test_precisions_and_metadata_without_production_bank(self):
        zero_only = {(0, 0): (arb(1), arb(0), arb(0), arb(1))}
        with ctx.workprec(80):
            for precision in (256, 384):
                seen = []
                def fake_cap(*args, **kwargs):
                    seen.append(ctx.prec)
                    return {}
                with patch.object(regions, 'packet_caps_arb', side_effect=fake_cap), \
                        patch.object(regions, 'packet_from_table', return_value=zero_only):
                    result = regions.regional_bank(precision)
                self.assertEqual(ctx.prec, 80)
                self.assertEqual(seen, [precision])
                self.assertEqual(len(result['regions']), 36)
                metadata = result['metadata']
                self.assertEqual(metadata['mark'], '993/1000')
                self.assertEqual(metadata['precision_bits'], precision)
                self.assertEqual(metadata['global_normalization_count'], 1)
                self.assertEqual(metadata['carry_free_stride'], 37)
                self.assertTrue(metadata['independent_fair_swaps'])
                self.assertFalse(metadata['source_count_in_regional_entries'])
                self.assertFalse(metadata['post_refresh_bank_used'])
                self.assertFalse(metadata['maximum17_guard_modified'])
                self.assertFalse(metadata['whole_code_certificate'])
                self.assertFalse(metadata['distance_certificate'])

    def test_tiny_regional_replay_agrees_at_two_precisions(self):
        # Only degree2 and three packets: this is not the production bank.
        signatures = []
        for precision in (256, 384):
            with ctx.workprec(precision):
                caps = regions.packet_caps_arb(regions.MARK, 2)
                coefficients = regions.packet_from_table(caps, 2, 1)
                powered = regions.power_arb(regions.encode_arb(coefficients, 2, 1), 3, 2, 1)
                result, _ = regions.normalize_arb(powered, 2, 1, pairs=96)
                self.assertEqual(result[0, 0][0], arb(1))
                self.assertTrue(result[1, 0][0].is_zero())
                signatures.append(tuple(None if value.is_zero() else
                    int((value.log()/arb(2).log()*10**12).upper().ceil().unique_fmpz())
                    for key in sorted(result) for value in result[key]))
        self.assertEqual(signatures[0], signatures[1])

    def test_guards_and_precision_restoration(self):
        for precision in (128, 257, True):
            with self.assertRaises(ValueError):
                regions.regional_bank(precision)
        for mark in (.993, F(99, 100), F(0)):
            with self.assertRaises(ValueError):
                regions.regional_bank(mark=mark)
        with self.assertRaises(ValueError):
            regions.packet_caps_arb(F(1), 19)
        with self.assertRaises(ValueError):
            regions.encode_arb({(0, 2): (1, 0, 0, 1)}, 4, 1)
        bad = tuple(arb_poly([0, 0, 0, 1]) for _ in range(4))
        with self.assertRaises(ValueError):
            regions.product_arb(bad, bad, 2, 0)
        with ctx.workprec(80):
            with patch.object(regions, 'packet_caps_arb', side_effect=TimeoutError('tiny stop')):
                with self.assertRaises(TimeoutError):
                    regions.regional_bank(256)
            self.assertEqual(ctx.prec, 80)


if __name__ == '__main__':
    main()

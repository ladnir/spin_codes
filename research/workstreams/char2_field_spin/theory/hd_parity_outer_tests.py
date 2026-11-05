from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import unittest

from hd_parity_outer import (
    nonzero_sum_count, nonzero_nonone_sum_count, support_branches,
    support_counts, forced_support_law, direct, delta_full_pairs,
    pair_profile, profile_classes, primitive_counts)
from hd_outer_pair_layers import one_layer_support_profile, one_layer_weight_law


def xor(values):
    answer = 0
    for value in values:
        answer ^= value
    return answer


class ParityOuterTests(unittest.TestCase):
    def test_additive_counts_GF4_exact(self):
        for length in range(6):
            counts = Counter(xor(values) for values in product(range(1, 4), repeat=length))
            self.assertEqual(nonzero_sum_count(length, q=4), counts[0])
            self.assertEqual(nonzero_sum_count(length, target_zero=False, q=4), counts[1])
            punctured = Counter(xor(values) for values in product((2, 3), repeat=length))
            self.assertEqual(nonzero_nonone_sum_count(length, q=4), punctured[0])
            self.assertEqual(nonzero_nonone_sum_count(length, target_zero=False, q=4), punctured[1])

    def test_every_GF4_allowed_and_forced_support_in_four_coordinates(self):
        for mask in range(16):
            indices = [i for i in range(4) if mask >> i & 1]
            for mode in ('allowed', 'forced'):
                expected = {'quiet': Counter(), 'nonquiet': Counter()}
                domain = range(4) if mode == 'allowed' else range(1, 4)
                for values in product(domain, repeat=len(indices)):
                    source = [0]*4
                    for i, value in zip(indices, values):
                        source[i] = value
                    output = sum(bool(value) << i for i, value in enumerate(direct(source)))
                    expected['quiet' if xor(source) == 0 else 'nonquiet'][output] += 1
                self.assertEqual(support_branches(mask, mode=mode, q=4, size=4), expected)
                total = expected['quiet']+expected['nonquiet']
                self.assertEqual(support_counts(mask, mode=mode, q=4, size=4), total)

    def test_involution_zero_and_full_allowed_support(self):
        for values in product(range(4), repeat=4):
            self.assertEqual(direct(direct(values)), values)
        self.assertEqual(support_counts(0), {0: 1})
        counts = support_counts(15, q=4, size=4)
        self.assertEqual(counts, {mask: 3**mask.bit_count() for mask in range(16)})
        # Quiet and nonquiet can have the SAME support; do not overwrite.
        branches = support_branches(15, q=4, size=4)
        self.assertTrue(set(branches['quiet']).intersection(branches['nonquiet']))

    def test_allowed_support_is_not_forced_support(self):
        allowed = support_counts(3)
        forced = support_counts(3, mode='forced')
        self.assertEqual(sum(allowed.values()), 256**2)
        self.assertEqual(sum(forced.values()), 255**2)
        self.assertEqual(allowed[0], 1)
        self.assertNotIn(0, forced)
        self.assertEqual(forced_support_law(3), {3: F(1, 255), 65535: F(254, 255)})
        self.assertEqual(forced_support_law(1), {65534: F(1)})

    def test_quiet_pairs_cross_layer_exact_profiles(self):
        histogram = Counter()
        for a, b in combinations(range(16), 2):
            mask = (1 << a) | (1 << b)
            profile = one_layer_support_profile(mask)
            singles = profile.get(1, 0)+profile.get(2, 0)
            doubles = profile.get(3, 0)
            histogram[singles, doubles] += 1
        self.assertEqual(histogram, {(32, 0): 56, (28, 2): 64})
        self.assertEqual(one_layer_weight_law(3)[64], 1)
        mask = (1 << 0) | (1 << 8)
        law = one_layer_weight_law(mask)
        self.assertEqual(min(i for i, p in enumerate(law) if p), 58)
        self.assertEqual(law[58], F(2, 255)**2)
        columns = sum(1 << i for i in (0, 1, 2, 4, 6, 8, 10, 12, 14))
        law = one_layer_weight_law(mask, columns)
        self.assertEqual(min(i for i, p in enumerate(law) if p), 30)

    def test_profile_partition_and_complement_identity(self):
        classes = profile_classes((2,), ordered=False, max_masks=120)
        self.assertEqual(sum(classes.values()), 120)
        self.assertEqual(len(classes), 2)
        self.assertEqual(sorted(classes.values()), [56, 64])
        for mask in (0, 1, 3, 0x137, 0xaaaa, 0xff00, 0xffff):
            a = mask.bit_count()
            for delta in range(1, 16):
                self.assertEqual(delta_full_pairs(mask ^ 65535, delta), 8-a+delta_full_pairs(mask, delta))
            self.assertEqual(sum(n for _, n in pair_profile(mask, ordered=False)), 8)

    def test_interface_guards_and_cost(self):
        self.assertEqual(primitive_counts()['binary_xors'], 31)
        with self.assertRaises(ValueError):
            support_counts(1, mode='uniform')
        with self.assertRaises(ValueError):
            support_counts(1, size=3)
        with self.assertRaises(ValueError):
            profile_classes((2,), max_masks=119)


if __name__ == '__main__':
    unittest.main()

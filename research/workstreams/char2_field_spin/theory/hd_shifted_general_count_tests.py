import itertools
import unittest
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb

from hd_shifted_general_count import (
    ProfileBox, box_count_upper, exact_joint_count, exact_profiles,
    feasible_witness, joint_count_upper, joint_is_feasible, joint_run_count_upper,
    matching_histograms, oriented_pairing_count, transfer_value,
)
from hd_shifted_pairs_count import exact_profiles as weight_profiles


def exhaustive_profiles(pairs):
    result = Counter()
    for bits in itertools.product((0, 1), repeat=2 * pairs):
        even = sum(bits[::2])
        odd = sum(bits[1::2])
        a = sum(bits[i] * bits[i + 1] for i in range(0, 2 * pairs, 2))
        b = sum(bits[i] * bits[(i + 1) % (2 * pairs)]
                for i in range(1, 2 * pairs, 2))
        result[even, odd, a, b] += 1
    return result


class GeneralShiftedCountTests(unittest.TestCase):
    def test_small_dp_and_formula_against_exhaustive(self):
        for pairs in range(1, 7):
            expected = exhaustive_profiles(pairs)
            self.assertEqual(expected, exact_profiles(pairs))
            # Include impossible classes, not merely the DP's positive keys.
            for profile in itertools.product(range(pairs + 1), repeat=4):
                actual = exact_joint_count(pairs, *profile)
                self.assertEqual(expected[profile], actual, (pairs, profile))
                self.assertEqual(bool(actual), joint_is_feasible(pairs, *profile),
                                 (pairs, profile))
                self.assertLessEqual(actual, joint_count_upper(pairs, *profile))
                refined = joint_run_count_upper(pairs, *profile)
                self.assertLessEqual(actual, refined)
                self.assertLessEqual(refined, joint_count_upper(pairs, *profile))

    def test_witness_is_a_valid_transition_allocation(self):
        for pairs in range(1, 8):
            for even, odd, a, b in exact_profiles(pairs):
                runs, selected = feasible_witness(pairs, even, odd, a, b)
                capacities = (pairs - even - runs, runs, runs, even - runs)
                selections = (odd - a - b + selected, b - selected,
                              a - selected, selected)
                self.assertTrue(all(0 <= value <= capacity
                                    for value, capacity in zip(selections, capacities)))
                self.assertTrue(runs == 0 if even in (0, pairs)
                                else 1 <= runs <= min(even, pairs - even))

    def test_marginals_and_count_symmetries(self):
        for pairs in range(1, 9):
            profiles = exact_profiles(pairs)
            by_parity, by_matching, by_weight = (defaultdict(int) for _ in range(3))
            for (even, odd, a, b), count in profiles.items():
                self.assertEqual(count, profiles[odd, even, a, b])
                self.assertEqual(count, profiles[even, odd, b, a])
                self.assertEqual(count, profiles[pairs-even, pairs-odd,
                                                pairs-even-odd+a, pairs-even-odd+b])
                by_parity[even, odd] += count
                by_matching[even, odd, a] += count
                by_weight[even + odd, a, b] += count
                hist0, hist1 = matching_histograms(pairs, even, odd, a, b)
                self.assertEqual(sum(hist0), pairs)
                self.assertEqual(sum(hist1), pairs)
                self.assertGreaterEqual(min(hist0 + hist1), 0)
            for (even, odd), count in by_parity.items():
                self.assertEqual(count, comb(pairs, even) * comb(pairs, odd))
            for (even, odd, full), count in by_matching.items():
                self.assertEqual(count, oriented_pairing_count(pairs, even, odd, full))
            self.assertEqual(dict(by_weight), weight_profiles(pairs))
            self.assertEqual(sum(profiles.values()), 1 << (2 * pairs))

    def test_transfer_and_box_bounds(self):
        markers_list = ((1, 1, 1, 1), (2, Fraction(1, 3), 3, Fraction(1, 2)),
                        (Fraction(1, 4), 2, Fraction(2, 3), 5))
        for pairs in range(1, 6):
            profiles = exact_profiles(pairs)
            for markers in markers_list:
                expected = sum(count * Fraction(markers[0]) ** even
                               * Fraction(markers[1]) ** odd
                               * Fraction(markers[2]) ** a
                               * Fraction(markers[3]) ** b
                               for (even, odd, a, b), count in profiles.items())
                self.assertEqual(expected, transfer_value(pairs, markers))
            root = ProfileBox(*((0, pairs),) * 4)
            self.assertEqual(box_count_upper(pairs, root), 1 << (2 * pairs))
            boxes = [root]
            for _ in range(4):
                boxes = [child for box in boxes for child in box.split()]
            covered = Counter()
            for box in boxes:
                count = 0
                for profile, multiplicity in profiles.items():
                    if all(lo <= value <= hi for value, (lo, hi)
                           in zip(profile, box.ranges)):
                        covered[profile] += 1
                        count += multiplicity
                for markers in markers_list:
                    self.assertLessEqual(count, box_count_upper(pairs, box, markers))
            self.assertEqual(set(covered), set(profiles))
            self.assertEqual(set(covered.values()), {1})

    def test_boundaries_and_large_full_pair_reference(self):
        for pairs in (1, 2, 5, 4096):
            self.assertEqual(exact_joint_count(pairs, 0, 0, 0, 0), 1)
            self.assertEqual(exact_joint_count(pairs, pairs, pairs, pairs, pairs), 1)
            self.assertEqual(exact_joint_count(pairs, pairs, 0, 0, 0), 1)
            self.assertEqual(exact_joint_count(pairs, 0, pairs, 0, 0), 1)
            self.assertEqual(exact_joint_count(pairs, pairs, 1, 1, 1), pairs)
            self.assertEqual(exact_joint_count(pairs, 1, 0, 0, 0), pairs)
            self.assertEqual(joint_count_upper(pairs, -1, 0, 0, 0), 0)
        # Marginal histograms and total cyclic runs alone admit this false
        # class. The exact feasibility test rejects its mixed parity.
        self.assertFalse(joint_is_feasible(4, 2, 2, 0, 0))
        self.assertEqual(exact_joint_count(4096, 512, 512, 512, 511), 4096)
        count = exact_joint_count(4096, 512, 512, 512, 384)
        self.assertEqual(count, joint_count_upper(4096, 512, 512, 512, 384))

    def test_count_symmetry_does_not_swap_noncommuting_operators(self):
        # Initial row (1,0), terminal column (1,1), one alternating block.
        # Matching exchange changes the value despite preserving counts.
        left = ((1, 2), (0, 1))
        right = ((1, 0), (3, 1))
        forward = sum(left[0][k] * sum(right[k]) for k in range(2))
        reverse = sum(right[0][k] * sum(left[k]) for k in range(2))
        self.assertEqual((forward, reverse), (9, 3))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            exact_profiles(0)
        with self.assertRaises(ValueError):
            transfer_value(3, (1, 1, 0, 1))
        with self.assertRaises(TypeError):
            transfer_value(3, (1, 1, 0.5, 1))
        with self.assertRaises(ValueError):
            box_count_upper(3, ProfileBox((-1, 2), (0, 3), (0, 3), (0, 3)))
        with self.assertRaises(ValueError):
            ProfileBox((0, 0), (0, 0), (0, 0), (0, 0)).split()


if __name__ == "__main__":
    unittest.main()

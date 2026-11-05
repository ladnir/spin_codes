from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import unittest

from hd_outer_pair_layers import (
    BASELINE_DELTAS, VARIED_DELTAS, apply_pair, apply_component,
    pair_mask_law, component_mask_law, weight_law, evaluate_component,
    partner, components, closed_rowset, support_profile, primitive_counts,
    field_multiply, one_layer_partner, backward_support_kernel, backward_activity_kernel)
from hd_outer_pair_layers import one_layer_support_profile, one_layer_weight_law


def support(values):
    return sum(bool(value) << i for i, value in enumerate(values))


class PairLayerTests(unittest.TestCase):
    def test_exact_pair_law_for_every_fixed_GF4_input(self):
        for a, b in product(range(4), repeat=2):
            counts = Counter(support(apply_pair(a, b, alpha, modulus=0b111))
                             for alpha in range(1, 4))
            self.assertEqual(tuple(F(counts[i], 3) for i in range(4)), pair_mask_law(support((a, b)), 4))

    def test_two_stage_law_every_fixed_GF4_input(self):
        # Includes correlated representatives; no iid input-amplitude premise.
        for values in product(range(4), repeat=4):
            counts = Counter(support(apply_component(values, keys, modulus=0b111))
                             for keys in product(range(1, 4), repeat=4))
            self.assertEqual(tuple(F(counts[i], 81) for i in range(16)), component_mask_law(support(values), 4))

    def test_key_sharing_is_not_the_proved_ensemble(self):
        counts = Counter(support(apply_component((1, 1, 1, 1), (key,)*4, modulus=0b111))
                         for key in range(1, 4))
        self.assertNotEqual(tuple(F(counts[i], 3) for i in range(16)), component_mask_law(15, 4))

    def test_final_diagonal_gives_joint_uniform_nonzero_amplitudes(self):
        value = (1, 2, 0, 3)
        outputs = Counter(tuple(field_multiply(x, diagonal, modulus=0b111)
                                for x, diagonal in zip(value, keys))
                          for keys in product(range(1, 4), repeat=4))
        self.assertEqual(len(outputs), 3**3)
        self.assertEqual(set(outputs.values()), {3})
        self.assertTrue(all(support(y) == support(value) for y in outputs))

    def test_production_law_mass_and_sharp_local_floors(self):
        minimum = {1: 4, 2: 2, 3: 2, 4: 1}
        for mask in range(16):
            law = component_mask_law(mask)
            self.assertEqual(sum(law), 1)
            if mask:
                self.assertEqual(law[0], 0)
                self.assertEqual(min(out.bit_count() for out, probability in enumerate(law) if probability), minimum[mask.bit_count()])
        self.assertEqual(weight_law(15)[1], F(4, 255**3))
        self.assertEqual(weight_law(3)[2], F(2, 255))
        self.assertEqual(weight_law(5)[2], F(4, 255**2))

    def test_labelled_polynomial_not_only_weight(self):
        atoms = (F(1, 2), F(1, 3), F(1, 5), F(1, 7))
        for mask in range(16):
            expected = F(0)
            for output, probability in enumerate(component_mask_law(mask, 4)):
                monomial = F(1)
                for i, atom in enumerate(atoms):
                    if output >> i & 1:
                        monomial *= atom
                expected += probability*monomial
            self.assertEqual(evaluate_component(mask, atoms, 4), expected)
        self.assertEqual(evaluate_component(1, atoms), atoms[0]*atoms[1]*atoms[2]*atoms[3])
        self.assertNotEqual(component_mask_law(3), component_mask_law(5))

    def test_backward_kernel_retains_original_activity(self):
        atoms = (F(1, 2), F(1, 3), F(1, 5), F(1, 7))
        output_weights = []
        for mask in range(16):
            term = F(1)
            for i, atom in enumerate(atoms):
                if mask >> i & 1:
                    term *= atom
            output_weights.append(term)
        # Sum all GF4 input arrays; each fixed support has3^weight values.
        inputs = tuple(3**mask.bit_count() for mask in range(16))
        labels = (2, 0, 3, 1)
        got = backward_activity_kernel(inputs, output_weights, q=4, column_labels=labels)
        expected = [F(0) for _ in range(16)]
        for values in product(range(4), repeat=4):
            flag = sum(bool(value) << labels[i] for i, value in enumerate(values))
            for keys in product(range(1, 4), repeat=4):
                output = support(apply_component(values, keys, modulus=0b111))
                expected[flag] += output_weights[output]/81
        self.assertEqual(got, tuple(expected))
        self.assertEqual(backward_support_kernel((1,)*16), (1,)*16)
        # A nonzero input flag is never inferred from the final support.
        self.assertTrue(component_mask_law(15, 4)[1] > 0)

    def test_all_superblock_coordinates_exactly_once(self):
        for schedule in (BASELINE_DELTAS, VARIED_DELTAS):
            pieces = components(schedule)
            self.assertEqual(len(pieces), 64)
            self.assertEqual(sorted(p for part in pieces for p in part), list(range(256)))
            for point in range(256):
                for stage in (0, 1):
                    self.assertEqual(partner(partner(point, stage, schedule), stage, schedule), point)
            self.assertTrue(all(len({(p % 16)//4 for p in part}) == 1 for part in pieces))
        for point in range(256):
            self.assertEqual(one_layer_partner(one_layer_partner(point)), point)

    def test_preferred_schedule_subset_spans(self):
        def span(generators):
            values = {0}
            for generator in generators:
                values |= {v ^ generator for v in tuple(values)}
            return values
        for selected in combinations(range(8, 16), 5):
            self.assertEqual(span(selected), set(range(16)))
        for selected in combinations(VARIED_DELTAS, 3):
            self.assertEqual(span((1,)+selected), set(range(16)))

    def test_all140_affine_planes_one_layer_profiles(self):
        histogram = Counter()
        for rows in combinations(range(16), 4):
            if rows[0] ^ rows[1] ^ rows[2] ^ rows[3]:
                continue
            mask = sum(1 << row for row in rows)
            profile = one_layer_support_profile(mask)
            singles, doubles = profile.get(1, 0)+profile.get(2, 0), profile.get(3, 0)
            histogram[singles, doubles] += 1
        self.assertEqual(histogram, {(64, 0): 28, (48, 8): 112})
        law = one_layer_weight_law(15)
        self.assertEqual(law[128], 1)
        # The aligned plane{0,1,8,9} attains the104 full-column floor.
        mask = sum(1 << i for i in (0, 1, 8, 9))
        law = one_layer_weight_law(mask)
        self.assertEqual(min(i for i, p in enumerate(law) if p), 104)
        self.assertEqual(law[104], F(2, 255)**8)
        # Two aligned fully active pairs plus five singly active pairs.
        columns = sum(1 << i for i in (0, 1, 2, 3, 4, 6, 8, 10, 12))
        law = one_layer_weight_law(mask, columns)
        self.assertEqual(min(i for i, p in enumerate(law) if p), 48)

    def test_baseline_quiet_plane_counterrisk_and_varied_wiring(self):
        from hd_bit_syndrome_outer import direct
        quiet = tuple(int(i < 4) for i in range(16))
        self.assertEqual(direct(quiet), quiet)
        self.assertTrue(closed_rowset(15, BASELINE_DELTAS))
        self.assertEqual(support_profile(15, BASELINE_DELTAS), {0: 48, 15: 16})
        self.assertGreater(evaluate_component(15, (F(1, 2),)*4), F(1, 16))
        self.assertFalse(closed_rowset(15, VARIED_DELTAS))
        self.assertEqual(support_profile(15, VARIED_DELTAS), {0: 36, 3: 24, 15: 4})
        # Closure under the varied generators reaches every row from0.
        orbit = {0}
        for generator in (1, 2, 4, 8):
            orbit |= {v ^ generator for v in tuple(orbit)}
        self.assertEqual(orbit, set(range(16)))
        self.assertTrue(closed_rowset(0, VARIED_DELTAS))
        self.assertTrue(closed_rowset(65535, VARIED_DELTAS))

    def test_counts_and_invalid_parameters(self):
        counts = primitive_counts()
        self.assertEqual((counts['variable_products'], counts['binary_xors'], counts['fresh_nonzero_alpha_keys']), (512, 512, 256))
        other = primitive_counts(precompute_c_alpha=False)
        self.assertEqual((other['variable_products'], other['fixed_products']), (256, 256))
        with self.assertRaises(ValueError):
            pair_mask_law(3, 2)
        with self.assertRaises(ValueError):
            apply_pair(1, 1, 0)
        with self.assertRaises(ValueError):
            components((1, 1, 1, 1))


if __name__ == '__main__':
    unittest.main()

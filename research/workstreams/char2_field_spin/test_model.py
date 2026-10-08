"""Exact validation of the baseline model; no distance or timing claims."""

from fractions import Fraction
from collections import Counter
import itertools
import random
import unittest

try:
    from . import model
except ImportError:
    import model


def enumerated_transition(field, input_value):
    """Count every state/coefficient assignment, compressing equal histograms.

    The coefficient-product histograms and XOR convolution count all triples
    exactly. Caching identical histograms avoids a cubic GF256 multiplication
    loop; it does not assume a uniform next-state distribution.
    """
    nonzero = range(1, field.order)

    def product_histogram(value):
        counts = [0] * field.order
        for coefficient in nonzero:
            counts[field.mul(coefficient, value)] += 1
        return tuple(counts)

    input_histogram = product_histogram(input_value)
    cache = {}
    result = []
    for state_values in ([0], nonzero):
        counts = [[0, 0], [0, 0]]
        denominator = len(state_values) * (field.order - 1) ** 2
        for state in state_values:
            state_histogram = product_histogram(state)
            if state_histogram not in cache:
                convolution = [0] * field.order
                for left, left_count in enumerate(state_histogram):
                    for right, right_count in enumerate(input_histogram):
                        convolution[left ^ right] += left_count * right_count
                cache[state_histogram] = convolution
            distribution = cache[state_histogram]
            emission = int(input_value ^ state != 0)
            counts[0][emission] += distribution[0]
            counts[1][emission] += sum(distribution[1:])
            # Conditional on a nonzero successor, it must be uniform as well.
            if any(distribution[1:]):
                assert len(set(distribution[1:])) == 1
        result.append(tuple(tuple(Fraction(c, denominator) for c in pair)
                            for pair in counts))
    return tuple(result)


def directly_enumerated_transition(field, input_value):
    result = []
    for states in ([0], range(1, field.order)):
        counts = [[0, 0], [0, 0]]
        for state, a, b in itertools.product(states, range(1, field.order),
                                            range(1, field.order)):
            output = input_value ^ state
            successor = field.mul(a, state) ^ field.mul(b, input_value)
            counts[int(successor != 0)][int(output != 0)] += 1
        denominator = len(states) * (field.order - 1) ** 2
        result.append(tuple(tuple(Fraction(c, denominator) for c in pair)
                            for pair in counts))
    return tuple(result)


class FieldTests(unittest.TestCase):
    def test_gf4_table(self):
        expected = [[0, 0, 0, 0], [0, 1, 2, 3],
                    [0, 2, 3, 1], [0, 3, 1, 2]]
        for left, right in itertools.product(range(4), repeat=2):
            self.assertEqual(model.GF4.mul(left, right), expected[left][right])

    def test_gf256_known_product_and_inverses(self):
        self.assertEqual(model.GF256.mul(0x57, 0x83), 0xC1)
        for value in range(1, 256):
            self.assertEqual(model.GF256.mul(value, model.GF256.inverse(value)), 1)
        with self.assertRaises(ZeroDivisionError):
            model.GF256.inverse(0)

    def test_distributivity(self):
        for left, a, b in itertools.product(range(4), repeat=3):
            self.assertEqual(model.GF4.mul(left, a ^ b),
                             model.GF4.mul(left, a) ^ model.GF4.mul(left, b))


class TransitionTests(unittest.TestCase):
    def test_gf4_all_inputs_direct_and_histogram(self):
        for value in range(4):
            expected = model.marked_support_transition(4, value != 0)
            self.assertEqual(directly_enumerated_transition(model.GF4, value), expected)
            self.assertEqual(enumerated_transition(model.GF4, value), expected)

    def test_gf256_inputs_and_conditional_uniformity(self):
        for value in (0, 1, 0x57, 0xFF):
            self.assertEqual(enumerated_transition(model.GF256, value),
                             model.marked_support_transition(256, value != 0))

    def test_rows_are_probability_distributions(self):
        for order, active in itertools.product((4, 256), (False, True)):
            for row in model.marked_support_transition(order, active):
                self.assertEqual(sum(sum(polynomial) for polynomial in row), 1)


class QuotientRankTests(unittest.TestCase):
    def test_gl2_row_law_over_gf4(self):
        """Test row forms after a uniform GL(2) map, not just one message."""
        field = model.GF4
        matrices = [entries for entries in itertools.product(range(4), repeat=4)
                    if field.mul(entries[0], entries[3]) ^
                    field.mul(entries[1], entries[2])]
        self.assertEqual(len(matrices), 180)
        fixed_row = (1, 2)
        rows = []
        for m00, m01, m10, m11 in matrices:
            rows.append((field.mul(fixed_row[0], m00) ^ field.mul(fixed_row[1], m10),
                         field.mul(fixed_row[0], m01) ^ field.mul(fixed_row[1], m11)))
        row_counts = Counter(rows)
        self.assertEqual(len(row_counts), 15)
        self.assertNotIn((0, 0), row_counts)
        self.assertEqual(set(row_counts.values()), {12})

        for rank in (1, 2):
            # Right multiplication by Q=diag(1,0), or Q=I, respectively.
            images = Counter((left, right if rank == 2 else 0) for left, right in rows)
            zero_probability = Fraction(4 ** (2 - rank) - 1, 4 ** 2 - 1)
            nonzero_probability = Fraction(4 ** (2 - rank), 4 ** 2 - 1)
            self.assertEqual(Fraction(images[(0, 0)], len(matrices)), zero_probability)
            expected_space = ((left, right)
                              for left in range(4)
                              for right in (range(4) if rank == 2 else (0,)))
            for image in expected_space:
                if image != (0, 0):
                    self.assertEqual(Fraction(images[image], len(matrices)),
                                     nonzero_probability)
            self.assertEqual(sum(images.values()), 180)


class EncoderTests(unittest.TestCase):
    def test_exhaustive_gf4_length_two(self):
        field = model.GF4
        vectors = list(itertools.product(range(4), repeat=2))
        for coefficients in itertools.product(range(1, 4), repeat=4):
            a, b = coefficients[:2], coefficients[2:]
            delta, gamma = model.independent_scales(field, a, b)
            images = set()
            for values in vectors:
                output = model.recursive_encode(field, values, a, b)
                images.add(tuple(output))
                self.assertEqual(model.recursive_decode(field, output, a, b), list(values))
                self.assertEqual(model.scaled_prefix_encode(field, values, delta, gamma), output)
                self.assertEqual(model.scaled_suffix_transpose(field, values, delta, gamma),
                                 model.recursive_transpose(field, values, a, b))
            self.assertEqual(len(images), 16)

    def test_random_lengths_adjoint_and_inverse(self):
        rng = random.Random(20261003)
        for field in (model.GF4, model.GF256):
            for length in (0, 1, 2, 3, 8, 17, 64):
                for _ in range(8):
                    x = [rng.randrange(field.order) for _ in range(length)]
                    z = [rng.randrange(field.order) for _ in range(length)]
                    a = [rng.randrange(1, field.order) for _ in range(length)]
                    b = [rng.randrange(1, field.order) for _ in range(length)]
                    delta, gamma = model.independent_scales(field, a, b)
                    output = model.recursive_encode(field, x, a, b)
                    transpose = model.scaled_suffix_transpose(field, z, delta, gamma)
                    self.assertEqual(output, model.scaled_prefix_encode(field, x, delta, gamma))
                    self.assertEqual(transpose, model.recursive_transpose(field, z, a, b))
                    self.assertEqual(model.recursive_decode(field, output, a, b), x)
                    self.assertEqual(model.field_dot(field, output, z),
                                     model.field_dot(field, x, transpose))

    def test_independent_scale_distribution(self):
        # Retain delta_n here: the full map of nonzero coefficients is bijective.
        field, length = model.GF4, 2
        images = set()
        for values in itertools.product(range(1, 4), repeat=2 * length):
            a, b = values[:length], values[length:]
            delta, gamma = model.independent_scales(field, a, b)
            delta_final = field.mul(delta[-1], a[-1])
            images.add(tuple(delta[1:] + [delta_final] + gamma))
        self.assertEqual(len(images), 3 ** (2 * length))

    def test_rejects_zero_coefficients_and_mismatched_lengths(self):
        with self.assertRaises(ValueError):
            model.recursive_encode(model.GF4, [1], [0], [1])
        with self.assertRaises(ValueError):
            model.scaled_prefix_encode(model.GF4, [1], [1], [])


class GVCalibrationTests(unittest.TestCase):
    def test_asymptotic_half_rate_reference_points(self):
        expected = {8: 0.3804717309726672, 32: 0.46883761881173835,
                    64: 0.48438599319413844, 128: 0.4921888754327439}
        for m, distance in expected.items():
            self.assertAlmostEqual(model.gv_distance(0.5, m), distance, places=13)
            self.assertAlmostEqual(model.qary_entropy(distance, m), 0.5, places=14)

    def test_asymptotic_ten_percent_rate_reference_points(self):
        expected = {8: 0.8414461328406041, 32: 0.8853438877013847,
                    64: 0.8926719438501675, 128: 0.8963359719250837}
        for m, rate in expected.items():
            self.assertAlmostEqual(model.gv_rate(0.1, m), rate, places=13)
            self.assertAlmostEqual(model.gv_distance(rate, m), 0.1, places=13)

    def test_endpoints_and_invalid_domains(self):
        for m in (8, 32, 64, 128):
            self.assertEqual(model.qary_entropy(0, m), 0)
            self.assertEqual(model.gv_distance(1, m), 0)
            self.assertEqual(model.gv_distance(0, m), 1 - 2 ** (-m))
        with self.assertRaises(ValueError):
            model.gv_distance(-0.1, 8)
        with self.assertRaises(ValueError):
            model.qary_entropy(0.1, 0)


if __name__ == "__main__":
    unittest.main()

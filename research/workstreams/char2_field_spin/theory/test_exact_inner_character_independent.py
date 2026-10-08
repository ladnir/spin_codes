"""Independent GF4 enumeration of the complete marked inner transfer."""

from collections import Counter
from fractions import Fraction
from itertools import product
import unittest

from exact_inner_character import build_histogram, moment


def gf4_product(a, b):
    # Independent two-bit polynomial product and reduction by X^2+X+1.
    polynomial = (a if b & 1 else 0) ^ ((a << 1) if b & 2 else 0)
    return polynomial ^ (7 if polynomial & 4 else 0)


def direct_counts(state, width):
    columns = []
    for point in range(width):
        powers = [1]
        for _ in range(1, state):
            powers.append(gf4_product(powers[-1], point))
        columns.append(powers)
    inputs = []
    for x in product(range(4), repeat=width):
        syndrome = [0] * state
        for coordinate in range(width):
            for power in range(state):
                syndrome[power] ^= gf4_product(columns[coordinate][power], x[coordinate])
        inputs.append((x, tuple(syndrome), sum(value != 0 for value in x)))
    counts = [[Counter(), Counter()], [Counter(), Counter()]]
    for refreshed in product(range(4), repeat=state):
        nonzero_state = int(any(refreshed))
        emitted = []
        for column in columns:
            value = 0
            for coefficient, state_value in zip(column, refreshed):
                value ^= gf4_product(coefficient, state_value)
            emitted.append(value)
        for x, syndrome, input_weight in inputs:
            return_is_nonzero = int(syndrome != refreshed)
            output_weight = sum(a != b for a, b in zip(x, emitted))
            counts[nonzero_state][return_is_nonzero][input_weight, output_weight] += 1
    return counts


def evaluate_counts(counts, state, v, z):
    return [
        [sum(count * v**input_weight * z**output_weight
             for (input_weight, output_weight), count in entry.items())
         / (4**state - 1 if row else 1)
         for entry in entries]
        for row, entries in enumerate(counts)
    ]


class IndependentCharacterTransferTest(unittest.TestCase):
    def test_all_states_inputs_and_signed_factor_regimes(self):
        # Width-three cases also test maps with no full-field moment sums.
        for state, width in ((1, 4), (2, 4), (3, 4), (2, 3), (3, 3)):
            counts = direct_counts(state, width)
            histogram = build_histogram(4, state, width)
            for v, z in ((Fraction(1, 7), Fraction(2, 5)),
                         (Fraction(2, 3), Fraction(1, 5)),
                         (Fraction(0), Fraction(1, 3)),
                         (Fraction(3, 2), Fraction(1))):
                with self.subTest(state=state, width=width, v=v, z=z):
                    self.assertEqual(moment(v, z, histogram),
                                     evaluate_counts(counts, state, v, z))


if __name__ == "__main__":
    unittest.main()

"""Independent exact checks for the GF(q) pre-output-refresh alternative.

These checks concern uniform GL(2,q) refresh and independently uniform
nonzero input scalars, not the native-field small-set three-shear variant.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arithmetic import multiply


def masks(labels):
    for h in range(len(labels) + 1):
        yield from combinations(labels, h)


def syndrome(labels, values, modulus):
    low = high = 0
    for j, value in zip(labels, values):
        low ^= value
        high ^= multiply(j, value, modulus)
    return low, high


def nonzero_solutions(labels, target, q, modulus):
    """Closed count for distinct columns (1,j), independently audited."""
    n = len(labels)
    if target == (0, 0):
        correction = q * q - 1 - n * q * (q - 1)
    else:
        low, high = target
        on_line = int(low != 0 and any(
            multiply(j, low, modulus) == high for j in labels))
        correction = -1 + n * q - on_line * q * q
    numerator = (q - 1) ** n + (-1) ** n * correction
    assert numerator % (q * q) == 0
    return numerator // (q * q)


def enumerated_law(labels, active, q, modulus, incoming_active):
    """Enumerate every R and every active X, independently of the formula."""
    states = [(x, y) for x, y in product(range(q), repeat=2)
              if (x, y) != (0, 0)] if incoming_active else [(0, 0)]
    count = Counter()
    for r in states:
        for values in product(range(1, q), repeat=len(active)):
            word = dict(zip(active, values))
            weight = sum((word.get(j, 0) ^ r[0]
                          ^ multiply(j, r[1], modulus)) != 0 for j in labels)
            sx = syndrome(active, values, modulus)
            count[weight, sx == r] += 1
    denominator = len(states) * (q - 1) ** len(active)
    return {key: F(value, denominator) for key, value in count.items()}


def expanded_return_polynomial(labels, active, q, modulus):
    """Subset expansion, retaining actual lane labels rather than only h."""
    polynomial = [0] * (len(labels) + 1)
    for r0, r1 in product(range(q), repeat=2):
        if (r0, r1) == (0, 0):
            continue
        shift = {j: r0 ^ multiply(j, r1, modulus) for j in labels}
        inactive_weight = sum(shift[j] != 0 for j in labels if j not in active)
        for selected in masks(active):
            if any(shift[j] == 0 for j in selected):
                continue
            fixed = syndrome(selected, [shift[j] for j in selected], modulus)
            target = r0 ^ fixed[0], r1 ^ fixed[1]
            remaining = tuple(j for j in active if j not in selected)
            ways = nonzero_solutions(remaining, target, q, modulus)
            power = inactive_weight + len(active) - len(selected)
            for degree in range(len(selected) + 1):
                polynomial[power + degree] += (
                    ways * comb(len(selected), degree) * (-1) ** degree)
    denominator = (q * q - 1) * (q - 1) ** len(active)
    return tuple(F(value, denominator) for value in polynomial)


class JointTransferAudit(unittest.TestCase):
    def test_nonzero_solution_count(self):
        for q, modulus, max_support in ((4, 0x7, 4), (8, 0xB, 3)):
            for h in range(max_support + 1):
                for labels in combinations(range(q), h):
                    counts = Counter(syndrome(labels, values, modulus)
                                     for values in product(range(1, q), repeat=h))
                    for target in product(range(q), repeat=2):
                        self.assertEqual(nonzero_solutions(labels, target, q, modulus),
                                         counts[target])

    def test_subset_expansion_every_gf4_support(self):
        labels = tuple(range(4))
        for active in masks(labels):
            law = enumerated_law(labels, active, 4, 0x7, True)
            returned = tuple(law.get((weight, True), F(0)) for weight in range(5))
            self.assertEqual(expanded_return_polynomial(labels, active, 4, 0x7),
                             returned)
            self.assertEqual(sum(law.values()), 1)

    def test_inactive_transition(self):
        labels = tuple(range(4))
        for active in masks(labels):
            law = enumerated_law(labels, active, 4, 0x7, False)
            h = len(active)
            returning = F(nonzero_solutions(active, (0, 0), 4, 0x7), 3 ** h)
            self.assertEqual(law.get((h, True), 0), returning)
            self.assertEqual(law.get((h, False), 0), 1 - returning)
            self.assertTrue(all(weight == h for weight, _ in law))

    def test_joint_geometry_is_not_weight_only(self):
        labels = tuple(range(4))
        first = enumerated_law(labels, (0,), 4, 0x7, True)
        second = enumerated_law(labels, (1,), 4, 0x7, True)
        self.assertEqual(first.get((3, True)), F(1, 15))
        self.assertEqual(second.get((4, True)), F(1, 15))
        self.assertNotEqual(first, second)
        # The unconditional output marginal nevertheless agrees.
        for weight in range(5):
            self.assertEqual(sum(first.get((weight, returned), 0)
                                 for returned in (False, True)),
                             sum(second.get((weight, returned), 0)
                                 for returned in (False, True)))

    def test_uniform_sl2_is_enough_for_fixed_nonzero_state(self):
        q, modulus = 4, 0x7
        matrices = [(a, b, c, d) for a, b, c, d in product(range(q), repeat=4)
                    if multiply(a, d, modulus) ^ multiply(b, c, modulus) == 1]
        self.assertEqual(len(matrices), q * (q * q - 1))
        for x, y in product(range(q), repeat=2):
            if (x, y) == (0, 0):
                continue
            images = Counter((multiply(a, x, modulus) ^ multiply(b, y, modulus),
                              multiply(c, x, modulus) ^ multiply(d, y, modulus))
                             for a, b, c, d in matrices)
            self.assertEqual(len(images), q * q - 1)
            self.assertEqual(set(images.values()), {q})

    def test_subset_check_matrix_depends_on_size_parity_and_xor(self):
        modulus = 0x13
        for selected in masks(tuple(range(8))):
            matrix = [0, 0, 0, 0]
            label_xor = 0
            for j in selected:
                label_xor ^= j
                matrix[0] ^= 1
                matrix[1] ^= j
                matrix[2] ^= j
                matrix[3] ^= multiply(j, j, modulus)
            self.assertEqual(matrix, [len(selected) % 2, label_xor, label_xor,
                                      multiply(label_xor, label_xor, modulus)])


if __name__ == "__main__":
    unittest.main()

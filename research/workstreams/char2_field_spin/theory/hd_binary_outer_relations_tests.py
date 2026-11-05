from collections import Counter
from itertools import combinations, product
import random
import unittest

from hd_binary_outer_relations import (
    kernel_subspaces, injection_count, joint_counts,
    joint_counts_from_matrix, select_columns)
from hd_binary_syndrome_outer import matrix, two_row_joint_counts


def brute(forms, arity, q):
    result = Counter()
    for values in product(range(q), repeat=arity):
        mask = 0
        for j, form in enumerate(forms):
            value = 0
            for i, x in enumerate(values):
                if form >> i & 1:
                    value ^= x
            mask |= bool(value) << j
        result[sum(x != 0 for x in values), mask] += 1
    return result


class BinaryRelationsTests(unittest.TestCase):
    def test_complete_subspace_counts_and_closure(self):
        self.assertEqual([len(kernel_subspaces(a)) for a in range(5)], [1, 2, 5, 16, 67])
        for arity in range(5):
            for kernel in kernel_subspaces(arity):
                members = [v for v in range(1 << arity) if kernel >> v & 1]
                self.assertIn(0, members)
                self.assertEqual(len(members) & (len(members)-1), 0)
                self.assertTrue(all(kernel >> (u ^ v) & 1 for u in members for v in members))

    def test_injection_multiplicity_and_small_fields(self):
        self.assertEqual([injection_count(d, 4) for d in range(5)], [1, 3, 6, 0, 0])
        for arity in range(5):
            for q in (2, 4, 256):
                identity = tuple(1 << i for i in range(arity))
                counts = joint_counts(identity, arity, q=q)
                self.assertEqual(sum(counts.values()), q**arity)
                self.assertEqual(counts[0, 0], 1)

    def test_GF4_all_inputs_for_selected_a1_to_a4(self):
        rng = random.Random(20261009)
        for arity in range(5):
            fixtures = [tuple(range(1 << arity)), (0,)*5]
            fixtures += [tuple(rng.randrange(1 << arity) for _ in range(16)) for _ in range(4)]
            for forms in fixtures:
                self.assertEqual(joint_counts(forms, arity, q=4), brute(forms, arity, 4))

    def test_all_two_row_interfaces_and_source_order(self):
        for mixer in ('parity', 'binary_syndrome'):
            for a, b in combinations(range(16), 2):
                restricted = select_columns(matrix(mixer), (a, b))
                self.assertEqual(joint_counts_from_matrix(restricted),
                                 two_row_joint_counts(a, b, mixer=mixer))
                reversed_sources = select_columns(matrix(mixer), (b, a))
                self.assertEqual(joint_counts_from_matrix(restricted),
                                 joint_counts_from_matrix(reversed_sources))

    def test_arbitrary_singular_matrix_and_labelled_outputs(self):
        forms = (0, 1, 1, 3, 7, 8, 15)
        counts = joint_counts(forms, 4, q=4)
        self.assertEqual(counts, brute(forms, 4, 4))
        reversed_counts = joint_counts(tuple(reversed(forms)), 4, q=4)
        expected = Counter()
        for (weight, mask), count in counts.items():
            remapped = sum(((mask >> i) & 1) << (len(forms)-1-i) for i in range(len(forms)))
            expected[weight, remapped] += count
        self.assertEqual(reversed_counts, expected)
        self.assertEqual(joint_counts((), 4), {(w, 0): n for (w, _), n in
                          sorted(joint_counts((0,), 4).items())})

    def test_bounds_and_input_validation(self):
        for args in (((0,), 5), ((2,), 1)):
            with self.assertRaises(ValueError):
                joint_counts(*args)
        with self.assertRaises(ValueError):
            joint_counts((1,), 1, q=3)
        with self.assertRaises(ValueError):
            joint_counts_from_matrix(((1, 0), (1,)))
        with self.assertRaises(ValueError):
            select_columns(matrix(), (2, 2))


if __name__ == '__main__':
    unittest.main()

from collections import Counter
from itertools import combinations, product
import random
import unittest

from hd_quaternary_outer_relations import (
    AES_EMBEDDING, multiply_gf4, embed_aes, pack_row,
    kernel_subspaces, kernel_bitsets, injection_count, category_data,
    column_masks, column_categories, joint_counts,
    row_forms_from_matrix, joint_counts_from_matrix)
from hd_lambda_outer import (
    matrix as lambda_matrix, joint_counts as lambda_counts,
    column_branches as lambda_branches)
from hd_outer_pair_layers import field_multiply


def brute(forms, arity):
    result = Counter()
    for values in product(range(4), repeat=arity):
        mask = 0
        for j, row in enumerate(forms):
            value = 0
            for i, x in enumerate(values):
                value ^= multiply_gf4(row >> (2*i) & 3, x)
            mask |= bool(value) << j
        result[sum(x != 0 for x in values), mask] += 1
    return result


class QuaternaryRelationsTests(unittest.TestCase):
    def test_embedding_is_explicit_field_homomorphism(self):
        self.assertEqual(AES_EMBEDDING, (0, 1, 0xbc, 0xbd))
        for a, b in product(range(4), repeat=2):
            self.assertEqual(embed_aes(a ^ b), embed_aes(a) ^ embed_aes(b))
            self.assertEqual(embed_aes(multiply_gf4(a, b)),
                             field_multiply(embed_aes(a), embed_aes(b)))

    def test_canonical_subspaces_counts_and_closure(self):
        self.assertEqual([len(kernel_subspaces(a)) for a in range(5)], [1, 2, 7, 44, 529])
        for arity in range(5):
            self.assertEqual(len(set(kernel_bitsets(arity))), len(kernel_subspaces(arity)))
            for rows, kernel in zip(kernel_subspaces(arity), kernel_bitsets(arity)):
                self.assertEqual(kernel.bit_count(), 4**len(rows))
                pivots = []
                for row in rows:
                    pivot = next(i for i in range(arity) if row >> (2*i) & 3)
                    pivots.append(pivot)
                    self.assertEqual(row >> (2*pivot) & 3, 1)
                self.assertEqual(pivots, sorted(pivots))
                self.assertTrue(all((row >> (2*pivot) & 3) == int(i == j)
                                    for i, row in enumerate(rows) for j, pivot in enumerate(pivots)))
                members = [v for v in range(4**arity) if kernel >> v & 1]
                for v in members:
                    for row in rows:
                        self.assertTrue(kernel >> (v ^ row) & 1)
                    scaled = sum(multiply_gf4(v >> (2*i) & 3, 2) << (2*i) for i in range(arity))
                    self.assertTrue(kernel >> scaled & 1)

    def test_counts_and_stable_input_categories(self):
        self.assertEqual([injection_count(d, 4) for d in range(5)], [1, 3, 0, 0, 0])
        for arity in range(5):
            for q in (4, 16, 256):
                data = category_data(arity, q)
                self.assertEqual(sum(data['category_counts']), q**arity)
                by_weight = Counter()
                by_mask = Counter()
                for mask, weight, count in zip(data['category_input_masks'],
                                               data['category_input_weights'], data['category_counts']):
                    self.assertEqual(weight, mask.bit_count())
                    by_weight[weight] += count
                    by_mask[mask] += count
                from math import comb
                self.assertEqual(dict(by_weight), {w: comb(arity, w)*(q-1)**w for w in range(arity+1)})
                self.assertEqual(dict(by_mask), {mask: (q-1)**mask.bit_count() for mask in range(1 << arity)})
            self.assertEqual(category_data(arity, 4)['kernel_rrefs'], category_data(arity, 256)['kernel_rrefs'])

    def test_GF4_exhaustive_allowed_inputs_through_a4(self):
        rng = random.Random(4316)
        for arity in range(5):
            fixtures = [tuple(range(4**arity)), (0,)*4]
            fixtures += [tuple(rng.randrange(4**arity) for _ in range(16)) for _ in range(3)]
            for forms in fixtures:
                self.assertEqual(joint_counts(forms, arity, q=4), brute(forms, arity))
                categories = column_categories(forms, arity, q=4)
                self.assertEqual(len(categories), len(kernel_subspaces(arity)))
                self.assertEqual(tuple(w for w, _, _ in categories), category_data(arity, 4)['category_input_weights'])

    def test_all120_lambda_pairs_and_parameter_support_transfer(self):
        for source in tuple((i,) for i in range(16))+tuple(combinations(range(16), 2)):
            for lam in (0xbc, 0xbd):
                selected = tuple(tuple(row[i] for i in source) for row in lambda_matrix(lam))
                result = joint_counts_from_matrix(selected, q=256, encoding='aes')
                self.assertEqual(result, lambda_counts(source, lambda_value=lam))
                self.assertEqual(result, lambda_counts(source, lambda_value=2))
                self.assertEqual(lambda_branches(source, lambda_value=lam), lambda_branches(source, lambda_value=2))
            toy = tuple(tuple(row[i] for i in source) for row in lambda_matrix(2, 4))
            self.assertEqual(joint_counts_from_matrix(toy, q=4), lambda_counts(source, q=4))

    def test_ordered_labels_and_nonbinary_coefficients(self):
        rows = ((1, 2, 0), (3, 1, 1), (0, 2, 3), (1, 0, 0))
        forms = row_forms_from_matrix(rows)
        self.assertEqual(forms, tuple(pack_row(row) for row in rows))
        masks = column_masks(forms, 3)
        reversed_masks = column_masks(tuple(reversed(forms)), 3)
        expected = tuple(sum(((mask >> i) & 1) << (len(rows)-1-i) for i in range(len(rows))) for mask in masks)
        self.assertEqual(reversed_masks, expected)
        self.assertEqual(joint_counts(forms, 3, q=4), brute(forms, 3))
        aes = tuple(tuple(embed_aes(value) for value in row) for row in rows)
        self.assertEqual(row_forms_from_matrix(aes, encoding='aes'), forms)

    def test_explicit_boundaries_and_reject_wrong_AES_parameter(self):
        with self.assertRaises(ValueError):
            kernel_subspaces(5)
        with self.assertRaises(ValueError):
            joint_counts((4,), 1)
        with self.assertRaises(ValueError):
            category_data(2, 8)
        with self.assertRaises(ValueError):
            row_forms_from_matrix(((2,),), encoding='aes')
        with self.assertRaises(ValueError):
            row_forms_from_matrix(((1, 2), (1,)))
        with self.assertRaises(ValueError):
            multiply_gf4(4, 0)


if __name__ == '__main__':
    unittest.main()

"""Tiny exact circuit checks only; no moment, provider, or native runs."""
import random
import unittest

import hd_lambda_double_transvection as target


class SpecializedCircuitTests(unittest.TestCase):
    def test_all_binary_basis_vectors_in_both_fields(self):
        for q, lam in ((4, 2), (256, 188)):
            for i in range(16):
                x = tuple(int(i == j) for j in range(16))
                self.assertEqual(target.circuit(x, lam, q), target.direct(x, lam, q))
                self.assertEqual(target.circuit(target.circuit(x, lam, q), lam, q), x)

    def test_all_GF4_scalars_on_basis_and_deterministic_random_cases(self):
        rng = random.Random(163)
        for i in range(16):
            for scalar in range(4):
                x = tuple(scalar if i == j else 0 for j in range(16))
                self.assertEqual(target.circuit(x, 2, 4), target.direct(x, 2, 4))
        for q, lam in ((4, 2), (256, 188)):
            for _ in range(32):
                x = tuple(rng.randrange(q) for _ in range(16))
                actual = target.circuit(x, lam, q)
                self.assertEqual(actual, target.direct(x, lam, q))
                self.assertEqual(actual, target.core.circuit(x, target.SUPPORTS, lam, q))

    def test_count_is64_and_fixed_product1(self):
        counts = target.primitive_counts()
        self.assertEqual((counts['binary_xors'], counts['fixed_products'], counts['variable_products']), (64, 1, 0))
        self.assertEqual(counts['extra_xors_over_Hlambda'], 7)
        self.assertEqual(target.core.primitive_counts(target.SUPPORTS)['extra_xors_over_Hlambda'], 9)
        self.assertFalse(counts['native_performance_claim'])
        self.assertFalse(counts['old_receipts_transferred'])

    def test_full_matrix_identity_and_rank(self):
        h = target.matrix()
        self.assertEqual(h, tuple(zip(*h)))
        self.assertEqual(set(x for row in h for x in row), {188, 189})
        self.assertEqual(target.core.syndrome_extension(target.SUPPORTS)['correction_rank_if_independent'], 7)
        self.assertEqual(target.CONSTRUCTION, 'g16-hlambda-bc-commuting-quiet-transvections-000f-00c3-research')
        self.assertEqual(target.core.pair_profile(target.SUPPORTS)['two_active_minimum_output_weight'], 6)
        self.assertEqual(target.checked_dependencies()['hd_lambda_multi_transvection.py'], target.CORE_SHA256)

    def test_invalid_inputs(self):
        for values in ((0,)*15, (True,)+(0,)*15, (256,)+(0,)*15):
            with self.assertRaises(ValueError):
                target.circuit(values)
        for lam, q in ((0, 256), (1, 256), (2, 8)):
            with self.assertRaises(ValueError):
                target.circuit((0,)*16, lam, q)


if __name__ == '__main__':
    unittest.main()

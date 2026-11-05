from collections import Counter
from itertools import product
import unittest

import hd_lambda_double_exact6 as exact
import hd_quaternary_outer_relations as old


class ExactSixTests(unittest.TestCase):
    def test_dependencies(self):
        exact.dependencies()

    def test_small_category_identity(self):
        for a in range(1, 5):
            basis, twice, masses = exact.relation_arrays(a)
            observed = Counter()
            for i, count in enumerate(masses):
                rows = tuple(int(row) for row in basis[:, i] if row)
                observed[rows] += int(count)
            expected = Counter({rows: count for rows, count in zip(old.kernel_subspaces(a),
                old.category_data(a)['category_counts']) if count})
            self.assertEqual(observed, expected)
            for form in range(4**a):
                flags = exact.zero_membership(form, basis, twice)
                bits = {rows: kernel for rows, kernel in zip(old.kernel_subspaces(a), old.kernel_bitsets(a))}
                self.assertEqual(list(flags), [bool(bits[tuple(int(row) for row in basis[:, i] if row)] >> form & 1)
                                              for i in range(len(masses))])

    def test_six_rows_toy_exhaustive(self):
        source = tuple(range(6))
        observed, _ = exact.exact_columns(source, q=4)
        inputs = list(product(range(4), repeat=6))
        for column in range(16):
            expected = Counter()
            forms = exact.scheduled_forms(source, column)
            for values in inputs:
                mask = 0
                for i, form in enumerate(forms):
                    output = 0
                    for k, value in enumerate(values):
                        output ^= old.multiply_gf4(form >> (2*k) & 3, value)
                    mask |= bool(output) << i
                expected[sum(bool(x) for x in values), mask] += 1
            self.assertEqual(observed[column], dict(expected))

    def test_gf256_against_small_exact_provider(self):
        for source in ((0,), (1, 9), (0, 3, 7), (2, 5, 8, 13)):
            columns, _ = exact.exact_columns(source, columns=(0, 5, 15))
            for j, counts in zip((0, 5, 15), columns):
                self.assertEqual(counts, old.joint_counts(exact.scheduled_forms(source, j), len(source)))

    def test_invalid_parameters(self):
        for args in ((0, 256), (7, 256), (6, 8), (True, 4)):
            with self.assertRaises(ValueError): exact.relation_arrays(*args)
        with self.assertRaises(ValueError): exact.exact_columns((0, 0))
        with self.assertRaises(ValueError): exact.exact_columns((0, 1), columns=(0, 0))


if __name__ == '__main__':
    unittest.main()

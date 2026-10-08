"""Small exact checks of outer averaging, not a whole-code certificate."""

from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product
from math import factorial, log2
import unittest

import outer_freshness_transfer as transfer


class OuterFreshnessTests(unittest.TestCase):
    def test_positive_row_automaton(self):
        for length in range(1, 7):
            for parity in range(length):
                fresh = [F(j + 2, j + 3) for j in range(length)]
                dependent = [F(j + 4, j + 2) for j in range(length)]
                expected = 0
                for pattern in transfer.valid_row_patterns(length, parity):
                    term = 1
                    for j, kind in enumerate(pattern):
                        if kind:
                            term *= fresh[j] if kind == 1 else dependent[j]
                    expected += term
                    if any(pattern):
                        self.assertEqual(pattern.count(1),
                                         sum(c != 0 for c in pattern) - parity)
                        self.assertEqual(pattern.count(2), parity)
                self.assertEqual(transfer.row_polynomial(fresh, dependent, parity),
                                 expected)

    def test_profile_coefficient_probability(self):
        rows, length, parity = 2, 3, 1
        valid = set(transfer.valid_row_patterns(length, parity))
        coefficients = Counter()
        for row_patterns in product(valid, repeat=rows):
            key = tuple((sum(p[j] == 1 for p in row_patterns),
                         sum(p[j] == 2 for p in row_patterns))
                        for j in range(length))
            coefficients[key] += 1
        perms = list(permutations(range(rows)))
        for columns in product(list(product(range(3), repeat=rows)), repeat=length):
            key = tuple((column.count(1), column.count(2)) for column in columns)
            good = 0
            for region_perms in product(perms, repeat=length):
                routed_rows = [tuple(columns[j][region_perms[j][r]]
                                     for j in range(length)) for r in range(rows)]
                good += all(row in valid for row in routed_rows)
            denominator = 1
            for fresh, dependent in key:
                denominator *= transfer.ternary_count(rows, fresh, dependent)
            self.assertEqual(F(good, factorial(rows)**length),
                             F(coefficients[key], denominator))

    @staticmethod
    def toy_matrix(pattern):
        # Placement dependence is intentional: counting only f,d must average.
        a, b = pattern
        return [[F(1 + a + 2*b, 5), F(1 + a*b, 7)],
                [F(2 + 2*a + b, 11), F(3 + a + b, 13)]]

    def test_packet_normalization(self):
        table = {}
        for pattern in product(range(3), repeat=2):
            key = pattern.count(1), pattern.count(2)
            matrix = table.setdefault(key, [[F(0), F(0)], [F(0), F(0)]])
            term = self.toy_matrix(pattern)
            for i in range(2):
                for j in range(2):
                    matrix[i][j] += term[i][j]
        for key, matrix in table.items():
            patterns = [p for p in product(range(3), repeat=2)
                        if (p.count(1), p.count(2)) == key]
            self.assertEqual(len(patterns), transfer.ternary_count(2, *key))
            for i in range(2):
                for j in range(2):
                    self.assertEqual(matrix[i][j] / len(patterns),
                                     sum(self.toy_matrix(p)[i][j] for p in patterns)
                                     / len(patterns))

    def test_profile_free_envelope_against_all_routes(self):
        rows, length, parity, width = 4, 2, 1, 2
        tau, x = F(2), F(2, 3)
        fresh_marks, dependent_marks = [F(3, 2), F(4, 3)], [F(2), F(5, 4)]
        table = {}
        for pattern in product(range(3), repeat=width):
            key = pattern.count(1), pattern.count(2)
            matrix = table.setdefault(key, [[F(0), F(0)], [F(0), F(0)]])
            term = self.toy_matrix(pattern)
            for i in range(2):
                for j in range(2):
                    matrix[i][j] += term[i][j]
        exact = F(0)
        valid = list(transfer.valid_row_patterns(length, parity))
        perms = list(permutations(range(rows)))
        for row_patterns in product(valid, repeat=rows):
            if not any(any(p) for p in row_patterns):
                continue
            fresh = sum(p.count(1) for p in row_patterns)
            weight = sum(sum(c != 0 for c in p) for p in row_patterns)
            route_total = F(0)
            for region_perms in product(perms, repeat=length):
                matrix = transfer.identity(2)
                for j, perm in enumerate(region_perms):
                    column = [row_patterns[r][j] for r in perm]
                    for start in range(0, rows, width):
                        matrix = transfer.matmul(matrix,
                            self.toy_matrix(column[start:start+width]))
                route_total += sum(matrix[0])
            exact += tau**fresh * x**weight * route_total / len(perms)**length
        bound = transfer.exact_bound(table, rows=rows, length=length,
            parity=parity, width=width, tau=tau, x=x, y=F(1),
            fresh_marks=fresh_marks, dependent_marks=dependent_marks,
            input_max=rows*length, output_max=rows*length)
        # Undo only the existence/weight-cut prefactors, leaving H^L-1 times
        # the matrix bound.  The exact left side is a weighted support sum.
        moment_bound = bound * tau * x**(rows*length)
        self.assertLessEqual(exact, moment_bound)
        self.assertGreater(exact, 0)
        scalar_bound = transfer.exact_bound(table, rows=rows, length=length,
            parity=parity, width=width, tau=tau, x=x, y=F(1),
            fresh_marks=fresh_marks, dependent_marks=dependent_marks,
            input_max=rows*length, output_max=rows*length,
            continuation=[F(1), F(3, 2)])
        self.assertLessEqual(exact, scalar_bound*tau*x**(rows*length))

    def test_exact_packet_adapter(self):
        import numpy as np
        from rank_packet_operator import packet_matrices
        exact = transfer.exact_rank_packet_table(F(3), F(2, 5), q=4, width=4, state=1)
        floating = packet_matrices(3.0, 0.4, q=4, t=4, s=1)
        for (fresh, dep), matrix in exact.items():
            self.assertTrue(np.allclose(np.array(matrix, dtype=float),
                                        floating[fresh, dep], rtol=2e-14, atol=2e-14))

    def test_floating_scaling_against_exact(self):
        params = dict(rows=8, length=3, parity=1, width=4,
                      tau=F(3), x=F(2, 3), y=F(2, 5),
                      fresh_marks=[F(3, 2), F(4, 3), F(5, 4)],
                      dependent_marks=[F(2), F(5, 4), F(3, 2)],
                      input_max=17, output_max=4)
        exact_table = transfer.exact_rank_packet_table(params['tau'], params['y'],
                                                     q=4, width=4, state=1)
        exact = transfer.exact_bound(exact_table, **params)
        float_params = {k: ([float(v) for v in value] if isinstance(value, list)
                            else float(value) if isinstance(value, F) else value)
                        for k, value in params.items()}
        actual = transfer.floating_log_bound(q=4, state=1, **float_params)
        self.assertAlmostEqual(actual, log2(exact), places=11)
        # A non-power-of-two packet count tests accumulated matrix scales.
        params['rows'] = 20
        exact = transfer.exact_bound(exact_table, **params)
        float_params['rows'] = 20
        actual = transfer.floating_log_bound(q=4, state=1, **float_params)
        self.assertAlmostEqual(actual, log2(exact), places=10)
        for law in ('unscaled', 'eta-scaled'):
            exact_table = transfer.exact_rank_packet_table(params['tau'], params['y'],
                                        q=4, width=4, state=1, law=law)
            exact = transfer.exact_bound(exact_table, continuation=[F(1), F(3, 2)],
                                         **params)
            actual = transfer.floating_log_bound(q=4, state=1, law=law,
                                continuation=[1., 1.5], **float_params)
            self.assertAlmostEqual(actual, log2(exact), places=10)

    def test_empty_source_removed_before_envelope(self):
        self.assertEqual(transfer.row_polynomial([F(0)]*3, [F(2)]*3, 1), 1)
        table = {(0, 0): [[F(2), F(3)], [F(5), F(7)]]}
        envelope = transfer.packet_envelope(table, 2, F(1), F(1))
        self.assertEqual(envelope, table[0, 0])
        # This deliberately loose matrix allows fictitious empty-input paths.
        # Subtracting only one AFTER multiplication would not remove them.
        self.assertGreater(sum(transfer.matpow(envelope, 4)[0]), 1)

    def test_worst_modal_conditioning(self):
        from math import factorial
        for rows in range(1, 17):
            for cap in range(rows+1):
                costs = []
                for f in range(cap+1):
                    for d in range(cap-f+1):
                        z = rows-f-d
                        costs.append(F(rows**rows*factorial(f)*factorial(d)*factorial(z),
                                       factorial(rows)*f**f*d**d*z**z))
                self.assertEqual(transfer.exact_conditioning_cost(rows, cap), max(costs))

    def test_bernstein_thinning_identity(self):
        for width in range(1, 6):
            coefficients = {(f, d): F(1+f*f+3*d, 2+f+d)
                            for f in range(width+1) for d in range(width-f+1)}
            density = F(2, 7)
            thin = transfer.thinned_coefficients(coefficients, width, density)
            for a, b in ((F(0), F(0)), (F(1), F(0)), (F(1, 3), F(1, 4))):
                def evaluate(table, u, v):
                    return sum(transfer.ternary_count(width, f, d)*value
                               *u**f*v**d*(1-u-v)**(width-f-d)
                               for (f,d), value in table.items())
                exact = evaluate(coefficients, density*a, density*b)
                self.assertEqual(exact, evaluate(thin, a, b))
                self.assertLessEqual(exact, max(thin.values()))

    def test_active_row_iid_comparison(self):
        rows, active, width = 4, 2, 2
        potential = [F(1), F(3, 2)]
        marks = [(F(3, 2), F(5, 4)), (F(4, 3), F(7, 5))]
        table = {}
        for pattern in product(range(3), repeat=width):
            key = pattern.count(1), pattern.count(2)
            matrix = table.setdefault(key, [[F(0), F(0)], [F(0), F(0)]])
            term = self.toy_matrix(pattern)
            for i in range(2):
                for j in range(2):
                    matrix[i][j] += term[i][j]
        bound = F(potential[0], min(potential))*transfer.exact_conditioning_cost(rows, active)**2
        for region, (v, w) in enumerate(marks):
            rho = 0
            for r in range(2):
                coefficients = {key: sum(row[r][k]*potential[k] for k in range(2))
                    /(potential[r]*transfer.ternary_count(width, *key)
                      *v**key[0]*w**key[1]) for key, row in table.items()}
                thin = transfer.thinned_coefficients(coefficients, width, F(active, rows))
                rho = max(rho, *(value for (f, d), value in thin.items()
                                 if (d == 0 if region == 0 else f == 0)))
            bound *= rho**(rows//width)
        perms = list(permutations(range(rows)))
        row_patterns = [(1, 2)]*active+[(0, 0)]*(rows-active)
        exact = 0
        for route in product(perms, repeat=2):
            matrix = transfer.identity(2)
            for region, perm in enumerate(route):
                column = [row_patterns[r][region] for r in perm]
                for start in range(0, rows, width):
                    matrix = transfer.matmul(matrix,
                        self.toy_matrix(column[start:start+width]))
            exact += sum(matrix[0])
        exact /= len(perms)**2*marks[0][0]**active*marks[1][1]**active
        self.assertLessEqual(exact, bound)


if __name__ == '__main__':
    unittest.main()

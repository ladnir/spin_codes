"""Exact positive mass, noncommuting order, floors, scope and sharing tests."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_two_sum_affine_sparse_certificate as target
from hd_two_sum_affine_outer import column_action_exact, sparse_column_coefficients


def convert(regions):
    return [[tuple(target.rational(x) for x in matrix) for matrix in table] for table in regions]


def act(matrix, vector):
    a, b, c, d = matrix
    x, y = vector
    return a*x+b*y, c*x+d*y


def encloses_fraction(ball, value):
    def exact(endpoint):
        mantissa, exponent = endpoint.man_exp()
        return F(int(mantissa))*F(2)**int(exponent)
    return exact(ball.lower()) <= value <= exact(ball.upper())


class AffineSparseCertificateTests(unittest.TestCase):
    def test_shared_recurrence_matches_exact_labelled_noncommuting_words(self):
        regions = [[(F(1+i+m), F(2+i), F(1+m*m), F(3+2*i+m)) for m in range(4)]
                   for i in range(8)]
        potential = F(2, 5)
        with ctx.workprec(256):
            for total in range(4):
                vectors = target.hole_vectors_arb(convert(regions), target.rational(potential), total)
                for double, vector in enumerate(vectors):
                    expected = [F(0), F(0)]
                    for holes in product(range(-1, 8), repeat=double):
                        probability = F(247**holes.count(-1), 255**double)
                        v = (F(1), potential)
                        for i in reversed(range(8)):
                            v = act(regions[i][total-holes.count(i)], v)
                        for state in range(2):
                            expected[state] += probability*v[state]
                    for state in range(2):
                        self.assertTrue(encloses_fraction(vector[state], expected[state]))
                    self.assertEqual(tuple(expected), column_action_exact(regions, (1, potential), total-double, double))
            self.assertNotEqual(column_action_exact(regions, (1, potential), 1, 2),
                                column_action_exact(list(reversed(regions)), (1, potential), 1, 2))

    def test_hole_mass_normalizes_and_caps_follow_mixture(self):
        regions = [[(F(1+i), F(2), F(3), F(4+i))]*4 for i in range(8)]
        potential = F(1, 3)
        expected = (F(1), potential)
        for table in reversed(regions):
            expected = act(table[0], expected)
        with ctx.workprec(256):
            matrices = convert(regions)
            caps = target.category_caps_arb(matrices, target.rational(potential), 3)
            for total in range(4):
                for vector in target.hole_vectors_arb(matrices, target.rational(potential), total):
                    for state in range(2):
                        self.assertTrue(encloses_fraction(vector[state], expected[state]))
            exactcap = max(expected[0], expected[1]/potential)
            for cap in caps.values():
                self.assertGreaterEqual(cap, target.rational(exactcap).upper())

    def test_source_matching_envelopes_keep_upper_endpoint_direction(self):
        raw = {(m, f): tuple(arb(10*m+f+entry, '0.01') for entry in range(4))
               for m in range(4) for f in range(m//2+1)}
        with ctx.workprec(256):
            tables = target.regional_envelopes(raw, 3)
            self.assertEqual(tables[0][3], tuple(value.upper() for value in raw[3, 0]))
            self.assertEqual(tables[1][3], tuple(value.upper() for value in raw[3, 1]))
            self.assertEqual(tables[1][2], tuple(value.upper() for value in raw[2, 1]))

    def test_positive_polynomial_and_exact_support_floor(self):
        caps = {(s, d): F(2+s+3*d, 7) for s in range(3) for d in range(3) if s+d <= 2}
        q, n1, n2, columns, dimension = 4, 1, 1, 3, 2
        exact = sparse_column_coefficients(caps, n1, n2, q=q)
        coefficients = [F(0)]*(3*columns+1)
        for word in product(range(len(exact)), repeat=columns):
            value = F(1)
            for degree in word:
                value *= exact[degree]
            coefficients[sum(word)] += value
        expected = sum(coefficients[6:])/(q-1)**3
        with ctx.workprec(256):
            directed = {key: target.rational(value) for key, value in caps.items()}
            actual = target.sparse_coefficients_arb(directed, n1, n2, q=q)
            for left, right in zip(actual, exact):
                self.assertTrue(encloses_fraction(left, right))
            result = target.sparse_tail_arb(directed, n1, n2, columns=columns, dimension=dimension, q=q)
            self.assertTrue(encloses_fraction(result, expected))

    def test_frozen_scope_and_mock_replay_share_every_required_bank(self):
        counts = target.checked_proposal()
        self.assertEqual(len(counts), 81)
        self.assertEqual(sum(counts.values()), comb(8192, 17))
        with patch.object(target, 'arb_pair_regions', return_value=None) as raw, \
                patch.object(target, 'regional_envelopes', return_value=[None]*9) as envelopes, \
                patch.object(target, 'category_caps_arb', return_value={}) as caps, \
                patch.object(target, 'sparse_tail_arb', return_value=arb(1)) as tails:
            result = target.receipt(256)
        self.assertEqual(raw.call_count, 1)
        self.assertEqual(envelopes.call_count, 1)
        self.assertEqual(caps.call_count, 81)
        self.assertEqual(tails.call_count, 1458)
        self.assertEqual(len(result['entries']), 1458)
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['complete_dense_cohort_coverage'])
        self.assertFalse(result['saved_float_estimates_used'])
        with ctx.workprec(512):
            p = arb(2)**target.rational(target.LOG_POTENTIAL)
            expected = (comb(8192, 17)*36**17*target.rational(target.MARK)**(-target.LIMIT)/(255*p))
            self.assertEqual(F(result['failure_log2_upper']), target.ceiling(expected.log()/arb(2).log()))

    def synthetic_receipt(self, precision=256):
        counts = target.checked_proposal()
        data = dict(target.EXPECTED_SCOPE, precision_bits=precision,
            frozen_proposal_sha256=target.PROPOSAL_SHA256, reviewed_helper_sha256=target.HELPER_PINS,
            source_sha256=sha256(Path(target.__file__).read_bytes()).hexdigest(),
            entries=[dict(target.cell_scope(source, count, n1), failure_log2_upper='-100')
                     for source, count in counts.items() for n1 in range(18)])
        return data, counts

    def test_checked_cells_fail_closed_on_scope_duplicates_and_counts(self):
        data, counts = self.synthetic_receipt()
        self.assertEqual(len(target.checked_rows(data, 256, counts)), 1458)
        changed = deepcopy(data)
        changed['whole_code_certificate'] = True
        with self.assertRaises(ValueError):
            target.checked_rows(changed, 256, counts)
        for change in ('duplicate', 'missing', 'rowchoices', 'sourcecount'):
            changed = deepcopy(data)
            if change == 'duplicate':
                changed['entries'][-1] = changed['entries'][0]
            elif change == 'missing':
                changed['entries'].pop()
            elif change == 'rowchoices':
                changed['entries'][0]['pair_row_choices'] *= 28
            else:
                changed['entries'][0]['source_supports'] += 1
            with self.assertRaises(ValueError):
                target.checked_rows(changed, 256, counts)

    def test_independent_union_and_budget_guard(self):
        result = target.rounded_union_512([F(-100)]*4)
        self.assertGreaterEqual(result, -98)
        self.assertLessEqual(result, -98+F(1, 10**12))
        for precision in (64, 512, True):
            with self.assertRaises(ValueError):
                target.receipt(precision)
        for budget in (0, 601, True):
            with self.assertRaises(ValueError):
                target.receipt(max_seconds=budget)
        with patch.object(target, 'arb_pair_regions', return_value=None), \
                patch.object(target, 'regional_envelopes', return_value=[None]*9), \
                patch.object(target, 'perf_counter', side_effect=[0., 601.]):
            with self.assertRaises(TimeoutError):
                target.receipt()


if __name__ == '__main__':
    unittest.main()

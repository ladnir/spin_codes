"""Exact small ordered-matrix controls; no production regional bank."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, ctx

import hd_two_sum_gf16_single_dense_certificate as target
from hd_two_sum_affine_sparse_certificate_tests import encloses_fraction


def flat(matrix):
    return tuple(target.rational(x) for x in matrix.ravel())


class SingleDenseCertificateTests(unittest.TestCase):
    def tables(self):
        return [np.array([[[1+i, 2+i*i], [3, 1+2*i]], [[2+i, 1], [1+3*i, 4+i*i]]], dtype=object)
                for i in range(4)]

    def test_every_labelled_mask_ordered_weight_coefficient(self):
        tables = self.tables()
        with ctx.workprec(256):
            got = target.weight_matrices_arb([tuple(map(flat, table)) for table in tables])
            expected = [np.zeros((2, 2), dtype=object) for _ in range(5)]
            for bits in product((0, 1), repeat=4):
                value = np.eye(2, dtype=object)
                for table, bit in zip(tables, bits):
                    value = value @ table[bit]
                expected[sum(bits)] += value
            for actual, exact in zip(got, expected):
                self.assertTrue(all(encloses_fraction(x, F(y)) for x, y in zip(actual, exact.ravel())))
            reverse = target.weight_matrices_arb([tuple(map(flat, table)) for table in tables[::-1]])
            self.assertFalse(got[1][0].overlaps(reverse[1][0]))

    def test_positive_matrix_power_matches_all_noncommuting_words(self):
        branches = np.array([[[2, 1], [1, 3]], [[1, 4], [2, 1]]], dtype=object)
        with ctx.workprec(256):
            actual, degrees = target.matrix_tail_arb(tuple(map(flat, branches)), 5, 3)
            expected, per_degree = 0, [0]*6
            for bits in product((0, 1), repeat=5):
                value = np.eye(2, dtype=object)
                for bit in bits:
                    value = value @ branches[bit]
                if sum(bits) >= 3:
                    contribution = value[0, 0]+value[0, 1]
                    expected += contribution
                    per_degree[sum(bits)] += contribution
            self.assertTrue(encloses_fraction(actual, F(expected)))
            for got, exact in zip(degrees, per_degree[3:]):
                self.assertTrue(encloses_fraction(got, F(exact)))

    def test_initial_zero_row_and_no_endpoint_potential(self):
        with ctx.workprec(256):
            branches = ((arb(2), arb(0), arb(0), arb(2)), (arb(3), arb(0), arb(0), arb(3)))
            actual, _ = target.matrix_tail_arb(branches, 5, 3)
            exact = sum(comb(5, j)*2**(5-j)*3**j for j in range(3, 6))
            self.assertTrue(encloses_fraction(actual, F(exact)))
            enlarged = tuple((a, b, c, arb(2)**100) for a, b, c, _ in branches)
            again, _ = target.matrix_tail_arb(enlarged, 5, 3)
            self.assertTrue(encloses_fraction(again, F(exact)))
            self.assertEqual(actual, again)

    def test_holder_weights_are_mixed_before_matrix_power(self):
        matrices = [tuple(map(arb, row)) for row in ((1, 2, 3, 4), (2, 1, 4, 3), (5, 6, 7, 8))]
        got = target.matrix_branches_arb(matrices, {0: arb(2), 1: arb(3), 2: arb(5)})
        expected = ((2, 4, 6, 8), (31, 33, 47, 49))
        for matrix, exact in zip(got, expected):
            self.assertTrue(all(encloses_fraction(x, F(y)) for x, y in zip(matrix, exact)))

    def test_three_banks_fourteen_arities_and_exact_once_only_prefix(self):
        self.assertEqual(target.checked_proposal(), target.SCHEDULE)
        with patch.object(target, 'arity_weights_arb', return_value=({0: arb(1)}, 'toy')), \
                patch.object(target, 'arb_pair_regions', return_value={(0, 0): None, (1, 0): None}) as bank, \
                patch.object(target, 'weight_matrices_arb', return_value=None) as weights, \
                patch.object(target, 'matrix_branches_arb', return_value=None), \
                patch.object(target, 'matrix_tail_arb', return_value=(arb(1), [])) as tails:
            data = target.receipt()
        self.assertEqual([call.args for call in bank.call_args_list], [(z, 1) for z in target.MARKS])
        self.assertEqual(weights.call_count, 3)
        self.assertEqual(tails.call_count, 14)
        self.assertEqual(data['endpoint_potential_factors'], 0)
        self.assertFalse(data['common_potential_used'])
        self.assertFalse(data['whole_code_certificate'])
        with ctx.workprec(512):
            expected = sum((8192*target.rational(target.SCHEDULE[a])**(-734003)/255
                            for a in target.ARITIES), arb(0))
            self.assertEqual(F(data['failure_log2_upper']), target.ceiling(expected.log()/arb(2).log()))

    def fixture(self):
        return dict(target.EXPECTED_SCOPE, precision_bits=256,
            frozen_proposal_sha256=target.PROPOSAL_SHA256, reviewed_helper_sha256=target.HELPER_PINS,
            proposal_source_sha256=target.SOURCE_PINS,
            source_sha256=sha256(Path(target.__file__).read_bytes()).hexdigest(),
            entries=[dict(active_constituents=a, mark=str(target.SCHEDULE[a]), fixed_row_sets=comb(16, a),
                exact_moments_sha256=target.exact_moments(a)[2], root_upper_checks_passed=True,
                positive_weight_states=sum(n > 0 for n in target.exact_moments(a)[0]),
                failure_log2_upper='-100') for a in target.ARITIES])

    def test_scope_selection_and_complete_partition_fail_closed(self):
        data = self.fixture()
        self.assertEqual(len(target.checked_rows(data, 256)), 14)
        for key, value in (('whole_code_certificate', True), ('endpoint_potential_factors', 1),
                           ('initial_state', 1), ('cohort_counts', [0, 1, 0])):
            with self.assertRaises(ValueError):
                target.checked_rows(dict(data, **{key: value}), 256)
        changed = deepcopy(data)
        changed['entries'][0]['mark'] = '1199/1200'
        with self.assertRaises(ValueError):
            target.checked_rows(changed, 256)
        changed['entries'][0] = changed['entries'][1]
        with self.assertRaises(ValueError):
            target.checked_rows(changed, 256)
        with patch.object(target, 'PROPOSAL_SHA256', '0'*64), self.assertRaises(ValueError):
            target.checked_proposal()

    def test_fresh_rounded_union_and_dual_disagreement(self):
        rows = {a: F(-100-a) for a in target.ARITIES}
        expected = target.rounded_union_512(list(rows.values()))
        encoded = json.dumps(dict(failure_log2_upper=str(expected-F(1, 10**12)), elapsed_seconds=0.)).encode()
        with patch.object(Path, 'read_bytes', return_value=encoded), \
                patch.object(target, 'checked_rows', return_value=rows):
            data = target.verify()
        self.assertEqual(F(data['rounded_arity_union_log2_upper']), expected)
        with patch.object(Path, 'read_bytes', return_value=encoded), \
                patch.object(target, 'checked_rows', side_effect=[rows, dict(rows, bad=F(-100))]), \
                self.assertRaises(ValueError):
            target.verify()

    def test_zero_matrices_guards_and_budget(self):
        zero = (arb(0),)*4
        value, _ = target.matrix_tail_arb((zero, zero), 3, 2)
        self.assertTrue(value.is_zero())
        for p in (128, 512, True):
            with self.assertRaises(ValueError):
                target.receipt(p)
        for budget in (True, 0, 121, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                target.receipt(max_seconds=budget)
        with patch.object(target, 'arb_pair_regions') as bank, self.assertRaises(TimeoutError):
            target.receipt(max_seconds=1e-12)
        bank.assert_not_called()
        with self.assertRaises(ValueError):
            target.weight_matrices_arb([((arb(-1),)*4, zero)])


if __name__ == '__main__':
    unittest.main()

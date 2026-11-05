"""Exact tiny tests for the g16 cardinality replay; no production bank builds."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_two_sum_gf16_arity_certificate as target
from hd_two_sum_affine_sparse_certificate_tests import encloses_fraction
from hd_xor_dense_holder_certificate import dyadic_power_dominates


class Gf16ArityCertificateTests(unittest.TestCase):
    def tables(self):
        return [[(F(1+i+m), F(2+i), F(1+m*m), F(3+2*i+m)) for m in range(4)]
                for i in range(3)]

    def directed(self, tables):
        return [[tuple(target.rational(x) for x in matrix) for matrix in table] for table in tables]

    def test_every_tiny_labelled_hole_and_dense_mask_assignment(self):
        tables, p, q = self.tables(), F(2, 3), 8
        with ctx.workprec(256):
            vectors = target.weight_vectors_arb(self.directed(tables), target.rational(p), 2, q=q)
            for (single, double), weights in vectors.items():
                expected = [[F(0), F(0)] for _ in range(4)]
                for mask in range(8):
                    for holes in product(range(-1, 3), repeat=double):
                        probability = F(4**holes.count(-1), 7**double)
                        value = (F(1), p)
                        for i in reversed(range(3)):
                            a, b, c, d = tables[i][single+double-holes.count(i)+((mask >> i) & 1)]
                            x, y = value
                            value = a*x+b*y, c*x+d*y
                        for state in range(2):
                            expected[mask.bit_count()][state] += probability*value[state]
                for w in range(4):
                    for state in range(2):
                        self.assertTrue(encloses_fraction(weights[w][state], expected[w][state]))
            backward = target.weight_vectors_arb(self.directed(tables[::-1]), target.rational(p), 2, q=q)
            self.assertFalse(vectors[0, 2][1][0].overlaps(backward[0, 2][1][0]))

    def test_identity_has_one_binomial_and_zero_is_preserved(self):
        with ctx.workprec(256):
            identity = [[(arb(1), arb(0), arb(0), arb(1))]*4]*4
            vectors = target.weight_vectors_arb(identity, target.rational(F(2, 3)), 2, q=8)
            for weights in vectors.values():
                for w, vector in weights.items():
                    self.assertTrue(encloses_fraction(vector[0], F(comb(4, w))))
                    self.assertTrue(encloses_fraction(vector[1], F(2, 3)*comb(4, w)))
            zeros = [[(arb(0),)*4]*4]*3
            values = target.weight_vectors_arb(zeros, arb(1), 2, q=8)
            self.assertTrue(all(x.is_zero() for table in values.values() for vector in table.values() for x in vector))

    def test_all14_exact_sixteenth_root_endpoints(self):
        with ctx.workprec(256):
            for arity in target.ARITIES:
                nums, denominator, digest = target.exact_moments(arity)
                weights, actual_digest = target.arity_weights_arb(arity)
                self.assertEqual(actual_digest, digest)
                self.assertEqual(denominator, 255**(8*arity))
                self.assertEqual(nums[0], comb(16, arity))
                self.assertEqual(set(weights), {w for w, n in enumerate(nums) if n})
                for w, upper in weights.items():
                    self.assertTrue(dyadic_power_dominates(upper, nums[w], denominator, 16))
            # Full arity is exactly ambient uniform comparison, including its density.
            nums, den, _ = target.exact_moments(16)
            self.assertEqual(nums[5]*255**48, den)

    def test_one_shared_bank_and_fourteen_tails_with_exact_prefix(self):
        metadata = target.checked_proposal()
        with patch.object(target, 'arity_weights_arb', return_value=({0: arb(1)}, 'toy')), \
                patch.object(target, 'arb_pair_regions', return_value=None) as raw, \
                patch.object(target, 'regional_envelopes', return_value=[None]*9), \
                patch.object(target, 'weight_vectors_arb', return_value=None) as vectors, \
                patch.object(target, 'holder_branch_caps_arb', return_value=({}, {})) as quotient, \
                patch.object(target, 'column_coefficients_arb', return_value=[arb(1)]), \
                patch.object(target, 'bivariate_tail_arb', return_value=(arb(1), [])) as tails:
            result = target.receipt(256)
        raw.assert_called_once_with(target.MARK, 17)
        self.assertEqual(vectors.call_count, 1)
        self.assertEqual(len(vectors.call_args.args[0]), 16)
        self.assertEqual(quotient.call_count, 14)
        self.assertEqual(tails.call_count, 14)
        for call in tails.call_args_list:
            self.assertEqual(call.kwargs, dict(columns=16, x_floor=288, y_floor=9))
        self.assertEqual(result['fixed_dense_row_sets'], 65399)
        self.assertEqual(result['sparse_pair_assignments'], 120**16)
        self.assertFalse(result['extra_target_mask_binomial'])
        self.assertFalse(result['extra_dense_beta_factor'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['mixer_matrix_sha256'], metadata['mixer_matrix_sha256'])
        with ctx.workprec(512):
            p = arb(2)**target.rational(target.LOG_POTENTIAL)
            expected = 14*target.SOURCE_COUNT*17*120**16*target.rational(target.MARK)**(-734003)/(255*p*255**256)
            self.assertEqual(F(result['failure_log2_upper']), target.ceiling(expected.log()/arb(2).log()))

    def receipt_fixture(self):
        return dict(target.EXPECTED_SCOPE, **target.checked_proposal(), precision_bits=256,
            frozen_proposal_sha256=target.PROPOSAL_SHA256, reviewed_helper_sha256=target.HELPER_PINS,
            source_sha256=sha256(Path(target.__file__).read_bytes()).hexdigest(),
            entries=[dict(active_constituents=a, fixed_row_sets=comb(16, a),
                exact_moments_sha256=target.exact_moments(a)[2], root_upper_checks_passed=True,
                positive_weight_states=sum(n > 0 for n in target.exact_moments(a)[0]),
                failure_log2_upper='-100') for a in target.ARITIES])

    def test_scope_partition_and_pins_fail_closed(self):
        data = self.receipt_fixture()
        self.assertEqual(len(target.checked_rows(data, 256)), 14)
        for key, value in (('whole_code_certificate', True), ('sparse_pair_assignments', 120),
                ('source_matching_counts', [0, 8]), ('holder_exponent', 32),
                ('extra_target_mask_binomial', True), ('mark', '1/2')):
            with self.assertRaises(ValueError):
                target.checked_rows(dict(data, **{key: value}), 256)
        changed = deepcopy(data)
        changed['entries'][-1] = changed['entries'][0]
        with self.assertRaises(ValueError):
            target.checked_rows(changed, 256)
        with patch.object(target, 'PROPOSAL_SHA256', '0'*64), self.assertRaises(ValueError):
            target.checked_proposal()

    def test_verifier_freshly_sums_rounded_classes_and_rejects_disagreement(self):
        rows = {a: F(-100-a, 1) for a in target.ARITIES}
        expected = target.rounded_union_512(list(rows.values()))
        data = dict(failure_log2_upper=str(expected-F(1, 10**12)), elapsed_seconds=0.)
        encoded = json.dumps(data).encode()
        with patch.object(target, 'checked_proposal', return_value={}), \
                patch.object(target, 'checked_rows', return_value=rows), \
                patch.object(Path, 'read_bytes', return_value=encoded):
            verified = target.verify()
        self.assertEqual(F(verified['rounded_class_union_log2_upper']), expected)
        self.assertEqual(verified['independent_sum_precision_bits'], 512)
        with patch.object(target, 'checked_proposal', return_value={}), \
                patch.object(target, 'checked_rows', side_effect=[rows, dict(rows, **{'bad': F(-100)})]), \
                patch.object(Path, 'read_bytes', return_value=encoded), self.assertRaises(ValueError):
            target.verify()

    def test_budget_precision_and_negative_matrix_guards(self):
        for precision in (True, 128, 512):
            with self.assertRaises(ValueError):
                target.receipt(precision)
        for budget in (True, 0, 601, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                target.receipt(max_seconds=budget)
        with patch.object(target, 'arb_pair_regions') as bank, self.assertRaises(TimeoutError):
            target.receipt(max_seconds=1e-12)
        bank.assert_not_called()
        for sparse in (True, -1, 17):
            with self.assertRaises(ValueError):
                target.weight_vectors_arb(self.directed(self.tables()), arb(1), sparse, q=8)
        invalid = self.directed(self.tables())
        invalid[0][0] = (arb(-1), arb(0), arb(0), arb(1))
        with self.assertRaises(ValueError):
            target.weight_vectors_arb(invalid, arb(1), 2, q=8)


if __name__ == '__main__':
    unittest.main()

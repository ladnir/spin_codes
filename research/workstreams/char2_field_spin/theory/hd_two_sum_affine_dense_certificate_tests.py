"""Small exact controls for labelled affine holes and the six-arity replay."""
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

import hd_two_sum_affine_dense_certificate as target
from hd_two_sum_affine_outer import column_action_exact
from hd_two_sum_affine_sparse_certificate_tests import encloses_fraction
from hd_xor_dense_holder_certificate import dyadic_power_dominates


class AffineDenseCertificateTests(unittest.TestCase):
    def tables(self):
        return [[(F(1+i+m), F(2+i), F(1+m*m), F(3+2*i+m)) for m in range(4)]
                for i in range(3)]

    def test_all_tiny_labelled_masks_against_exact_ordered_mixture(self):
        tables, p = self.tables(), F(2, 3)
        with ctx.workprec(256):
            directed = [[tuple(target.rational(value) for value in matrix) for matrix in table] for table in tables]
            vectors = target.dense_mask_vectors_arb(directed, target.rational(p), 2, q=8)
            for (single, double), masks in vectors.items():
                for mask, vector in masks.items():
                    shifted = [[tables[i][m+((mask >> i) & 1)] for m in range(3)] for i in range(3)]
                    expected = column_action_exact(shifted, (F(1), p), single, double, q=8)
                    for state in range(2):
                        self.assertTrue(encloses_fraction(vector[state], expected[state]))
            reversed_vectors = target.dense_mask_vectors_arb(list(reversed(directed)), target.rational(p), 2, q=8)
            self.assertFalse(vectors[0, 2][0][0].overlaps(reversed_vectors[0, 2][0][0]))

    def test_every_labelled_assignment_includes_dense_offset_and_nohole_mass(self):
        tables, p = self.tables(), F(2, 3)
        with ctx.workprec(256):
            directed = [[tuple(target.rational(value) for value in matrix) for matrix in table] for table in tables]
            vectors = target.dense_mask_vectors_arb(directed, target.rational(p), 2, q=8)
            for mask in range(8):
                expected = [F(0), F(0)]
                for holes in product(range(-1, 3), repeat=2):
                    weight = F(4**holes.count(-1), 7**2)
                    vector = (F(1), p)
                    for i in reversed(range(3)):
                        a, b, c, d = tables[i][2-holes.count(i)+((mask >> i) & 1)]
                        x, y = vector
                        vector = a*x+b*y, c*x+d*y
                    for state in range(2):
                        expected[state] += weight*vector[state]
                for state in range(2):
                    self.assertTrue(encloses_fraction(vectors[0, 2][mask][state], expected[state]))
            identity = [[(arb(1), arb(0), arb(0), arb(1))]*4]*3
            unit = target.dense_mask_vectors_arb(identity, target.rational(p), 2, q=8)
            for masks in unit.values():
                for vector in masks.values():
                    self.assertTrue(encloses_fraction(vector[0], F(1)))
                    self.assertTrue(encloses_fraction(vector[1], p))

    def test_new_affine_moments_and_exact_root_endpoint_checks(self):
        masks = (7, 11)
        numerators, denominator = target.affine_class_moments(masks)
        with ctx.workprec(256):
            weights, digest = target.class_weights_arb(masks)
            self.assertEqual(digest, sha256(json.dumps([list(numerators), denominator], separators=(',', ':')).encode()).hexdigest())
            for mask, upper in weights.items():
                self.assertTrue(dyadic_power_dominates(upper, numerators[mask], denominator))
            self.assertEqual(numerators[0]*255**48, 2*denominator)

    def test_one_bank_six_classes_and_count_factors(self):
        metadata = target.checked_proposal()
        with patch.object(target, 'class_weights_arb', return_value=({0: arb(1)}, 'toy')), \
                patch.object(target, 'arb_pair_regions', return_value=None) as raw, \
                patch.object(target, 'regional_envelopes', return_value=[None]*9), \
                patch.object(target, 'dense_mask_vectors_arb', return_value=None) as vectors, \
                patch.object(target, 'holder_branch_caps_arb', return_value=({}, {})) as quotient, \
                patch.object(target, 'column_coefficients_arb', return_value=[arb(1)]), \
                patch.object(target, 'bivariate_tail_arb', return_value=(arb(1), [])) as tails:
            result = target.receipt(256)
        self.assertEqual(raw.call_count, 1)
        self.assertEqual(vectors.call_count, 1)
        self.assertEqual(quotient.call_count, 6)
        self.assertEqual(tails.call_count, 6)
        for call in tails.call_args_list:
            self.assertEqual(call.kwargs, dict(columns=32, x_floor=544, y_floor=17))
        self.assertEqual(result['fixed_row_sets'], 219)
        self.assertEqual(result['source_colorings'], 17)
        self.assertEqual(result['sparse_pair_assignments'], 28**16)
        self.assertFalse(result['extra_dense_beta_factor'])
        self.assertFalse(result['old_xor_mask_counts_used'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['affine_matrix_sha256'], metadata['affine_matrix_sha256'])
        with ctx.workprec(512):
            p = arb(2)**target.rational(target.LOG_POTENTIAL)
            expected = 6*target.SOURCE_COUNT*17*28**16*target.rational(target.MARK)**(-734003)/(255*p*255**512)
            self.assertEqual(F(result['failure_log2_upper']), target.ceiling(expected.log()/arb(2).log()))

    def test_strict_receipt_scope_and_complete_arity_partition(self):
        numerators, denominator = (1,)*256, 1
        digest = sha256(json.dumps([list(numerators), denominator], separators=(',', ':')).encode()).hexdigest()
        data = dict(target.EXPECTED_SCOPE, **target.checked_proposal(), precision_bits=256,
            frozen_proposal_sha256=target.PROPOSAL_SHA256, reviewed_helper_sha256=target.HELPER_PINS,
            source_sha256=sha256(Path(target.__file__).read_bytes()).hexdigest(),
            entries=[dict(active_constituents=arity, active_row_masks=list(masks), fixed_row_sets=comb(8, arity),
                exact_moments_sha256=digest, root_upper_checks_passed=True,
                positive_mask_weights=256, failure_log2_upper='-100') for arity, masks in target.ARITY_CLASSES])
        with patch.object(target, 'affine_class_moments', return_value=(numerators, denominator)):
            self.assertEqual(len(target.checked_rows(data, 256)), 6)
            for key, value in (('whole_code_certificate', True), ('sparse_pair_assignments', 28),
                               ('source_matching_counts', [0, 8]), ('mark', '1/2')):
                changed = dict(data, **{key: value})
                with self.assertRaises(ValueError):
                    target.checked_rows(changed, 256)
            changed = deepcopy(data)
            changed['entries'][-1] = changed['entries'][0]
            with self.assertRaises(ValueError):
                target.checked_rows(changed, 256)

    def test_precision_and_budget_guards(self):
        for precision in (128, 512, True):
            with self.assertRaises(ValueError):
                target.receipt(precision)
        for budget in (0, 601, True):
            with self.assertRaises(ValueError):
                target.receipt(max_seconds=budget)


if __name__ == '__main__':
    unittest.main()

"""Tiny exact/Arb sums and mocked selected receipts; no coefficient gate."""

from fractions import Fraction as F
from hashlib import sha256
import json
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_gamma_exact_refinement_certificate as cert
from hd_field128_gamma_near_face_certificate_tests import fixture, as_arb, multiply, power, tiny_bank


class DirectedRefinementTests(unittest.TestCase):
    def test_disabled_refinement_replays_frozen_directed_sum(self):
        with ctx.workprec(256):
            for pairs, weight in ((3, 3), (4, 4), (5, 6)):
                regions = as_arb(fixture(pairs, weight))
                old, unused = cert.frozen.transfer_sum(pairs, weight, regions, repetitions=8, max_singletons=2)
                value, report = cert.transfer_sum(pairs, weight, regions, repetitions=8,
                                                  max_singletons=2, max_exact_profiles=0)
                self.assertEqual(old, value)
                self.assertEqual(report['profile_pairs_evaluated'], len(regions)**2)

    def test_directed_contains_independent_rational_full_sum(self):
        with ctx.workprec(256):
            for pairs, weight, repetitions in ((4, 4, 8), (5, 6, 8), (3, 3, 256)):
                rows = fixture(pairs, weight)
                counts = cert.refinement.ExactRefinedCounts(pairs, weight, exact_singletons=2*pairs)
                oracle = F(0)
                for a, left in rows.items():
                    for b, right in rows.items():
                        matrix = power(multiply(left, right), repetitions//2)
                        oracle += counts.count(a, b)*(matrix[0]+matrix[1])
                value, report = cert.transfer_sum(pairs, weight, as_arb(rows), counts=counts, repetitions=repetitions)
                self.assertTrue(value >= cert.base.rational(oracle))
                self.assertTrue(value < cert.base.rational(oracle)*(1+arb(2)**-180))
                self.assertIsNone(report['pending_even_full_range'])

    def test_complete_mock_receipt_context_and_exact_selection(self):
        bank, saved, contexts = tiny_bank(), [], []
        counts = cert.refinement.ExactRefinedCounts(cert.base.PAIRS, 2, selected_profiles=[(0, 0)], exact_singletons=0)
        with patch.object(cert.frozen, 'regional_matrices', return_value=as_arb(fixture(cert.base.PAIRS, 2))):
            result = cert.weight_receipt(bank, 2, counts=counts, checkpoint=saved.append, context_checkpoint=contexts.append)
        self.assertEqual(bank.pin_checks, 2)
        self.assertTrue(result['complete_selected_weight'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['profile_pairs_evaluated'], 4)
        self.assertEqual(result['count_provenance']['selected_canonical_profiles'], [[0, 0]])
        self.assertFalse(result['count_provenance']['float_log_count_path_used'])
        context = contexts[0]['context']
        digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        self.assertEqual(result['context_sha256'], digest)
        self.assertTrue(all(row['context_sha256'] == digest and not row['complete_selected_weight'] for row in saved))
        self.assertIn('../HD_FIELD128_EXACT_PROFILE_COUNTS.md', result['source_sha256'])

    def test_timeout_preserves_rows_no_selected_weight_claim(self):
        class Budget:
            def __init__(self):
                from time import monotonic
                self.started, self.expired = monotonic(), False
            def check(self):
                if self.expired:
                    raise TimeoutError('test deadline')
        budget, saved = Budget(), []
        def stop(row):
            saved.append(row)
            budget.expired = True
        with self.assertRaises(cert.IncompleteTransfer) as caught:
            cert.transfer_sum(4, 4, as_arb(fixture(4, 4)), budget=budget, checkpoint=stop)
        self.assertEqual(caught.exception.partial['completed_even_rows'], 1)
        self.assertEqual(caught.exception.partial['pending_even_full_range'], [1, 2])
        with patch.object(cert.frozen, 'regional_matrices') as build:
            result = cert.weight_receipt(tiny_bank(), 2, max_seconds=1e-12)
        build.assert_not_called()
        self.assertFalse(result['complete_selected_weight'])
        self.assertIsNone(result['log2_upper_rational'])
        self.assertIsNone(result['context_sha256'])

    def test_strict_scope_proof_pin_and_checkpoint_mutation(self):
        with self.assertRaises(ValueError):
            cert.weight_receipt(object(), 2)
        class Other(cert.provider.AliasBank):
            pass
        with self.assertRaises(ValueError):
            cert.weight_receipt(object.__new__(Other), 2)
        pins = cert.source_pins()
        cert.check_source_pins(pins)
        pins['../HD_FIELD128_EXACT_PROFILE_COUNTS.md'] = 'changed'
        with self.assertRaises(RuntimeError):
            cert.check_source_pins(pins)
        def fail(row):
            raise RuntimeError('mock byte change')
        with patch.object(cert.frozen, 'regional_matrices', return_value=as_arb(fixture(cert.base.PAIRS, 2))):
            with self.assertRaisesRegex(RuntimeError, 'byte change'):
                cert.weight_receipt(tiny_bank(), 2, checkpoint=fail)


if __name__ == '__main__':
    unittest.main()

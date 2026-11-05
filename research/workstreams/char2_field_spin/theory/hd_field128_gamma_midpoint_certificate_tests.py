"""Tiny exact/Arb midpoint sums and mocked receipts, without coefficients."""

from fractions import Fraction as F
from hashlib import sha256
import json
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_gamma_midpoint_certificate as cert
from hd_field128_gamma_near_face_certificate_tests import fixture, as_arb, multiply, power, tiny_bank
from hd_shifted_pairs_count import exact_profiles


class MidpointCertificateTests(unittest.TestCase):
    def test_positive_directed_sum_contains_exact_oracle_and_improves_old(self):
        with ctx.workprec(256):
            for pairs, weight, repetitions in ((4, 4, 8), (5, 6, 8), (3, 3, 256), (2, 0, 8), (2, 4, 8)):
                rows = fixture(pairs, weight)
                counts = cert.midpoint.MidpointCounts(pairs, weight)
                actual = exact_profiles(pairs)
                oracle, true_sum = F(0), F(0)
                for a, left in rows.items():
                    for b, right in rows.items():
                        matrix = power(multiply(left, right), repetitions//2)
                        value = matrix[0]+matrix[1]
                        oracle += counts.count(a, b)*value
                        true_sum += actual.get((weight, a, b), 0)*value
                directed, report = cert.transfer_sum(pairs, weight, as_arb(rows), counts=counts, repetitions=repetitions)
                old, unused = cert.frozen.transfer_sum(pairs, weight, as_arb(rows), repetitions=repetitions, max_singletons=64)
                self.assertTrue(directed >= cert.base.rational(oracle))
                self.assertTrue(directed >= cert.base.rational(true_sum))
                self.assertTrue(directed < cert.base.rational(oracle)*(1+arb(2)**-180))
                self.assertTrue(directed <= old)
                self.assertEqual(report['profile_pairs_evaluated'], len(rows)**2)
                self.assertIsNone(report['pending_even_full_range'])

    def test_zero_regions_and_count_object_reuse(self):
        counts = cert.midpoint.MidpointCounts(4, 4)
        zeros = {f: (arb(0),)*4 for f in counts.full}
        # Guard __init__ rather than replacing the exact producer class.
        with patch.object(cert.midpoint.MidpointCounts, '__init__', side_effect=AssertionError('No rebuild')):
            value, report = cert.transfer_sum(4, 4, zeros, counts=counts)
        self.assertTrue(value.is_zero())
        self.assertEqual(report['completed_even_rows'], 3)
        self.assertTrue(all(row['row_support']['zero'] for row in report['rows']))

    def test_complete_mock_receipt_context_all_profiles_and_both_proof_pins(self):
        bank, saved, contexts = tiny_bank(), [], []
        counts = cert.midpoint.MidpointCounts(cert.base.PAIRS, 2)
        with patch.object(cert.frozen, 'regional_matrices', return_value=as_arb(fixture(cert.base.PAIRS, 2))):
            result = cert.weight_receipt(bank, 2, counts=counts, checkpoint=saved.append, context_checkpoint=contexts.append)
        self.assertEqual(bank.pin_checks, 2)
        self.assertTrue(result['complete_selected_weight'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['profile_pairs_evaluated'], 4)
        self.assertTrue(result['count_provenance']['all_profiles_retained'])
        self.assertFalse(result['count_provenance']['float_log_count_path_used'])
        context = contexts[0]['context']
        digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        self.assertEqual(result['context_sha256'], digest)
        self.assertTrue(all(row['context_sha256'] == digest and not row['complete_selected_weight'] for row in saved))
        for note in ('../HD_FIELD128_MIDPOINT_COUNTS.md', '../HD_FIELD128_EXACT_PROFILE_COUNTS.md'):
            self.assertIn(note, result['source_sha256'])

    def test_timeout_preserves_rows_without_selected_weight_bound(self):
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
        partial = dict(caught.exception.partial)
        partial['pending_even_full_range'] = [1, 1]
        with patch.object(cert.frozen, 'regional_matrices', return_value=as_arb(fixture(cert.base.PAIRS, 2))), \
             patch.object(cert, 'transfer_sum', side_effect=cert.IncompleteTransfer(TimeoutError('stop'), partial)):
            result = cert.weight_receipt(tiny_bank(), 2)
        self.assertEqual(result['completed_even_rows'], 1)
        self.assertIsNone(result['support_sum'])
        self.assertIsNone(result['log2_upper_rational'])

    def test_strict_scope_count_cutoff_pins_and_checkpoint_mutation(self):
        with self.assertRaises(ValueError):
            cert.weight_receipt(object(), 2)
        class Other(cert.provider.AliasBank):
            pass
        with self.assertRaises(ValueError):
            cert.weight_receipt(object.__new__(Other), 2)
        counts = cert.midpoint.MidpointCounts(4, 4, max_singletons=2)
        with self.assertRaises(ValueError):
            cert.transfer_sum(4, 4, as_arb(fixture(4, 4)), counts=counts, max_singletons=3)
        with self.assertRaises(ValueError):
            cert.transfer_sum(4, 3, as_arb(fixture(4, 3)), counts=counts)
        pins = cert.source_pins()
        cert.check_source_pins(pins)
        pins['../HD_FIELD128_MIDPOINT_COUNTS.md'] = 'changed'
        with self.assertRaises(RuntimeError):
            cert.check_source_pins(pins)
        def fail(row):
            raise RuntimeError('mock byte change')
        with patch.object(cert.frozen, 'regional_matrices', return_value=as_arb(fixture(cert.base.PAIRS, 2))):
            with self.assertRaisesRegex(RuntimeError, 'byte change'):
                cert.weight_receipt(tiny_bank(), 2, checkpoint=fail)


if __name__ == '__main__':
    unittest.main()

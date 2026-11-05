"""Tiny exact arithmetic and mocked receipt tests; no numerical proof gates."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_gamma_near_face_certificate as cert
from hd_shifted_pairs_count import exact_profiles


def multiply(a, b):
    return tuple(a[2*i]*b[j]+a[2*i+1]*b[2+j] for i in range(2) for j in range(2))


def power(a, n):
    result = (F(1), F(0), F(0), F(1))
    for unused in range(n):
        result = multiply(result, a)
    return result


def fixture(pairs, weight):
    return {f: (F(1, 256+f), F(2, 257+2*f), F(3, 258+f), F(1, 259+f))
            for f in cert.tail.full_counts(pairs, weight)}


def as_arb(rows):
    return {f: tuple(cert.base.rational(x) for x in row) for f, row in rows.items()}


def tiny_bank():
    # Exact producer class with test-only instance methods; no coefficients built.
    bank = object.__new__(cert.provider.AliasBank)
    bank.precision = 256
    bank.constant, bank.slope = arb(0), arb(0)
    bank.pin_checks = 0
    def check_pins():
        bank.pin_checks += 1
    def metadata():
        return dict(construction='mock-all-pair-bank-no-gate', all_pair_before_alpha=True,
                    state_dimension=30, witness={'mock': True}, directed=True, whole_code_certificate=False)
    bank.check_pins, bank.metadata = check_pins, metadata
    return bank


class NearFaceCertificateTests(unittest.TestCase):
    def test_cutoff_zero_exactly_replays_frozen_transfer(self):
        with ctx.workprec(256):
            for pairs in range(1, 6):
                for weight in range(2*pairs+1):
                    regions = as_arb(fixture(pairs, weight))
                    old, old_rows = cert.frozen.transfer_sum(pairs, weight, regions, repetitions=8)
                    new, new_rows = cert.transfer_sum(pairs, weight, regions, repetitions=8, max_singletons=0)
                    self.assertEqual(new, old)
                    self.assertEqual(new_rows['profile_pairs_evaluated'], old_rows['profile_pairs_evaluated'])

    def test_exact_independent_sum_and_count_dominance(self):
        with ctx.workprec(256):
            for pairs, weight, repetitions in ((3, 3, 8), (4, 4, 8), (4, 5, 8), (3, 4, 256)):
                rows = fixture(pairs, weight)
                counts = cert.near.NearFaceCounts(pairs, weight, max_singletons=pairs)
                actual = exact_profiles(pairs)
                oracle, true_sum = F(0), F(0)
                for first, left in rows.items():
                    for second, right in rows.items():
                        matrix = power(multiply(left, right), repetitions//2)
                        value = matrix[0]+matrix[1]
                        count = counts.count(first, second)
                        true_count = actual.get((weight, first, second), 0)
                        self.assertGreaterEqual(count, true_count)
                        oracle += count*value
                        true_sum += true_count*value
                directed, report = cert.transfer_sum(pairs, weight, as_arb(rows), repetitions=repetitions,
                                                     max_singletons=pairs)
                old, unused = cert.frozen.transfer_sum(pairs, weight, as_arb(rows), repetitions=repetitions)
                self.assertTrue(directed >= cert.base.rational(oracle))
                self.assertTrue(directed >= cert.base.rational(true_sum))
                self.assertTrue(directed < cert.base.rational(oracle)*(1+arb(2)**-180))
                self.assertTrue(directed <= old)
                self.assertIsNone(report['pending_even_full_range'])

    def test_compact_complete_receipt_and_context(self):
        rows, contexts, progress, bank = [], [], [], tiny_bank()
        with patch.object(cert, 'regional_matrices', return_value=as_arb(fixture(cert.base.PAIRS, 2))):
            result = cert.weight_receipt(bank, 2, checkpoint=rows.append,
                                         context_checkpoint=contexts.append, progress=progress.append)
        self.assertEqual(bank.pin_checks, 2)
        self.assertEqual(len(contexts), 1)
        self.assertEqual([row['even_full'] for row in rows], [0, 1])
        self.assertTrue(result['complete_selected_weight'])
        self.assertTrue(result['covers_all_even_odd_splits'])
        self.assertIsNotNone(result['support_sum'])
        self.assertEqual(result['profile_pairs_evaluated'], 4)
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['count_provenance']['float_log_count_path_used'])
        context = contexts[0]['context']
        digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        self.assertEqual(digest, result['context_sha256'])
        for row in rows:
            self.assertEqual(row['context_sha256'], digest)
            self.assertNotIn('witness', row)
            self.assertNotIn('source_sha256', row)
            self.assertFalse(row['complete_selected_weight'])
            self.assertIsNone(row.get('log2_upper_rational'))

    def test_partial_rows_preserved_without_full_bound(self):
        class StopAfterRow:
            def __init__(self):
                from time import monotonic
                self.started, self.expired = monotonic(), False
            def check(self):
                if self.expired:
                    raise TimeoutError('test stop after row')
        budget, rows = StopAfterRow(), []
        def checkpoint(row):
            rows.append(row)
            budget.expired = True
        with self.assertRaises(cert.IncompleteTransfer) as caught:
            cert.transfer_sum(4, 4, as_arb(fixture(4, 4)), budget=budget, checkpoint=checkpoint)
        partial = caught.exception.partial
        self.assertEqual(len(rows), 1)
        self.assertEqual(partial['completed_even_rows'], 1)
        self.assertEqual(partial['pending_even_full_range'], [1, 2])
        public_partial = dict(partial, pending_even_full_range=[1, 1], profile_pairs_evaluated=2)
        bank = tiny_bank()
        with patch.object(cert, 'regional_matrices', return_value=as_arb(fixture(cert.base.PAIRS, 2))), \
             patch.object(cert, 'transfer_sum', side_effect=cert.IncompleteTransfer(TimeoutError('stop'), public_partial)):
            result = cert.weight_receipt(bank, 2)
        self.assertEqual(bank.pin_checks, 2)
        self.assertEqual(result['status'], 'incomplete')
        self.assertIsNone(result['support_sum'])
        self.assertIsNone(result['log2_upper_rational'])
        self.assertFalse(result['covers_all_supports_at_weight'])

    def test_strict_scope_source_checks_timeout_and_immutable_writes(self):
        with self.assertRaises(ValueError):
            cert.weight_receipt(object(), 2)
        with self.assertRaises(ValueError):
            cert.weight_receipt(object.__new__(cert.frozen.provider.AliasBank), 2)
        class OtherProducer(cert.provider.AliasBank):
            pass
        with self.assertRaises(ValueError):
            cert.weight_receipt(object.__new__(OtherProducer), 2)
        for kwargs in ({'max_singletons': -1}, {'max_singletons': True}, {'failure_bits': -1}):
            with self.assertRaises(ValueError):
                cert.weight_receipt(tiny_bank(), 2, **kwargs)
        with patch.object(cert, 'regional_matrices') as build:
            result = cert.weight_receipt(tiny_bank(), 2, max_seconds=1e-12)
        build.assert_not_called()
        self.assertFalse(result['complete_selected_weight'])
        self.assertEqual(result['completed_even_rows'], 0)
        pins = cert.source_pins()
        for name in ('hd_field128_gamma_near_face_certificate.py', 'hd_field128_gamma_near_face_counts.py',
                     'hd_field128_gamma_all_pair_certificate.py', 'hd_field128_gamma_profile_certificate.py',
                     'hd_field128_lift_alias_cover.py'):
            self.assertIn(name, pins)
        cert.check_source_pins(pins)
        bad = dict(pins, **{'hd_field128_gamma_near_face_certificate.py': 'changed'})
        with self.assertRaises(RuntimeError):
            cert.check_source_pins(bad)
        def mutate(row):
            raise RuntimeError('mock pinned-source mutation')
        with patch.object(cert, 'regional_matrices', return_value=as_arb(fixture(cert.base.PAIRS, 2))):
            with self.assertRaisesRegex(RuntimeError, 'mutation'):
                cert.weight_receipt(tiny_bank(), 2, checkpoint=mutate)
        with TemporaryDirectory() as directory:
            path = Path(directory)/'row.json'
            cert.write_new(path, {'complete_selected_weight': False})
            with self.assertRaises(FileExistsError):
                cert.write_new(path, {})


if __name__ == '__main__':
    unittest.main()

"""Tiny scalar-oracle and mocked receipt tests; no construction gates."""

from fractions import Fraction as F
from hashlib import sha256
from contextlib import redirect_stdout
from io import StringIO
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_gamma_profile_certificate as cert


def multiply(left, right):
    return tuple(left[2*i]*right[j]+left[2*i+1]*right[2+j] for i in range(2) for j in range(2))


def power(matrix, exponent):
    result = (F(1), F(0), F(0), F(1))
    while exponent:
        if exponent & 1:
            result = multiply(result, matrix)
        exponent >>= 1
        if exponent:
            matrix = multiply(matrix, matrix)
    return result


def exact_regions(pairs, weight, q=256):
    # Small-field-sized scalar fixtures, not a claim about a gamma construction over GF256.
    return {f: (F(1, q+f), F(2, q+2*f+1), F(3, q+f+2), F(1, q+f+3))
            for f in cert.tail.full_counts(pairs, weight)}


def arb_regions(regions):
    return {f: tuple(cert.base.rational(value) for value in entries) for f, entries in regions.items()}


class TinyBank(cert.provider.AliasBank):
    """No coefficients or FFTs: only receipt-control metadata is mocked."""

    def __init__(self):
        self.precision = 256
        self.constant, self.slope = arb(0), arb(0)

    def check_pins(self):
        pass

    def metadata(self):
        return dict(construction='mock-gamma-bank-no-construction-gate', directed=True,
                    whole_code_certificate=False, state_dimension=24, witness={'mock': True})


class ProfileCertificateTests(unittest.TestCase):
    def test_cached_count_is_exact_original_bound_on_all_small_profiles(self):
        for pairs in range(1, 12):
            for weight in range(2*pairs+1):
                counts = cert.ProfileCounts(pairs, weight, cert.base.Budget())
                for first in counts.full:
                    self.assertEqual(counts.marginals[first], cert.tail.one_pairing_count(pairs, weight, first))
                    for second in counts.full:
                        self.assertEqual(counts.count(first, second),
                                         cert.profile_count_upper(pairs, weight, first, second))

    def test_small_q_scalar_oracle_and_fixed_power128(self):
        with ctx.workprec(256):
            for pairs, weight, repetitions in ((2, 1, 2), (3, 3, 8), (3, 4, 256), (2, 4, 256)):
                exact = exact_regions(pairs, weight)
                oracle = F(0)
                for first, left in exact.items():
                    for second, right in exact.items():
                        count = cert.profile_count_upper(pairs, weight, first, second)
                        value = power(multiply(left, right), repetitions//2)
                        oracle += count*(value[0]+value[1])
                result, report = cert.transfer_sum(pairs, weight, arb_regions(exact), repetitions=repetitions)
                self.assertTrue(result >= cert.base.rational(oracle))
                self.assertTrue(result < cert.base.rational(oracle)*(1+arb(2)**-180))
                self.assertEqual(report['profile_pairs_evaluated'], len(exact)**2)
                self.assertIsNone(report['pending_even_full_range'])
            left, right = tuple(map(arb, (2, 3, 5, 7))), tuple(map(arb, (11, 13, 17, 19)))
            fixed = cert.pair_row_sum128(left, right)
            generic = cert.base.matrix_power(cert.base.matrix_product(left, right), 128)
            oracle = power(multiply(tuple(map(F, (2, 3, 5, 7))), tuple(map(F, (11, 13, 17, 19)))), 128)
            exact = cert.base.rational(oracle[0]+oracle[1])
            self.assertTrue(fixed >= exact)
            self.assertTrue(generic[0]+generic[1] >= exact)
            self.assertTrue(fixed < exact*(1+arb(2)**-180))

    def test_transfer_is_below_existing_common_potential_run_bound(self):
        with ctx.workprec(256):
            for weight in (1, 2, 3, 4, 7, 8):
                regions = arb_regions(exact_regions(4, weight))
                value, unused = cert.transfer_sum(4, weight, regions, repetitions=8)
                for ratio in (F(1, 4), F(1), F(4)):
                    common = cert.tail.support_sum(4, weight, regions, ratio, method='runs', repetitions=8)
                    self.assertTrue(value <= common)

    def test_complete_weight_receipt_and_immediate_row_checkpoints(self):
        rows, progress = [], []
        weight = 2
        regions = arb_regions(exact_regions(cert.base.PAIRS, weight))
        with patch.object(cert, 'regional_matrices', return_value=regions):
            result = cert.weight_receipt(TinyBank(), weight, checkpoint=rows.append, progress=progress.append)
        self.assertEqual([row['even_full'] for row in rows], [0, 1])
        self.assertTrue(result['complete_selected_weight'])
        self.assertTrue(result['covers_all_supports_at_weight'])
        self.assertEqual(result['profile_pairs_evaluated'], 4)
        self.assertIsNotNone(result['log2_upper_rational'])
        self.assertFalse(result['common_potential_used'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(progress[-1]['completed_even_rows'], 2)
        for row in rows:
            self.assertFalse(row['complete_selected_weight'])
            self.assertFalse(row['covers_all_supports_at_weight'])
            self.assertEqual(row['source_sha256'], result['source_sha256'])
            self.assertIsNone(row.get('log2_upper_rational'))

    def test_timeout_after_row_keeps_checkpoints_but_never_returns_full_bound(self):
        rows = []
        original = cert.pair_row_sum128
        def expire_after_first(left, right):
            if rows:
                raise TimeoutError('mock stop after checkpoint')
            return original(left, right)
        regions = arb_regions(exact_regions(cert.base.PAIRS, 2))
        with patch.object(cert, 'regional_matrices', return_value=regions), \
             patch.object(cert, 'pair_row_sum128', side_effect=expire_after_first):
            result = cert.weight_receipt(TinyBank(), 2, checkpoint=rows.append)
        self.assertEqual(len(rows), 1)
        self.assertEqual(result['completed_even_rows'], 1)
        self.assertEqual(result['pending_even_full_range'], [1, 1])
        self.assertEqual(result['status'], 'incomplete')
        self.assertIsNone(result['log2_upper_rational'])
        self.assertIsNone(result['log2_upper_enclosure'])
        self.assertFalse(result['covers_all_even_odd_splits'])
        self.assertFalse(result['meets_requested_selected_weight_target'])
        self.assertEqual(result['failure']['error_type'], 'TimeoutError')

    def test_validation_timeout_before_regions_and_immutable_pins(self):
        with self.assertRaises(ValueError):
            cert.weight_receipt(object(), 2)
        for weight in (0, 8193, True):
            with self.assertRaises(ValueError):
                cert.weight_receipt(TinyBank(), weight)
        with self.assertRaises(ValueError):
            cert.transfer_sum(2, 2, {0: (arb(1),)*4})
        with patch.object(cert, 'regional_matrices') as build:
            result = cert.weight_receipt(TinyBank(), 2, max_seconds=1e-12)
        build.assert_not_called()
        self.assertIsNone(result['log2_upper_rational'])
        self.assertEqual(result['completed_even_rows'], 0)
        with TemporaryDirectory() as directory:
            path = Path(directory)/'even_full_0000.json'
            receipt = cert.write_new(path, {'mock': True})
            self.assertEqual(receipt['sha256'], sha256(path.read_bytes()).hexdigest())
            self.assertEqual(json.loads(path.read_bytes()), {'mock': True})
            with self.assertRaises(FileExistsError):
                cert.write_new(path, {})
        pins = cert.source_pins()
        for name in ('hd_field128_gamma_profile_certificate.py', 'hd_field128_full_rs_gamma_alias_certificate.py',
                     'hd_field128_lift_alias_cover.py', 'hd_shifted_pairs_count.py'):
            self.assertIn(name, pins)
        cert.check_source_pins(pins)
        changed = dict(pins)
        changed['hd_field128_gamma_profile_certificate.py'] = 'changed'
        with self.assertRaises(RuntimeError):
            cert.check_source_pins(changed)

    def test_incomplete_cli_writes_summary_and_exits_nonzero(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)/'new-output'
            argv = ['profile-certificate', '--r', '2', '--z', '1/2', '--alpha', '1/8',
                    '--marker', '1/2', '1/3', '--output-dir', str(output)]
            partial = dict(status='incomplete', complete_selected_weight=False,
                           completed_even_rows=1, pending_even_full_range=[1, 1],
                           log2_upper_rational=None, failure={'error_type': 'TimeoutError'})
            bank = TinyBank()
            with patch('sys.argv', argv), patch.object(cert.provider, 'AliasBank', return_value=bank), \
                 patch.object(cert, 'weight_receipt', return_value=partial), redirect_stdout(StringIO()):
                with self.assertRaises(SystemExit) as error:
                    cert.main()
            self.assertEqual(error.exception.code, 2)
            result = json.loads((output/'summary.json').read_bytes())
            self.assertFalse(result['complete_selected_weight'])
            self.assertIsNone(result['log2_upper_rational'])


if __name__ == '__main__':
    unittest.main()

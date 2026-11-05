"""Tiny count and mock-bank artifact audits; no coefficient construction."""

from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import ctx

import hd_field128_gamma_near_face_band_audit as audit
import hd_field128_gamma_near_face_counts as near
from hd_field128_gamma_near_face_band_tests import FakeCache, read_link


def synthetic_bank():
    bank = object.__new__(audit.band.provider.AliasBank)
    bank.precision = 256
    bank.construction = audit.EXPECTED['construction']
    bank.check_pins = lambda: None
    metadata = dict(audit.EXPECTED, bank_source_sha256=audit.band.provider.source_pins(),
        witness=dict(z='1/2', alpha='1/8', markers=[['1/3', '1/5']], precision=256,
                     packet_mantissa_bits=80, grid=128, state_dimension=30, inner_width=64,
                     gamma_extension_degree=16, gamma_base_field_order=256,
                     joint_rank_return_cap=True, all_pair_before_alpha=True))
    with ctx.workprec(256):
        bank.constant = audit.outer_prefix(metadata, 0)
        bank.slope = audit.outer_prefix(metadata, 1)-bank.constant
    bank.metadata = lambda: deepcopy(metadata)
    return bank


def fixture(directory, first=1, last=3):
    FakeCache.fail_weight = None
    with patch.object(audit.band, 'RegionalCache', FakeCache):
        result = audit.band.band_receipt(synthetic_bank(), first, last, output_dir=directory,
                                        max_weights=last-first+1, max_singletons=2)
    return result


class AuditTests(unittest.TestCase):
    def test_independent_counts_all_small_profiles(self):
        for pairs in range(1, 8):
            for weight in range(1, 2*pairs+1):
                for cutoff in (0, 1, 2, 2*pairs):
                    independent = audit.IndependentCounts(pairs, weight, cutoff)
                    reference = near.NearFaceCounts(pairs, weight, max_singletons=cutoff)
                    for a in independent.full:
                        for b in independent.full:
                            self.assertEqual(independent.count(a, b), reference.count(a, b))
        self.assertEqual(audit._runs(0, 0), [1])
        self.assertEqual(audit._runs(4, 0), [1])
        self.assertEqual(audit._runs(4, 4), [1, 0, 0, 0, 0])

    def test_complete_resolution_fresh_prefix_and_union(self):
        with TemporaryDirectory() as root:
            directory = Path(root)/'band'
            result = fixture(directory)
            checked = audit.audit_band(directory/'summary.json')
            self.assertEqual(checked['audit_result'], 'PASS')
            self.assertEqual(checked['weight_range'], [1, 3])
            self.assertEqual(checked['saved_band_log2_upper_rational'], result['log2_upper_rational'])
            self.assertFalse(checked['coefficient_replay'])
            self.assertFalse(checked['whole_code_certificate'])
            self.assertEqual([item['profile_pairs'] for item in checked['weight_audits']], [1, 4, 4])
            self.assertEqual(len(checked['checked_file_pins']), 5)
            self.assertIn(Path(audit.__file__).name, checked['source_sha256'])
            self.assertTrue(all(item['fresh_outer_prefix'] for item in checked['weight_audits']))

    def test_rows_reject_coverage_counts_prefix_and_understatement(self):
        with TemporaryDirectory() as root:
            result = fixture(Path(root)/'band', 2, 2)
            context, row = read_link(result['context']), read_link(result['weight_receipts'][0])
            def reject(change):
                bad = deepcopy(row)
                change(bad)
                with ctx.workprec(384), self.assertRaises((ValueError, KeyError)):
                    audit.audit_rows(bad, context)
            reject(lambda x: x['rows'].pop())
            reject(lambda x: x['rows'][0].update(profiles_with_positive_count=0))
            reject(lambda x: x.update(log2_outer_prefix_enclosure='0'))
            reject(lambda x: x.update(log2_upper_rational='-1000000000'))
            reject(lambda x: x['rows'][0]['completed_support'].update(log2_upper_rational='-1000000000'))
            reject(lambda x: x.update(profile_pairs_evaluated=3))
            with patch.object(audit.band.selected.near.NearFaceCounts, 'count', return_value=0), \
                 ctx.workprec(384), self.assertRaisesRegex(ValueError, 'producer exact profile counts differ'):
                audit.audit_rows(row, context)

    def test_resolver_rejects_changed_bytes_partial_and_missing_weight(self):
        with TemporaryDirectory() as root:
            directory = Path(root)/'band'
            result = fixture(directory)
            bad = deepcopy(result)
            bad['weight_receipts'].pop()
            incomplete = audit.band.write_new(Path(root)/'missing.json', bad)
            with self.assertRaisesRegex(ValueError, 'weight links'):
                audit.audit_band(incomplete['path'])
            bad = deepcopy(result)
            bad['complete_requested_range'] = False
            partial = audit.band.write_new(Path(root)/'partial.json', bad)
            with self.assertRaises(ValueError):
                audit.audit_band(partial['path'])
            Path(result['context']['path']).write_bytes(b'changed test-only bytes')
            with self.assertRaises(ValueError):
                audit.audit_band(directory/'summary.json')

    def test_context_scope_and_source_guards(self):
        with TemporaryDirectory() as root:
            result = fixture(Path(root)/'band', 1, 1)
            context = read_link(result['context'])
            for key, value in [('state_dimension', 31), ('target_distance', '9/20'),
                               ('bank_source_sha256', {}), ('source_sha256', {})]:
                bad = deepcopy(context)
                bad[key] = value
                with self.assertRaises(ValueError):
                    audit._verify_witness(bad)
            bad = deepcopy(context)
            bad['witness']['markers'] = [['0', '1']]
            with self.assertRaisesRegex(ValueError, 'marker'):
                audit._verify_witness(bad)

    def test_promotion_checks_audit_summary_binding_before_prior_events(self):
        with TemporaryDirectory() as root:
            root = Path(root)
            fixture(root/'band', 1, 1)
            checked = audit.audit_band(root/'band'/'summary.json')
            checked['saved_band_log2_upper_rational'] = '-1000'
            pin = audit.band.write_new(root/'bad-audit.json', checked)
            parent = dict(schema='field128-gamma-partial-cover-index-v9',
                          **{key: audit.EXPECTED[key] for key in ('construction', 'field_order',
                             'source_dimension', 'output_length', 'target_distance', 'whole_code_certificate')})
            parent_pin = audit.band.write_new(root/'parent.json', parent)
            with self.assertRaisesRegex(ValueError, 'audit/summary binding'):
                audit.extend_v9(parent_pin['path'], checked, pin)


if __name__ == '__main__':
    unittest.main()

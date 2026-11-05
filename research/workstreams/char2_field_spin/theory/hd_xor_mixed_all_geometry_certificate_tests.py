"""Exact assignment, scope and one-time count-transfer checks; no bank replay."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from flint import arb, ctx

from hd_shifted_coverage import rational_ceiling
from hd_xor_mixed_all_geometry_certificate import (ALL_MASKS, EXPECTED_SCOPE,
    PROPOSAL_SHA256, PRIMITIVES_SHA256, checked_classes, checked_rows,
    receipt, rounded_union_512, verify)


def fake_receipt():
    entries = [dict(active_row_mask=mask, active_constituents=mask.bit_count(),
        witness_label=label, mark=str(z), log2_potential=str(logp),
        common_pair_failure_log2_upper='-200', failure_log2_upper='-127')
        for label, z, logp, masks in checked_classes() for mask in masks]
    return dict(EXPECTED_SCOPE, precision_bits=256, frozen_proposal_sha256=PROPOSAL_SHA256,
        reviewed_primitives_sha256=PRIMITIVES_SHA256,
        source_sha256=sha256(Path(__file__).with_name('hd_xor_mixed_all_geometry_certificate.py').read_bytes()).hexdigest(),
        entries=entries)


class AllGeometryCertificateTests(unittest.TestCase):
    def test_exact_five_bank219_mask_assignment(self):
        classes = checked_classes()
        masks = [mask for _, _, _, items in classes for mask in items]
        self.assertEqual(sorted(masks), list(ALL_MASKS))
        self.assertEqual(len(set(masks)), 219)
        self.assertEqual([len(row[3]) for row in classes], [6, 68, 104, 35, 6])
        self.assertEqual(next((z, p) for label, z, p, masks in classes if 195 in masks),
                         (F('4487137477776021/4503599627370496'), F('-70552129/1048576')))
        self.assertEqual(next((z, p) for label, z, p, masks in classes if 224 in masks),
                         (F('8975046825465429/9007199254740992'), F('-2381845503/33554432')))

    def test_scope_rejected_before_regional_work(self):
        for key, value in (('cohort_counts', [1, 15, 1]), ('source_matching_counts', [0, 8]),
                           ('mu_range', [16, 16]), ('sparse_pair_assignments', 28),
                           ('exact_mass_multiplier', 28**16), ('field_order', 255),
                           ('sparse_support_floor', 543), ('dense_nonempty_column_floor', 16),
                           ('pair_swaps', 1), ('whole_code_certificate', True)):
            with patch('hd_xor_mixed_all_geometry_certificate.arb_pair_regions') as regional:
                with self.assertRaises(ValueError):
                    receipt(scope=dict(EXPECTED_SCOPE, **{key: value}))
                regional.assert_not_called()

    def test_receipt_partition_literals_arity_and_count_fields_are_checked(self):
        valid = fake_receipt()
        self.assertEqual(set(checked_rows(valid, 256)), set(ALL_MASKS))
        variants = []
        for key, value in (('source_matching_counts', [0, 8]), ('base_sparse_pair_assignments', 28**16),
                           ('exact_mass_multiplier', 28**16), ('fixed_row_sets', 218)):
            bad = deepcopy(valid)
            bad[key] = value
            variants.append(bad)
        for key, value in (('active_constituents', 8), ('mark', '99/100'),
                           ('witness_label', 'setup-selected'), ('log2_potential', '-60')):
            bad = deepcopy(valid)
            bad['entries'][0][key] = value
            variants.append(bad)
        bad = deepcopy(valid)
        bad['entries'][0] = bad['entries'][1]
        variants.append(bad)
        for bad in variants:
            with self.assertRaises(ValueError):
                checked_rows(bad, 256)

    def test_exact_positive_union_has_one_geometry_factor_and_no_bank_penalty(self):
        common, transferred = rounded_union_512([F(0), F(1), F(3)])
        with ctx.workprec(1024):
            self.assertEqual(common, rational_ceiling(arb(11).log()/arb(2).log()))
            self.assertEqual(transferred, rational_ceiling(arb(11*28**15).log()/arb(2).log()))
            with_extra_bank_factor = rational_ceiling(arb(5*11*28**15).log()/arb(2).log())
        self.assertGreater(with_extra_bank_factor-transferred, 2)

    def test_observed_dual_replays_and178_unchanged_cells(self):
        result = verify()
        self.assertEqual(result['rounded_cell_union_log2_upper'], '-61000916136179/1000000000000')
        self.assertEqual([row['sha256'] for row in result['inputs']], [
            '071a50922931379a0878c78a9c0ea594ca62ccca5e1e728123fa505316ec6642',
            '81da27a2ddcc6807aa5c908d299630fb1c5c812bdaafb236f2228943bc071df4'])
        for precision in (256, 384):
            directory = Path(__file__).parent
            old = json.loads((directory/f'hd_xor_mixed_mask_8_0_mu16_receipt{precision}.json').read_text())
            new = json.loads((directory/f'hd_xor_mixed_all_geometry_8_0_receipt{precision}.json').read_text())
            previous = {row['active_row_mask']: row for row in old['entries']}
            retained = [row for row in new['entries'] if not row['witness_label'].startswith('all_geometry_')]
            self.assertEqual(len(retained), 178)
            for row in retained:
                original = previous[row['active_row_mask']]
                self.assertEqual((row['mark'], row['log2_potential']),
                                 (original['mark'], original['log2_potential']))
                self.assertEqual(row['common_pair_failure_log2_upper'], original['failure_log2_upper'])


if __name__ == '__main__':
    unittest.main()

"""Tiny synthetic receipt tests; no regional replay or benchmark."""

from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

import hd_field128_lift_partial_index as index


def alias_receipt(first, bounds):
    rows = [dict(active_groups=weight, regional_histograms=weight // 2 - max(0, weight - index.PAIRS) + 1,
                 log2_upper_rational=str(bound), log2_upper_enclosure=str(F(bound).__floor__() - 1),
                 state_potential_ratio='1', covers_all_supports_at_weight=True, covers_all_even_odd_splits=True)
            for weight, bound in enumerate(bounds, first)]
    last = first + len(rows) - 1
    method = next(iter(index.ALIAS_METHODS))
    return dict(index.Validator(F(7, 20)).full(), coefficient_method=method,
                source_sha256=index.source_closure(index.ALIAS_METHODS[method]),
                coverage_driver_sha256=index.source_closure(index.COVER), support_method='runs',
                witness=dict(z='1/2', alpha='1/16', precision=256, packet_mantissa_bits=80),
                requested_weight_range=[first, last], weight_receipts=rows, weights_evaluated=len(rows),
                every_integer_weight_evaluated=True, complete_requested_range=True,
                completed_weight_range=[first, last], pending_weight_ranges=[],
                completed_union_log2_upper_rational='0')


def terminal_receipt(bound=-200):
    return dict(index.Validator(F(7, 20)).common(), weight_range=index.TERMINAL_RANGE,
                source_dimension=index.DIMENSION, output_length=index.OUTPUT_LENGTH,
                projective_divisor=str(index.FIELD_ORDER - 1), covers_all_supports_in_range=True,
                covers_all_even_odd_splits=True, no_support_subset_union=True,
                method='inverse-image-terminal-range', source_sha256=index.source_closure(index.TAIL),
                witness=dict(mark='1/2', precision=256), log2_upper_rational=str(bound),
                log2_upper_enclosure=str(bound - 1))


class PartialIndexTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix='field128-index-test-')
        self.addCleanup(self.directory.cleanup)
        self.base = Path(self.directory.name)

    def save(self, name, value):
        path = self.base / name
        path.write_text(json.dumps(value), encoding='utf-8')
        return path

    def test_missing_gap_not_interpolated(self):
        a = self.save('a.json', alias_receipt(1, [-100]))
        b = self.save('b.json', alias_receipt(3, [-110]))
        result = index.build_index([a, b])
        self.assertEqual(result['passing_weight_ranges'], [[1, 1], [3, 3]])
        self.assertIn(2, result['missing_weights'])
        self.assertFalse(result['all_weights_covered'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['union_event_count'], 2)

    def test_overlap_selects_exact_minimum_once(self):
        small = -F(100) - F(1, 10**30)
        a = self.save('a.json', alias_receipt(1, [-100, -90]))
        b = self.save('b.json', alias_receipt(1, [small]))
        result = index.build_index([a, b, a])
        row = result['selected_per_weight'][0]
        self.assertEqual(F(row['log2_upper_rational']), small)
        self.assertEqual(row['compatible_candidate_count'], 2)
        self.assertEqual(result['union_event_count'], 2)
        self.assertEqual(len(result['explicit_inputs']), 2)

    def test_terminal_is_one_event_and_replaces_overlap(self):
        a = self.save('a.json', alias_receipt(5854, [-500]))
        t = self.save('t.json', terminal_receipt(-200))
        u = self.save('u.json', terminal_receipt(-210))
        result = index.build_index([a, t, u])
        self.assertEqual(result['union_event_count'], 1)
        self.assertEqual(result['terminal_union_event_count'], 1)
        self.assertEqual(result['passing_weight_count'], 8192 - 5854 + 1)
        self.assertFalse(result['selected_per_weight'][0]['used_in_union'])
        self.assertEqual(F(result['terminal_event']['log2_upper_rational']), -210)
        self.assertLessEqual(F(result['selected_union_log2_upper_rational']), F(-209999999999999, 10**12))

    def test_above_target_not_missing_and_not_summed(self):
        path = self.save('a.json', alias_receipt(1, [-79, -81]))
        result = index.build_index([path])
        self.assertEqual(result['known_above_target_weights'], [1])
        self.assertNotIn(1, result['missing_weights'])
        self.assertEqual(result['union_event_count'], 1)

    def test_reject_wrong_field_distance_dimensions_and_direction(self):
        for key, value in [('field_order', 256), ('target_distance', '49/100'),
                           ('source_dimension', 42), ('directed', False)]:
            with self.subTest(key=key):
                receipt = alias_receipt(1, [-100])
                receipt[key] = value
                path = self.save(key + '.json', receipt)
                with self.assertRaisesRegex(ValueError, 'incompatible ' + key):
                    index.build_index([path])

    def test_reject_wrong_source_pin(self):
        receipt = alias_receipt(1, [-100])
        name = next(iter(receipt['source_sha256']))
        receipt['source_sha256'][name] = '0' * 64
        path = self.save('bad.json', receipt)
        with self.assertRaisesRegex(ValueError, 'source pins'):
            index.build_index([path])

    def test_reject_ceiling_below_enclosure(self):
        receipt = alias_receipt(1, [-100])
        receipt['weight_receipts'][0]['log2_upper_enclosure'] = '-99'
        path = self.save('bad.json', receipt)
        with self.assertRaisesRegex(ValueError, 'ceiling'):
            index.build_index([path])

    def test_batch_pins_child_and_inherits_checked_dimensions(self):
        child = alias_receipt(257, [-100, -110])
        child_path = self.save('child.json', child)
        pin = dict(path=str(child_path), sha256=sha256(child_path.read_bytes()).hexdigest())
        rows = [dict(row, band_id='one', receipt_sha256=pin['sha256']) for row in child['weight_receipts']]
        receipt = dict(index.Validator(F(7, 20)).common(), schema='field128-lift-alias-batch-v1',
                       pair_swaps=True, full_field_random_coefficients=True, old_numerical_ceilings_reused=False,
                       source_sha256=index.source_closure(index.BATCH), requested_weight_range=[257, 259],
                       weight_receipts=rows, covered_weight_count=2, duplicate_weight_count=0,
                       pending_weight_ranges=[[259, 259]], complete_requested_range=False,
                       completed_bands=[dict(id='one', weight_range=[257, 258], receipt=pin,
                                             complete_requested_range=True, completed_union_log2_upper_rational='0')])
        path = self.save('batch.json', receipt)
        result = index.build_index([path])
        self.assertEqual(result['passing_weight_ranges'], [[257, 258]])
        self.assertEqual(len(result['receipt_file_pins']), 2)
        altered = deepcopy(receipt)
        altered['weight_receipts'][0]['log2_upper_rational'] = '-1000'
        with self.assertRaisesRegex(ValueError, 'differ'):
            index.build_index([self.save('altered.json', altered)])
        child['source_dimension'] = 42
        self.save('child.json', child)
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            index.build_index([path])

    def test_wide_explicit_disjoint_ranges_keep_gap(self):
        witness = dict(z='1/2', alpha='1/160', precision=256, packet_mantissa_bits=80, paired=True)
        rows = [dict(active_groups=r, profiles=1, log2_upper_rational='-100', log2_upper_enclosure='-101',
                     covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
                     directed=True, whole_code_certificate=False,
                     witness=dict(witness, coefficient_maximum=128, coefficient_marker='1/32')) for r in (65, 67)]
        band = dict(index.Validator(F(7, 20)).full(), source_sha256=index.source_closure(index.WIDE),
                    coefficient_method='positive-carry-free-Arb-truncated-bivariate-cutoff256',
                    covers_all_supports_at_completed_weights=True, covers_all_even_odd_splits=True,
                    requested_weight_ranges=[[65, 65], [67, 67]], pending_weight_ranges=[],
                    weights_evaluated=2, every_requested_weight_evaluated=True, complete_requested_ranges=True,
                    coefficient_maximum=128, coefficient_marker='1/32', witness=witness,
                    weight_receipts=rows, old_witness_job=0)
        job = dict(old_numerical_ceilings_reused=False, requested_weight_ranges=band['requested_weight_ranges'],
                   coefficient_maximum=128, z=witness['z'], alpha=witness['alpha'])
        receipt = dict(index.Validator(F(7, 20)).full(), source_sha256=index.source_closure(index.WIDE),
                       plan=dict(old_numerical_ceilings_reused=False, directed=False, whole_code_certificate=False,
                                 requested_weight_ranges=band['requested_weight_ranges'], jobs=[job]),
                       band_receipts=[band], weight_receipts=[dict(row, job_index=0) for row in rows],
                       requested_weight_ranges=band['requested_weight_ranges'], completed_weight_count=2,
                       duplicate_weight_count=0, complete_requested_ranges=True, pending_weight_ranges=[])
        result = index.build_index([self.save('wide.json', receipt)])
        self.assertEqual(result['passing_weight_ranges'], [[65, 65], [67, 67]])
        self.assertIn(66, result['unresolved_weights'])

    def test_sparse_rows_match_pinned_original_or_repair(self):
        witness = dict(z='1/2', alpha='1/16', precision=256, packet_mantissa_bits=80,
                       paired=True, coefficient_maximum=1, coefficient_marker='1/4096')
        row = dict(active_groups=1, profiles=1, log2_upper_rational='-100', log2_upper_enclosure='-101',
                   covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
                   directed=True, whole_code_certificate=False, witness=witness)
        band = dict(index.Validator(F(7, 20)).full(), source_sha256=index.source_closure(index.BANDS),
                    weight_receipts=[row])
        original = dict(index.Validator(F(7, 20)).common(), source_sha256=index.source_closure(index.BANDS),
                        completed_bands=[band])
        original_path = self.save('original.json', original)
        pin = dict(path=str(original_path), sha256=sha256(original_path.read_bytes()).hexdigest())
        selected = dict(row, origin=dict(kind='input-band-row', band_index=0, row_index=0))
        receipt = dict(index.Validator(F(7, 20)).full(), source_sha256=index.source_closure(index.SPARSE),
                       repair_source_sha256=index.source_closure(index.CERTIFICATE),
                       selected_weight_receipts=[selected], requested_weight_range=[1, 1],
                       selected_weight_count=1, duplicate_weight_count=0, complete_requested_range=True,
                       every_weight_evaluated=True, covers_all_supports_in_range=True, covers_all_even_odd_splits=True,
                       input_provenance=pin, input_receipt=original, completed_repairs=[])
        result = index.build_index([self.save('repaired.json', receipt)])
        self.assertEqual(result['passing_weight_ranges'], [[1, 1]])
        self.assertEqual(len(result['receipt_file_pins']), 2)
        repair = dict(row, **{key: value for key, value in index.Validator(F(7, 20)).full().items() if key not in row},
                      source_sha256=index.source_closure(index.CERTIFICATE))
        repair['log2_upper_rational'], repair['log2_upper_enclosure'] = '-120', '-121'
        receipt['completed_repairs'] = [repair]
        receipt['selected_weight_receipts'] = [dict(repair, origin=dict(kind='own-weight-repair', repair_index=0))]
        result = index.build_index([self.save('own-repair.json', receipt)])
        self.assertEqual(result['selected_per_weight'][0]['log2_upper_rational'], '-120')
        receipt['selected_weight_receipts'][0]['log2_upper_rational'] = '-999'
        with self.assertRaisesRegex(ValueError, 'differs'):
            index.build_index([self.save('bad-repair.json', receipt)])

    def test_empty_and_unknown_inputs_rejected(self):
        with self.assertRaises(ValueError):
            index.build_index([])
        with self.assertRaisesRegex(ValueError, 'Unsupported'):
            index.build_index([self.save('unknown.json', {})])


if __name__ == '__main__':
    unittest.main()

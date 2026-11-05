"""Synthetic pinned-file contract tests; no coefficients or real index writes."""

from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_gamma_selected_index_v16 as mode


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    mode.band.write_new(path, value)
    return dict(path=str(path.resolve()), sha256=sha256(path.read_bytes()).hexdigest())


def parent_fixture(folder, *, source_conflict=False):
    event = write(folder/'prior.json', dict(synthetic_fixture=True))
    parent = dict(schema='field128-gamma-partial-cover-index-v15',
        construction=mode.midpoint.previous.EXPECTED['construction'], field_order=1 << 128,
        source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
        whole_code_certificate=False, directed_aggregation=True,
        entries=[dict(event, weight_range=[1, 1], log2_upper_rational='-100')],
        covered_weight_ranges=[[1, 1]], covered_weight_count=1,
        pending_weight_ranges=[[2, 8192]], pending_weight_count=8191, receipt_count=1,
        covered_events_union_log2_upper_rational='-99', producer_receipt_counts={},
        aggregation_source_sha256=({'hd_field128_gamma_midpoint_weight_bands.py': 'wrong'}
                                   if source_conflict else {}))
    pin = write(folder/'parent.json', parent)
    return Path(pin['path']), pin['sha256']


def band_fixture(folder, first=2, last=3, *, zero=False, mutate=None):
    weights = list(range(first, last+1))
    roots = mode._roots(first, last)
    descriptors = [root.descriptor() for root in roots]
    witness = dict(z='1/2', alpha='1/8', markers=[['1/2', '1/3']], precision=256,
                   grid=128, packet_mantissa_bits=80, state_dimension=30, inner_width=64,
                   gamma_extension_degree=16, gamma_base_field_order=256,
                   joint_rank_return_cap=True, all_pair_before_alpha=True)
    context = dict(mode.midpoint.previous.EXPECTED, directed=False,
        schema='gamma-midpoint-weight-band-context-v1', pairs=4096,
        repetitions=256, paired_region_power=128, requested_weight_range=[first, last],
        requested_weights=weights, source_sha256=mode.producer.source_pins(),
        bank_source_sha256=mode.producer.provider.source_pins(),
        normalized_every_original_weight=True, count_pointwise_minimum_used=False,
        count_only_complement_above_pairs=True, no_profile_pair_table=True,
        no_weight_interpolation=True, complete_requested_band=False,
        parity_weights={str(eps): [r for r in weights if r % 2 == eps] for eps in (0, 1)
                        if any(r % 2 == eps for r in weights)},
        witness=witness, witness_source=None, root_cells=descriptors)
    digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    nodes = []
    for i, root in enumerate(roots):
        # These are synthetic cell numerical receipts. No claim is made about
        # executing a production coefficient bank for this test fixture.
        receipt = dict(cell=root.descriptor(), context_sha256=digest, complete_cell=True,
            status='complete', count_upper=1, zero=zero,
            log2_upper_rational=None if zero else '-200', log2_upper_enclosure=None if zero else '-200',
            failure_bits=64, meets_cell_budget=True)
        nodes.append(dict(id=i, parent=None, cell=root.descriptor(), status='accepted', receipt=receipt))
    with ctx.workprec(384):
        value = None if zero else mode.producer.tail.log2arb(arb(len(roots))*arb(2)**-200)
        fresh = None if zero else mode.producer.base.rational_ceiling(value)
        saved = None if zero else fresh+F(1, 10**12)
    summary = dict(schema='gamma-midpoint-weight-band-cover-v1', status='complete',
        complete_requested_band=True, covers_all_supports_in_band=True, covers_all_matching_profiles=True,
        every_integer_weight_prepared=True, directed=True, numerical_backend='actual-all-pair-bank',
        no_weight_interpolation=True, whole_code_certificate=False, pending_frontier=[], stop_reason=None,
        batch_source_verified_at_entry=True, batch_source_verified_at_exit=True,
        requested_weight_range=[first, last], requested_weights=weights, context=context,
        context_sha256=digest, root_cells=descriptors, cell_failure_bits=64,
        accepted_leaves=list(range(len(nodes))), tree_nodes=len(nodes), evaluated_nodes=len(nodes),
        nodes=nodes, zero=zero, log2_upper_rational=None if zero else str(saved),
        log2_upper_enclosure=None if zero else str(value), requested_failure_bits=54,
        meets_requested_band_target=True)
    if mutate:
        mutate(summary)
    summary_pin = write(folder/'summary.json', summary)
    report = dict(schema='gamma-midpoint-weight-band-independent-audit-v1', audit_result='PASS',
        input_summary=summary_pin, complete_requested_band=True, covers_all_matching_profiles=True,
        integer_partition_checked=True, counts_recomputed=True, cell_log_enclosures_checked=True,
        directed_union_precision=384, coefficient_replay_in_this_audit=False,
        normalized_regional_replay=False, source_sha256=mode.audited.source_pins(),
        checked_source_sha256=mode.audited.source_pins(),
        whole_code_certificate=False, requested_weight_range=[first, last], requested_weights=weights,
        context_sha256=digest, checked_file_pins=[summary_pin],
        log2_upper_rational=None if zero else str(saved),
        saved_band_log2_upper_rational=None if zero else str(saved),
        reaggregated_cell_log2_upper_rational=None if zero else str(fresh))
    pin = write(folder/'audit.json', report)
    return Path(pin['path']), summary_pin


class IndexTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.folder = Path(self.temporary.name)
        self.parent, self.parent_hash = parent_fixture(self.folder/'old')

    def tearDown(self):
        self.temporary.cleanup()

    def extend(self, paths):
        with patch.object(mode, 'PARENT_SHA256', self.parent_hash):
            return mode.extend_v15(self.parent, paths)

    def test_two_disjoint_bands_counted_once_and_parent_unchanged(self):
        first, unused = band_fixture(self.folder/'first', 2, 3)
        second, unused = band_fixture(self.folder/'second', 5, 6)
        before = self.parent.read_bytes()
        result = self.extend([first, second])
        self.assertEqual(result['schema'], 'field128-gamma-partial-cover-index-v16')
        self.assertEqual(result['receipt_count'], 3)
        self.assertEqual(result['added_event_count'], 2)
        self.assertEqual(result['added_weight_count'], 4)
        self.assertEqual(result['covered_weight_ranges'], [[1, 3], [5, 6]])
        self.assertEqual(result['pending_weight_ranges'], [[4, 4], [7, 8192]])
        expected = mode.index.union_ceiling([F(row['log2_upper_rational']) for row in result['entries']], 384)
        self.assertEqual(result['covered_events_union_log2_upper_rational'], str(expected))
        self.assertFalse(result['whole_code_certificate'])
        self.assertTrue(result['new_integer_partitions_rechecked'])
        self.assertEqual(self.parent.read_bytes(), before)

    def test_zero_event_covers_weights_but_adds_no_positive_union(self):
        path, unused = band_fixture(self.folder/'zero', 2, 3, zero=True)
        result = self.extend([path])
        self.assertTrue(result['added_events_zero'])
        self.assertIsNone(result['added_events_union_log2_upper_rational'])
        self.assertEqual(result['covered_weight_count'], 3)
        self.assertTrue(result['entries'][1]['zero'])
        self.assertIsNone(result['entries'][1]['log2_upper_rational'])
        self.assertEqual(result['covered_events_union_log2_upper_rational'],
                         str(mode.index.union_ceiling([F(-100)], 384)))

    def test_overlaps_duplicates_and_explicit_input_guards(self):
        old_overlap, unused = band_fixture(self.folder/'overlap', 1, 2)
        fresh, unused = band_fixture(self.folder/'fresh', 2, 3)
        for paths in ([old_overlap], [fresh, fresh], [], [fresh]*33):
            with self.assertRaises(ValueError):
                self.extend(paths)
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            mode.extend_v15(self.parent, [fresh])

    def test_bad_geometry_partial_pin_scope_or_ceiling_rejected(self):
        changes = [
            lambda s: s.update(complete_requested_band=False),
            lambda s: s.update(requested_weights=[2]),
            lambda s: s['nodes'].pop(),
            lambda s: s['nodes'][0]['cell']['u_deficit_ranges'][0].__setitem__(0, 0),
            lambda s: s.update(context_sha256='wrong'),
            lambda s: s['context'].update(state_dimension=31),
            lambda s: s['context']['source_sha256'].update({'bad.py': 'wrong'}),
            lambda s: s.update(log2_upper_rational='-100000'),
            lambda s: s['nodes'][0]['receipt'].update(log2_upper_rational='-100000'),
            lambda s: s.update(directed=False),
        ]
        for i, change in enumerate(changes):
            path, unused = band_fixture(self.folder/f'bad{i}', mutate=change)
            with self.assertRaises((ValueError, KeyError)):
                self.extend([path])

    def test_audit_summary_hash_and_source_conflict_fail_closed(self):
        path, summary_pin = band_fixture(self.folder/'valid')
        read = Path.read_bytes
        target = Path(summary_pin['path'])
        def altered(p):
            raw = read(p)
            return raw+b' ' if p == target else raw
        with patch.object(Path, 'read_bytes', altered):
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                self.extend([path])
        self.parent, self.parent_hash = parent_fixture(self.folder/'conflict', source_conflict=True)
        with self.assertRaisesRegex(ValueError, 'Conflicting source pin'):
            self.extend([path])


if __name__ == '__main__':
    unittest.main()

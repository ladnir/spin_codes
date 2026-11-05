"""Checkpoint export guards; all expensive arithmetic replays are mocked."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import unittest
from unittest.mock import call, patch, sentinel

import hd_shifted_checkpoint_export as exporter
from hd_shifted_cover_ledger import Box
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION
from hd_shifted_swap_ledger import SwapLedger


def mock_receipt(first, last, *, bound='-100', endpoint=False):
    ranges = ((first, last), (0, 0), (0, 0)) if endpoint else (
        (first, last), (0, 4096), (0, 4096))
    return dict(ranges=ranges, log2_upper_rational=bound, directed=True,
                covers_all_even_odd_splits=True, construction=SWAP_CONSTRUCTION,
                field_order=256, forbidden_weight=734003,
                whole_code_certificate=False, test_fresh_receipt=True)


def mock_tail_source(*, terminal=True, first_bound='-101'):
    source = SwapLedger()
    source.accept(mock_receipt(4096, 4096, endpoint=True))
    # The first leaf crosses the requested boundary. Its whole probability
    # bound must survive restriction to 5120, without a volume discount.
    source.accept(mock_receipt(5119, 6143, bound=first_bound))
    source.accept(mock_receipt(6144, 7167, bound='-102'))
    source.accept(mock_receipt(7168, 7648, bound='-103'))
    if terminal:
        source.accept(mock_receipt(7649, 8192, bound='-104'))
    return source


class CheckpointExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = Path(exporter.__file__).parent
        cls.raw_plans = {name: (cls.directory/name).read_bytes()
                         for _, name, _, _, _ in exporter.PLANS}

    def test_established_plan_pins_and_assigned_ranges(self):
        self.assertEqual([(kind, low, high) for kind, _, low, high, _ in exporter.PLANS],
                         [('unsigned', 5120, 6143), ('unsigned', 6144, 7167),
                          ('unsigned', 7168, 7648)])
        for _, name, _, _, expected in exporter.PLANS:
            self.assertEqual(sha256(self.raw_plans[name]).hexdigest(), expected)

    def test_any_changed_pin_fails_before_any_replay(self):
        for _, changed, _, _, _ in exporter.PLANS:
            with self.subTest(changed=changed):
                def read_bytes(path):
                    raw = self.raw_plans[path.name]
                    return raw+b'\n' if path.name == changed else raw
                with patch.object(Path, 'read_bytes', read_bytes), \
                        patch.object(exporter, 'replay_swap_cover') as replay:
                    with self.assertRaisesRegex(ValueError, 'Established proposal changed'):
                        exporter.tail_checkpoint(self.directory)
                    replay.assert_not_called()

    def test_missing_plan_fails_before_replay(self):
        with patch.object(Path, 'read_bytes', side_effect=FileNotFoundError('missing')), \
                patch.object(exporter, 'replay_swap_cover') as replay:
            with self.assertRaises(FileNotFoundError):
                exporter.tail_checkpoint(self.directory)
            replay.assert_not_called()

    def test_sparse_replays_exactly_one_through_256_with_explicit_caps(self):
        def family(last, **kwargs):
            return mock_receipt(1, last)
        def batch(first, last, z, alpha, **kwargs):
            return mock_receipt(first, last)
        with patch.object(exporter, 'family_receipt', side_effect=family) as family_call, \
                patch.object(exporter, 'batch_receipt', side_effect=batch) as batch_call:
            result = exporter.sparse_checkpoint(384, sentinel.progress)
        family_call.assert_called_once_with(64, precision=384, progress=sentinel.progress)
        expected_calls = [call(low, high, z, alpha, maximum=maximum, precision=384,
                               progress=sentinel.progress)
                          for maximum, schedule in ((128, exporter.BATCH_SCHEDULE),
                                                     (256, exporter.BATCH_SCHEDULE_256))
                          for low, high, z, alpha in schedule]
        self.assertEqual(batch_call.call_args_list, expected_calls)
        weights = [weight for receipt in result['accepted']
                   for weight in range(receipt['ranges'][0][0], receipt['ranges'][0][1]+1)]
        self.assertEqual(weights, list(range(1, 257)))
        self.assertEqual(result['summary']['root'], ((1, 256), (0, 4096), (0, 4096)))
        self.assertTrue(result['summary']['complete_local_cover'])
        self.assertEqual(result['summary']['unresolved_profile_triples'], 0)
        self.assertTrue(result['fresh_arithmetic'])
        self.assertEqual(result['construction'], SWAP_CONSTRUCTION)
        self.assertEqual(result['precision'], 384)
        self.assertEqual(result['kind'], 'sparse')
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['summary']['whole_code_certificate'])
        self.assertTrue(all(item['test_fresh_receipt'] for item in result['accepted']))

    def test_sparse_gap_is_not_exported_as_complete(self):
        with patch.object(exporter, 'family_receipt', return_value=mock_receipt(1, 63)), \
                patch.object(exporter, 'batch_receipt',
                             side_effect=lambda first, last, *a, **k: mock_receipt(first, last)):
            with self.assertRaisesRegex(ArithmeticError, 'sparse checkpoint did not replay'):
                exporter.sparse_checkpoint()

    def test_sparse_weak_union_or_nondirected_receipt_fails_closed(self):
        with patch.object(exporter, 'family_receipt', return_value=mock_receipt(1, 64, bound='-39')), \
                patch.object(exporter, 'batch_receipt',
                             side_effect=lambda first, last, *a, **k: mock_receipt(first, last)):
            with self.assertRaises(ArithmeticError):
                exporter.sparse_checkpoint()
        bad = mock_receipt(1, 64)
        bad['directed'] = False
        with patch.object(exporter, 'family_receipt', return_value=bad), \
                patch.object(exporter, 'batch_receipt') as batch:
            with self.assertRaises(ValueError):
                exporter.sparse_checkpoint()
            batch.assert_not_called()

    def test_tail_exports_only_5120_through_8192_with_fresh_scope_flags(self):
        source = mock_tail_source()
        source.accepted[1].update(source_label='original-before-tail-import',
                                  witness=dict(markers=['1/3', '2/5'], precision=384))
        source_snapshot = deepcopy(source.accepted)
        provenance = [dict(label='established-plan',
                           replay=dict(witnesses=[dict(z='1/7', alpha='1/10')],
                                       accepted=[deepcopy(source.accepted[1])]))]
        provenance_snapshot = deepcopy(provenance)
        replay_result = dict(sources=provenance,
                             summary=dict(whole_code_certificate=True, complete_local_cover=True,
                                          accepted_union_log2_upper='-100000'))
        with patch.object(exporter, 'replay_swap_cover',
                          return_value=(source, replay_result)) as replay:
            result = exporter.tail_checkpoint(self.directory, 384, sentinel.progress)
        expected_inputs = [(name, json.loads(self.raw_plans[name]), low, high)
                           for _, name, low, high, _ in exporter.PLANS]
        replay.assert_called_once_with(expected_inputs, precision=384, threshold=-40,
                                       progress=sentinel.progress)
        self.assertEqual(result['summary']['root'], ((5120, 8192), (0, 4096), (0, 4096)))
        self.assertEqual([item['ranges'][0] for item in result['accepted']],
                         [(5120, 6143), (6144, 7167), (7168, 7648), (7649, 8192)])
        self.assertEqual(result['accepted'][0]['restricted_from_ranges'][0], (5119, 6143))
        self.assertEqual(result['accepted'][0]['log2_upper_rational'], '-101')
        self.assertEqual(result['accepted'][0]['source_label'], 'established-unsigned-tail')
        self.assertEqual(result['sources'], provenance_snapshot)
        self.assertEqual(result['source_accepted'], source_snapshot)
        self.assertEqual(result['source_accepted'][1]['source_label'], 'original-before-tail-import')
        self.assertEqual(result['source_accepted'][1]['witness'],
                         dict(markers=['1/3', '2/5'], precision=384))
        self.assertEqual(source.accepted[0]['ranges'][0], (4096, 4096))
        self.assertEqual(source.accepted[1]['ranges'][0], (5119, 6143))
        self.assertEqual(result['proposal_sha256'], {name: expected
                         for _, name, _, _, expected in exporter.PLANS})
        self.assertTrue(result['fresh_arithmetic'])
        self.assertTrue(result['summary']['complete_local_cover'])
        self.assertEqual(result['summary']['unresolved_profile_triples'], 0)
        self.assertLess(F(result['summary']['accepted_union_log2_upper']), -40)
        self.assertNotEqual(result['summary']['accepted_union_log2_upper'], '-100000')
        self.assertEqual(result['construction'], SWAP_CONSTRUCTION)
        self.assertEqual(result['precision'], 384)
        self.assertEqual(result['kind'], 'tail')
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['summary']['whole_code_certificate'])

    def test_tail_gap_or_weak_union_fails_closed(self):
        for source in (mock_tail_source(terminal=False), mock_tail_source(first_bound='-39')):
            with self.subTest(terminal=len(source.accepted) == 5):
                with patch.object(exporter, 'replay_swap_cover',
                                  return_value=(source, dict(sources=[]))):
                    with self.assertRaisesRegex(ArithmeticError, 'dense tail did not replay'):
                        exporter.tail_checkpoint(self.directory)


if __name__ == '__main__':
    unittest.main()

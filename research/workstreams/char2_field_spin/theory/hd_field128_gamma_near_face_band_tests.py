"""Tiny mocked-bank band tests; no FFT construction or numerical proof gates."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import arb

import hd_field128_gamma_near_face_band as band


def tiny_bank():
    bank = object.__new__(band.provider.AliasBank)
    bank.precision = 256
    bank.constant, bank.slope = arb(-11), arb(2)
    bank.construction = 'test-only-actual-class-fixture-no-gate'
    bank.check_pins = lambda: None
    bank.metadata = lambda: dict(construction=bank.construction, state_dimension=30,
                                witness={'test_only': True}, directed=True, whole_code_certificate=False)
    return bank


class FakeCache:
    created = 0
    visited = []
    fail_weight = None

    def __init__(self, bank):
        type(self).created += 1

    def regional(self, weight, budget):
        type(self).visited.append(weight)
        budget.check()
        if weight == type(self).fail_weight:
            raise TimeoutError('test normalization timeout')
        return {f: tuple(band.base.rational(F(n, 257+f)) for n in (1, 2, 3, 1))
                for f in band.tail.full_counts(band.base.PAIRS, weight)}


def read_link(link):
    raw = Path(link['path']).read_bytes()
    if sha256(raw).hexdigest() != link['sha256']:
        raise AssertionError('Artifact link hash mismatch')
    return json.loads(raw)


class BandTests(unittest.TestCase):
    def setUp(self):
        FakeCache.created, FakeCache.visited, FakeCache.fail_weight = 0, [], None

    def test_every_integer_weight_one_cache_and_compact_linked_receipts(self):
        with TemporaryDirectory() as parent, patch.object(band, 'RegionalCache', FakeCache):
            output, bank = Path(parent)/'band', tiny_bank()
            result = band.band_receipt(bank, 1, 3, output_dir=output, max_weights=3, max_singletons=0)
            self.assertEqual(FakeCache.created, 1)
            self.assertEqual(FakeCache.visited, [1, 2, 3])
            self.assertTrue(result['every_integer_weight_evaluated'])
            self.assertTrue(result['complete_requested_range'])
            self.assertEqual(result['pending_weight_ranges'], [])
            self.assertIsNone(result['partial_weight_receipt'])
            context = read_link(result['context'])
            self.assertEqual(context['count_provenance']['maximum_refined_singletons'], 0)
            self.assertIn('hd_field128_gamma_near_face_band.py', context['source_sha256'])
            rows = [read_link(link) for link in result['weight_receipts']]
            for weight, row in zip(range(1, 4), rows):
                self.assertEqual(row['context'], result['context'])
                self.assertEqual(row['active_groups'], weight)
                self.assertEqual(row['completed_even_rows'], weight//2+1)
                self.assertNotIn('witness', row)
                self.assertNotIn('source_sha256', row)
                self.assertTrue(row['complete_selected_weight'])
                self.assertFalse(row['whole_code_certificate'])
                self.assertEqual(arb(row['log2_outer_prefix_enclosure']), arb(-11+2*weight))
            union = band.cover.union_ceiling(rows, bank.precision)
            self.assertEqual(result['log2_upper_rational'], union)
            self.assertEqual(json.loads((output/'summary.json').read_bytes()), result)

    def test_weight_timeout_preserves_completed_union_and_pending_integer_range(self):
        FakeCache.fail_weight = 2
        with TemporaryDirectory() as parent, patch.object(band, 'RegionalCache', FakeCache):
            result = band.band_receipt(tiny_bank(), 1, 3, output_dir=Path(parent)/'band', max_weights=3)
            self.assertEqual(FakeCache.visited, [1, 2])
            self.assertEqual(result['weights_evaluated'], 1)
            self.assertEqual(result['pending_weight_ranges'], [[2, 3]])
            self.assertIsNone(result['log2_upper_rational'])
            self.assertIsNotNone(result['completed_union_log2_upper_rational'])
            partial = read_link(result['partial_weight_receipt'])
            self.assertIsNone(partial['support_sum'])
            self.assertIsNone(partial['log2_upper_rational'])
            self.assertEqual(partial['completed_even_rows'], 0)
            self.assertFalse(result['requested_range_meets_target'])

    def test_timeout_inside_profile_rows_keeps_partial_weight_artifact(self):
        def stop(row):
            if row['stage'] == 'completed-even-profile-row':
                raise TimeoutError('test stop after first even row')
        with TemporaryDirectory() as parent, patch.object(band, 'RegionalCache', FakeCache):
            result = band.band_receipt(tiny_bank(), 4, 5, output_dir=Path(parent)/'band', progress=stop)
            self.assertEqual(result['weights_evaluated'], 0)
            self.assertEqual(result['pending_weight_ranges'], [[4, 5]])
            self.assertIsNone(result['completed_union_log2_upper_rational'])
            partial = read_link(result['partial_weight_receipt'])
            self.assertEqual(partial['completed_even_rows'], 1)
            self.assertEqual(partial['pending_even_full_range'], [1, 2])
            self.assertEqual(len(partial['rows']), 1)
            self.assertFalse(partial['complete_selected_weight'])

    def test_scope_budget_and_exclusive_paths(self):
        with TemporaryDirectory() as parent:
            output = Path(parent)/'band'
            for kwargs in ({'max_weights': 1}, {'failure_bits': -1}, {'max_singletons': -1}):
                with self.assertRaises(ValueError):
                    band.band_receipt(tiny_bank(), 1, 2, output_dir=output, **kwargs)
            self.assertFalse(output.exists())
            with self.assertRaises(ValueError):
                band.band_receipt(object(), 1, 2, output_dir=output)
            class Subclass(band.provider.AliasBank):
                pass
            with self.assertRaises(ValueError):
                band.band_receipt(object.__new__(Subclass), 1, 2, output_dir=output)
            with patch.object(band, 'RegionalCache', FakeCache):
                result = band.band_receipt(tiny_bank(), 1, 2, output_dir=output, max_seconds=1e-12)
            self.assertEqual(FakeCache.visited, [])
            self.assertEqual(result['pending_weight_ranges'], [[1, 2]])
            self.assertFalse(result['complete_requested_range'])
            with self.assertRaises(FileExistsError):
                band.band_receipt(tiny_bank(), 1, 2, output_dir=output)

    def test_source_mutation_prevents_publishing_weight_or_summary(self):
        original = band.selected.check_source_pins
        changed = False
        def verify(pins):
            if changed:
                raise RuntimeError('mock pinned-source mutation')
            return original(pins)
        def mutate(row):
            nonlocal changed
            if row['stage'] == 'completed-even-profile-row':
                changed = True
        with TemporaryDirectory() as parent, patch.object(band, 'RegionalCache', FakeCache), \
             patch.object(band.selected, 'check_source_pins', side_effect=verify):
            output = Path(parent)/'band'
            with self.assertRaisesRegex(RuntimeError, 'mutation'):
                band.band_receipt(tiny_bank(), 1, 1, output_dir=output, progress=mutate)
            self.assertTrue((output/'context.json').exists())
            self.assertFalse((output/'weight_0001.json').exists())
            self.assertFalse((output/'summary.json').exists())

    def test_union_target_is_not_a_per_weight_target_and_links_are_verified(self):
        def fake_weight(bank, cache, weight, **kwargs):
            return dict(active_groups=weight, status='complete', complete_selected_weight=True,
                        log2_upper_rational='-64', completed_even_rows=1,
                        pending_even_full_range=None, failure=None)
        with TemporaryDirectory() as parent, patch.object(band, 'RegionalCache', FakeCache), \
             patch.object(band, '_weight', side_effect=fake_weight):
            output = Path(parent)/'band'
            result = band.band_receipt(tiny_bank(), 1, 2, output_dir=output, failure_bits=64)
            self.assertTrue(result['complete_requested_range'])
            self.assertEqual(result['weights_above_target'], [])
            self.assertFalse(result['requested_range_meets_target'])
            self.assertGreaterEqual(F(result['log2_upper_rational']), F(-63))
            Path(result['context']['path']).write_bytes(b'changed test-only context')
            with self.assertRaisesRegex(RuntimeError, 'artifact changed'):
                band._check_links((result['context'],))


if __name__ == '__main__':
    unittest.main()

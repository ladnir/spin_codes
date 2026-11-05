"""Tiny synthetic-atlas tests of the real bank interface; no distance gates."""

from dataclasses import FrozenInstanceError
from fractions import Fraction as F
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, ctx

import hd_field128_gamma_marker_subset as subset


def tiny_bank():
    # Use the strict public bank constructor with a tiny synthetic positive
    # coefficient fixture. These tests do not claim a construction receipt.
    composed = {(0, 0): tuple(map(arb, (2, 1, 1, 3))),
                (1, 0): tuple(map(arb, (1, 2, 3, 1))),
                (0, 1): tuple(map(arb, (3, 1, 2, 2)))}
    with patch.object(subset.provider, 'directed_parts', return_value=({}, composed)):
        return subset.provider.AliasBank(F(1, 2), F(1, 8),
                    ((F(1, 2), F(1, 3)), (F(2, 3), F(1, 2)), (F(3, 4), F(2, 3))),
                    state=30, grid=64, max_seconds=30)


class MarkerSubsetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parent = tiny_bank()
        cls.view = subset.MarkerSubsetBank(cls.parent, (2, 0))

    def test_caps_are_not_smaller_and_full_subset_replays(self):
        with ctx.workprec(self.parent.precision):
            all_cache = subset.cover.RegionalCache(self.parent)
            part_cache = subset.cover.RegionalCache(self.view)
            full_cache = subset.cover.RegionalCache(subset.MarkerSubsetBank(self.parent, (0, 1, 2)))
            for weight in range(1, 7):
                all_rows = all_cache.regional(weight, subset.cover.base.Budget())
                part_rows = part_cache.regional(weight, subset.cover.base.Budget())
                full_rows = full_cache.regional(weight, subset.cover.base.Budget())
                self.assertEqual(full_rows, all_rows)
                for key, row in all_rows.items():
                    self.assertTrue(all(a >= b for a, b in zip(part_rows[key], row)))

    def test_immutable_owned_selection_and_parent_unchanged(self):
        before = subset._parent_signature(self.parent)
        self.assertEqual(self.view.indices, (2, 0))
        self.assertEqual(self.view.markers, (self.parent.markers[2], self.parent.markers[0]))
        self.assertNotIsInstance(self.view, subset.provider.AliasBank)
        for atlas, index in zip(self.view.atlases, self.view.indices):
            self.assertFalse(np.shares_memory(atlas.upper, self.parent.atlases[index].upper))
            self.assertTrue(np.array_equal(atlas.upper, self.parent.atlases[index].upper))
            with self.assertRaises(ValueError):
                atlas.upper.setflags(write=True)
        with self.assertRaises(FrozenInstanceError):
            self.view.indices = (1,)
        with self.assertRaises(TypeError):
            self.view.pins['changed'] = 'bad'
        self.assertEqual(before, subset._parent_signature(self.parent))

    def test_exact_dynamic_metadata_and_complete_lineage(self):
        original = subset.provider.AliasBank.metadata(self.parent)
        result = self.view.metadata()
        self.assertEqual(result['witness']['markers'],
                         [original['witness']['markers'][index] for index in (2, 0)])
        self.assertEqual(result['parent_witness'], original['witness'])
        self.assertEqual(result['regional_diagnostics'],
                         [original['regional_diagnostics'][index] for index in (2, 0)])
        self.assertEqual(result['construction'], original['construction'])
        self.assertEqual(result['parent_bank_source_sha256'], self.parent.pins)
        self.assertEqual((result['subset_marker_count'], result['parent_marker_count']), (2, 3))
        self.assertTrue(result['all_pair_before_alpha'])
        self.assertFalse(result['whole_code_certificate'])
        for name in ('hd_field128_gamma_marker_subset.py', 'hd_field128_gamma_all_pair_certificate.py',
                     'hd_field128_lift_alias_cover.py'):
            self.assertIn(name, result['bank_source_sha256'])
        result['witness']['markers'].clear()
        result['parent_witness']['markers'].clear()
        result['regional_diagnostics'].clear()
        self.assertEqual(subset.provider.AliasBank.metadata(self.parent), original)
        self.assertEqual(len(self.view.metadata()['witness']['markers']), 2)

    def test_invalid_parent_and_indices(self):
        class Derived(subset.provider.AliasBank):
            pass

        for parent in (SimpleNamespace(), self.view, object(), object.__new__(Derived)):
            with self.assertRaises(TypeError):
                subset.MarkerSubsetBank(parent, (0,))
        for indices in ((), (0, 0), (-1,), (3,), (True,), (1.0,)):
            with self.assertRaises(ValueError):
                subset.MarkerSubsetBank(self.parent, indices)

    def test_parent_value_source_and_selected_data_mutations_rejected(self):
        for name, changed in (('markers', tuple(reversed(self.parent.markers))),
                              ('constant', arb(-1)), ('precision', 128)):
            with patch.object(self.parent, name, changed):
                with self.assertRaises(RuntimeError):
                    self.view.check_pins()
        # A mutation of an unselected fallback is also detected.
        altered = list(self.parent.numerators)
        altered[1] = tuple(map(arb, (1, 1, 1, 1)))
        with patch.object(self.parent, 'numerators', altered):
            with self.assertRaises(RuntimeError):
                self.view.check_pins()
        altered_atlases = list(self.parent.atlases)
        atlas = altered_atlases[1]
        data = atlas.upper.copy()
        data[0, 0, 0] += 1
        data.setflags(write=False)
        altered_atlases[1] = subset.provider.interval.RegionalAliases(data, atlas.factors)
        with patch.object(self.parent, 'atlases', altered_atlases):
            with self.assertRaises(RuntimeError):
                self.view.check_pins()
        with patch.object(subset.provider, 'source_pins', return_value={}):
            with self.assertRaises(RuntimeError):
                self.view.check_pins()
        # Deliberate bypass of frozen dataclass guards is detected at exit.
        bad = subset.MarkerSubsetBank(self.parent, (0,))
        object.__setattr__(bad, 'constant', arb(0))
        with self.assertRaises(RuntimeError):
            bad.check_pins()
        bad = subset.MarkerSubsetBank(self.parent, (0,))
        object.__setattr__(bad, 'indices', (1,))
        with self.assertRaisesRegex(RuntimeError, 'selection'):
            bad.check_pins()
        self.view.check_pins()

    def test_generic_driver_and_before_after_wrapper(self):
        before = subset._parent_signature(self.parent)
        result = subset.band_receipt(self.view, 1, 2, (F(1),),
                                     max_weights=2, max_seconds=30)
        self.assertTrue(result['complete_requested_range'])
        self.assertEqual(result['subset_marker_indices'], [2, 0])
        self.assertEqual(len(result['witness']['markers']), 2)
        self.assertTrue(result['subset_source_verified_at_entry'])
        self.assertTrue(result['subset_source_verified_at_exit'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(before, subset._parent_signature(self.parent))
        original = self.parent.constant

        def mutate_then_fail(*args, **kwargs):
            self.parent.constant = arb(0)
            raise TimeoutError('synthetic interruption')

        try:
            with patch.object(subset.cover, 'band_receipt', side_effect=mutate_then_fail):
                with self.assertRaisesRegex(RuntimeError, 'Parent bank'):
                    subset.band_receipt(self.view, 1, 1, (F(1),))
        finally:
            self.parent.constant = original
        self.view.check_pins()


if __name__ == '__main__':
    unittest.main()

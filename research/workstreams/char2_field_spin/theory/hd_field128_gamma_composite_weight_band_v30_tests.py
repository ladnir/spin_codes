"""Tiny exact domains and synthetic generic covers; no actual bank or files."""
from copy import deepcopy
from fractions import Fraction as F
import gc
from hashlib import sha256
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import weakref

from flint import arb, ctx
import hd_field128_gamma_composite_weight_band_v30 as composite
import hd_field128_gamma_portable_band_runner_v25 as controller

bands, base = controller.bands, controller.base


def table(pairs=4, first=3, last=5, index=0, *, constant=None, slope=None, repetitions=4, calls=None):
    def regional(r, budget):
        if calls is not None:
            calls.append(r)
        return {f: tuple(base.rational(F(x, 16)) for x in
            ((2+r, 1, 2+f, 3) if index == 0 else (3, 2+f, 1, 2+r)))
            for f in bands.tail.full_counts(pairs, r)}
    return bands.BandTables(pairs, first, last, regional,
        arb(-40 if index == 0 else -35) if constant is None else arb(constant),
        arb(2 if index == 0 else -1) if slope is None else arb(slope),
        repetitions=repetitions, budget=base.Budget())


def cells_in(tables):
    for r in tables.weights:
        u, parity = divmod(r, 2)
        limit = bands.deficit_limit(tables.pairs, r)
        for d in range(limit+1):
            for e in range(limit+1):
                yield bands.Rectangle(parity, ((u, u), (d, d), (e, e)))


class GenericEvaluator:
    def __init__(self, tables):
        self.tables, self._active = tables, False
    def check_pins(self):
        pass
    def context(self):
        return dict(requested_weight_range=[self.tables.weights[0], self.tables.weights[-1]],
                    requested_weights=list(self.tables.weights), synthetic_test=True)


class CompositeBandTests(unittest.TestCase):
    def test_complete_cell_minimum_own_constants_slopes_and_original_weights(self):
        calls = [[], []]
        parts = [table(index=i, calls=calls[i]) for i in range(2)]
        tables = composite.CompositeTables(parts, bands)
        self.assertEqual(calls, [[3, 4, 5], [3, 4, 5]])
        with ctx.workprec(256):
            for cell in list(tables.roots)+list(cells_in(tables)):
                separate = [t.receipt(cell, 48, base.Budget()) for t in parts]
                result = tables.receipt(cell, 48, base.Budget())
                zero = any(r['zero'] for r in separate)
                self.assertEqual(result['zero'], zero)
                self.assertEqual(result['count_upper'], separate[0]['count_upper'])
                if not zero:
                    self.assertEqual(F(result['log2_upper_rational']), min(F(r['log2_upper_rational']) for r in separate))
                self.assertEqual(result['component_bounds'][result['component_winner']]['log2_upper_rational'],
                                 result['log2_upper_rational'])
                self.assertFalse(result['entrywise_cross_witness_minimum'])

    def test_entrywise_splicing_would_be_strictly_wrong(self):
        def regional(row):
            return lambda r, budget: {f: tuple(map(arb, row)) for f in bands.tail.full_counts(2, r)}
        with ctx.workprec(256):
            parts = [bands.BandTables(2, 2, 2, regional(row), arb(c), arb(0), repetitions=4)
                     for row, c in (((1, 0, 0, 4), 0), ((4, 0, 0, 1), -6))]
            tables = composite.CompositeTables(parts, bands)
            result = tables.receipt(tables.roots[0], 40, base.Budget())
            mixed = bands.BandTables(2, 2, 2, regional((1, 0, 0, 1)), arb(-6), arb(0), repetitions=4)
            wrong = mixed.receipt(mixed.roots[0], 40, base.Budget())
            self.assertGreater(F(result['log2_upper_rational']), F(wrong['log2_upper_rational']))
            self.assertEqual(result['component_winner'], 0)

    def test_same_event_and_count_domains_required(self):
        parts = [table(index=i) for i in range(2)]
        with self.assertRaises(ValueError):
            composite.CompositeTables([parts[0], table(first=4)], bands)
        tables = composite.CompositeTables(parts, bands)
        actual = parts[1].count
        with patch.object(parts[1], 'count', new=lambda cell, budget: actual(cell, budget)+1):
            with self.assertRaises(ValueError):
                tables.receipt(tables.roots[0], 54, base.Budget())
        with patch.object(parts[1], 'receipt', side_effect=TimeoutError('unfinished witness')):
            with self.assertRaises(TimeoutError):
                tables.receipt(tables.roots[0], 54, base.Budget())

    def test_upper_half_tail_domain_and_every_weight_exactly_once(self):
        parts = [table(pairs=4, first=5, last=7, index=i, constant=-100, slope=0) for i in range(2)]
        tables = composite.CompositeTables(parts, bands)
        seen = set()
        for cell in cells_in(tables):
            u, d, e = (interval[0] for interval in cell.ranges)
            r = 2*u+cell.parity
            a, b = u-d, u-e
            self.assertIn(a, bands.tail.full_counts(4, r))
            self.assertIn(b, bands.tail.full_counts(4, r))
            self.assertNotIn((r, a, b), seen)
            seen.add((r, a, b))
            root = tables.root_by_parity[cell.parity]
            self.assertTrue(all(lo <= x[0] <= hi for x, (lo, hi) in zip(cell.ranges, root.ranges)))
        expected = {(r, a, b) for r in (5, 6, 7) for a in bands.tail.full_counts(4, r)
                    for b in bands.tail.full_counts(4, r)}
        self.assertEqual(seen, expected)
        evaluator = GenericEvaluator(tables)
        result = controller.cover(evaluator, max_nodes=32, cell_failure_bits=54, aggregate_bits=48)
        self.assertTrue(result['complete_requested_band'])
        self.assertEqual(result['requested_weights'], [5, 6, 7])
        self.assertFalse(result['directed'])
        self.assertEqual(result['numerical_backend'], 'synthetic-test-only')
        with self.assertRaises(ValueError):
            composite.production_envelope(evaluator, result, controller)

    def test_generic_partial_stays_partial_no_relabel(self):
        tables = composite.CompositeTables([table(index=i, constant=100, slope=0) for i in range(2)], bands)
        result = controller.cover(GenericEvaluator(tables), max_nodes=1, cell_failure_bits=54, aggregate_bits=48)
        self.assertFalse(result['complete_requested_band'])
        self.assertFalse(result['directed'])
        self.assertIsNone(result['log2_upper_rational'])
        self.assertTrue(result['pending_frontier'])
        with self.assertRaises(ValueError):
            composite.production_envelope(GenericEvaluator(tables), result, controller)

    def test_bank_cycles_reclaimed_and_requested_not_calibration_weights(self):
        refs, ranges = [], []
        detached = table(first=3, last=4)
        class Bank:
            def __init__(self, *args, **kwargs):
                if any(ref() is not None for ref in refs):
                    raise AssertionError('Previous bank retained')
                refs.append(weakref.ref(self))
            def metadata(self):
                return {}
            def check_pins(self):
                pass
        class Evaluator:
            def __init__(self, bank, first, last, **kwargs):
                self.bank, self.tables, self.cycle = bank, detached, self
                ranges.append((first, last))
            def context(self):
                return {'requested_weights': [3, 4]}
            def check_pins(self):
                pass
        fake = SimpleNamespace(cert=SimpleNamespace(AliasBank=Bank),
            bands=SimpleNamespace(Evaluator=Evaluator, BandTables=bands.BandTables))
        proposal = dict(witness=dict(z='1/2', alpha='1/2', markers=[['1', '1']]),
                        source=dict(path='tiny.json', sha256='0'*64), calibration_weight=2)
        enabled = gc.isenabled()
        gc.disable()
        try:
            for unused in range(2):
                tables, context, metadata = composite.prepare_component(proposal, 3, 4, fake, base.Budget())
                self.assertIs(tables, detached)
                self.assertEqual(metadata['calibration_weight'], 2)
                self.assertEqual(metadata['requested_weight_range'], [3, 4])
                self.assertTrue(all(ref() is None for ref in refs))
        finally:
            if enabled:
                gc.enable()
        self.assertEqual(ranges, [(3, 4), (3, 4)])

    def test_limits_frozen_sources_and_strict_production_identity(self):
        for args in ((0, 1, 1, 54, 48), (1, 65, 1, 54, 48), (8191, 8192, 1, 54, 48),
                     (3, 2, 1, 54, 48), (1, 2, 4097, 54, 48), (1, 2, 1, 48, 54), (1, 2, 1, 54, 39)):
            with self.assertRaises(ValueError):
                composite.validate_request(*args)
        composite.validate_request(540, 575, 4096, 54, 48)
        composite.validate_request(8128, 8191, 512, 64, 40)
        self.assertEqual(composite.MAX_OUTPUT_BYTES, 32 << 20)
        self.assertEqual(controller.RESERVE, 20 << 30)
        self.assertEqual(controller.MEMORY_BYTES, 4 << 30)
        self.assertEqual(composite.selected.HARD_SECONDS, 600)
        self.assertEqual(sha256((composite.HERE/composite.selected.NAME).read_bytes()).hexdigest(), composite.SELECTED_SHA)
        self.assertEqual(sha256((composite.HERE/(composite.CONTROLLER_NAME+'.py')).read_bytes()).hexdigest(), composite.CONTROLLER_SHA)
        with self.assertRaises(ValueError):
            composite.production_envelope(object(), {}, controller)


if __name__ == '__main__':
    unittest.main()

"""Tiny exhaustive geometry, exact transfer, and mocked-bank lifecycle tests."""

from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from flint import arb, ctx

import hd_field128_gamma_midpoint_weight_bands as mode
from hd_shifted_pairs_count import exact_profiles


def multiply(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h


def power(matrix, exponent):
    result = (F(1), F(0), F(0), F(1))
    for unused in range(exponent):
        result = multiply(result, matrix)
    return result


def rational_row(weight, full):
    return (F(weight+full+1, 32), F(full+2, 64),
            F(weight+2, 64), F(2*full+3, 64))


def regional(pairs):
    def get(weight, budget):
        return {f: tuple(mode.base.rational(x) for x in rational_row(weight, f))
                for f in mode.tail.full_counts(pairs, weight)}
    return get


def tables(pairs=4, first=2, last=6, precision=256):
    return mode.BandTables(pairs, first, last, regional(pairs), arb(-400), arb(4),
                           precision=precision, repetitions=4)


def cells(root):
    intervals = [[(low, high) for low in range(lo, hi+1) for high in range(low, hi+1)]
                 for lo, hi in root.ranges]
    return (mode.Rectangle(root.parity, tuple(rows)) for rows in product(*intervals))


def profiles_in(table, cell):
    (ul, uh), (dl, dh), (el, eh) = cell.ranges
    for r in table.weights:
        u, eps = divmod(r, 2)
        if eps != cell.parity or not ul <= u <= uh:
            continue
        limit = mode.deficit_limit(table.pairs, r)
        for d in range(dl, min(dh, limit)+1):
            for e in range(el, min(eh, limit)+1):
                yield r, u-d, u-e


def fake_evaluator(table, split_to_points=False):
    original = table.receipt
    if split_to_points:
        def receipt(cell, bits, budget):
            result = original(cell, bits, budget)
            if any(lo != hi for lo, hi in cell.ranges) and not result['zero']:
                result.update(log2_upper_rational='0', log2_upper_enclosure='0',
                              meets_cell_budget=False)
            return result
        table.receipt = receipt
    return SimpleNamespace(tables=table, _active=False, check_pins=Mock(),
        context=lambda: dict(requested_weight_range=[table.weights[0], table.weights[-1]],
                             requested_weights=list(table.weights), synthetic_test_only=True))


def production_fixture(witness_source=None):
    bank = object.__new__(mode.provider.AliasBank)
    bank.precision, bank.constant, bank.slope = 256, arb(-1000), arb(2)
    bank.state, bank.q, bank.distance = 30, 1 << 128, F(2, 5)
    bank.outer_measure = 'invertible-toeplitz'
    bank.pins, bank.check_pins = {'fixture': 'not-production-data'}, Mock()
    bank.metadata = lambda: dict(construction='synthetic-test-only', state_dimension=30,
        target_distance='2/5', witness=dict(precision=256), source_sha256=bank.pins)
    with patch.object(mode.base, 'PAIRS', 4), patch.object(mode.restored, 'RegionalCache') as cache:
        cache.return_value.regional.side_effect = regional(4)
        evaluator = mode.Evaluator(bank, 2, 6, witness_source=witness_source)
        calls = [call.args[0] for call in cache.return_value.regional.call_args_list]
    bank.check_pins.reset_mock()
    return evaluator, bank, calls


class BandTests(unittest.TestCase):
    def test_exhaustive_domains_counts_and_integer_splits(self):
        for p in range(1, 5):
            table = tables(p, 1, 2*p-1)
            exact = exact_profiles(p)
            seen = set()
            for root in table.roots:
                for point in profiles_in(table, root):
                    self.assertNotIn(point, seen)
                    seen.add(point)
                for cell in cells(root):
                    expected = actual = 0
                    for r, a, b in profiles_in(table, cell):
                        count = next(row[3].counts for row in table.count_tables[r % 2] if row[0] == r)
                        expected += count.midpoint_upper(a, b)
                        actual += exact.get((r, a, b), 0)
                    value = table.count(cell, mode.base.Budget())
                    self.assertEqual(value, expected)
                    self.assertGreaterEqual(value, actual)
                    for axis, (lo, hi) in enumerate(cell.ranges):
                        if lo < hi:
                            left, right = cell.split(axis, (lo+hi)//2)
                            self.assertEqual(value, table.count(left, mode.base.Budget())+
                                             table.count(right, mode.base.Budget()))
            expected = {(r, a, b) for r in table.weights
                        for a in mode.tail.full_counts(p, r) for b in mode.tail.full_counts(p, r)}
            self.assertEqual(seen, expected)
            self.assertTrue({point for point, value in exact.items() if 1 <= point[0] < 2*p and value} <= seen)

    def test_noncommuting_transfer_dominates_explicit_integer_sum(self):
        for precision in (128, 256, 384):
            table = tables(precision=precision)
            with ctx.workprec(precision):
                for cell in table.roots:
                    exact = F(0)
                    for r, a, b in profiles_in(table, cell):
                        count = next(row[3].counts for row in table.count_tables[r % 2] if row[0] == r)
                        raw = power(multiply(rational_row(r, a), rational_row(r, b)), 2)
                        exact += count.midpoint_upper(a, b)*F(2)**(-400+4*r)*(raw[0]+raw[1])
                    result = table.receipt(cell, 10, mode.base.Budget())
                    self.assertGreaterEqual(F(result['log2_upper_rational']),
                        mode.base.rational_ceiling(mode.tail.log2arb(mode.base.rational(exact))))
                # A singleton fixes the original r and chronological even/odd matrices.
                point = mode.Rectangle(0, ((2, 2), (0, 0), (1, 1)))
                result = table.receipt(point, 10, mode.base.Budget())
                count = next(row[3].counts for row in table.count_tables[0] if row[0] == 4)
                forward = power(multiply(rational_row(4, 2), rational_row(4, 1)), 2)
                backward = power(multiply(rational_row(4, 1), rational_row(4, 2)), 2)
                self.assertNotEqual(forward[0]+forward[1], backward[0]+backward[1])
                exact = count.midpoint_upper(2, 1)*F(2)**(-400+16)*(forward[0]+forward[1])
                ceiling = mode.base.rational_ceiling(mode.tail.log2arb(mode.base.rational(exact)))
                self.assertLessEqual(abs(F(result['log2_upper_rational'])-ceiling), F(2, 10**12))

    def test_zero_impossible_faces_and_upper_half_original_weights(self):
        table = tables(5, 4, 9)
        for r in (4, 6, 8):
            u = r//2
            cell = mode.Rectangle(0, ((u, u), (0, 0), (0, 0)))
            result = table.receipt(cell, 100, mode.base.Budget())
            self.assertTrue(result['zero'])
            self.assertEqual(result['count_upper'], 0)
        alternate = tables(4, 4, 4)
        cell = mode.Rectangle(0, ((2, 2), (2, 2), (2, 2)))
        self.assertEqual(alternate.count(cell, mode.base.Budget()), 2)
        # At original r=9, u=4 but only d=e=0 are marginally feasible.
        cell = mode.Rectangle(1, ((4, 4), (1, 1), (1, 1)))
        self.assertTrue(table.receipt(cell, 100, mode.base.Budget())['zero'])
        calls = []
        def get(r, budget):
            calls.append(r)
            return regional(5)(r, budget)
        mode.BandTables(5, 6, 9, get, arb(-400), arb(4), repetitions=4)
        self.assertEqual(calls, [6, 7, 8, 9])

    def test_complete_tree_partial_frontier_and_corruption(self):
        value = fake_evaluator(tables(3, 2, 4), split_to_points=True)
        saved, contexts = [], []
        result = mode.cover(value, max_nodes=1000, cell_failure_bits=10,
                            checkpoint=saved.append, context_checkpoint=contexts.append)
        self.assertTrue(result['complete_requested_band'])
        self.assertFalse(result['directed'])
        self.assertTrue(result['meets_requested_band_target'])
        self.assertEqual(value.check_pins.call_count, 2)
        self.assertTrue(any(row['action'] == 'integer-split' for row in saved))
        replay = {}
        for row in saved:
            self.assertFalse(row['complete_requested_band'])
            self.assertFalse(row['directed'])
            self.assertIsNone(row['log2_upper_rational'])
            replay.update((node['id'], node) for node in row['changed_nodes'])
        self.assertEqual([replay[i] for i in sorted(replay)], result['nodes'])
        for mutation in ('missing', 'duplicate', 'context'):
            bad = deepcopy(result['nodes'])
            if mutation == 'missing':
                bad.pop()
            elif mutation == 'duplicate':
                bad.append(deepcopy(bad[0]))
            else:
                next(row for row in bad if row['status'] == 'accepted')['receipt']['context_sha256'] = 'wrong'
            with self.assertRaises(ValueError):
                mode.audit_tree(bad, value.tables.roots, 10, result['context_sha256'])
        value = fake_evaluator(tables(3, 2, 4), split_to_points=True)
        partial = mode.cover(value, max_nodes=1, cell_failure_bits=10)
        self.assertFalse(partial['complete_requested_band'])
        self.assertTrue(partial['pending_frontier'])
        self.assertIsNone(partial['log2_upper_rational'])
        self.assertFalse(partial['directed'])

    def test_actual_entry_uses_all_weights_and_amortized_guards(self):
        value, bank, calls = production_fixture()
        self.assertEqual(calls, [2, 3, 4, 5, 6])
        saved = []
        result = value.run_cover(max_nodes=10, checkpoint=saved.append)
        self.assertTrue(result['directed'])
        self.assertTrue(result['complete_requested_band'])
        self.assertEqual(bank.check_pins.call_count, 2)
        self.assertEqual(result['requested_weights'], calls)
        self.assertEqual(result['context']['repetitions'], 256)
        self.assertIn(mode.NOTE, result['context']['source_sha256'])
        self.assertTrue(all('source_sha256' not in row for row in saved))
        self.assertTrue(all('source_sha256' not in node.get('receipt', {}) for node in result['nodes']))
        self.assertFalse(value._active)

    def test_source_witness_mutation_and_exception_guards(self):
        value, bank, unused = production_fixture()
        read, changed = Path.read_bytes, [False]
        def altered(path):
            raw = read(path)
            return raw+b'\n# synthetic mutation\n' if changed[0] and path.name == Path(mode.__file__).name else raw
        def mutate(row):
            changed[0] = True
        with patch.object(Path, 'read_bytes', altered):
            with self.assertRaisesRegex(RuntimeError, 'pinned'):
                value.run_cover(checkpoint=mutate)
        self.assertFalse(value._active)
        # The failed source-byte guard precedes the redundant bank guard.
        self.assertEqual(bank.check_pins.call_count, 1)
        path = Path(mode.__file__).resolve()
        source = dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest())
        value, bank, unused = production_fixture(source)
        changed[0] = False
        # Directly use a different already-existing file as witness so this
        # checks witness bytes separately from the module closure.
        witness = Path(__file__).resolve()
        value._witness_pin = witness, sha256(witness.read_bytes()).hexdigest()
        def witness_change(path):
            raw = read(path)
            return raw+b'changed' if changed[0] and path == witness else raw
        with patch.object(Path, 'read_bytes', witness_change):
            with self.assertRaisesRegex(RuntimeError, 'witness'):
                value.run_cover(checkpoint=mutate)
        self.assertFalse(value._active)
        value, bank, unused = production_fixture()
        def reenter(row):
            value.run_cover()
        with self.assertRaisesRegex(RuntimeError, 'already active'):
            value.run_cover(checkpoint=reenter)
        self.assertFalse(value._active)
        self.assertEqual(bank.check_pins.call_count, 2)

    def test_invalid_inputs_and_negative_nonfinite_caps(self):
        for args in ((4, 0, 2), (4, 2, 8), (100, 1, 65), (4, True, 2)):
            with self.assertRaises(ValueError):
                mode.checked_weights(*args)
        for args in ((2, ((0, 0),)*3), (0, ((0, 0),)*2), (0, ((False, 0),)*3)):
            with self.assertRaises(ValueError):
                mode.Rectangle(*args)
        for entry in (arb(-1), arb('nan')):
            def get(r, budget):
                return {0: (entry, arb(0), arb(0), arb(0))}
            with self.assertRaises(ArithmeticError):
                mode.BandTables(2, 1, 1, get, arb(0), arb(0))
        with self.assertRaises(ValueError):
            mode.Evaluator(object(), 1, 2)
        value = fake_evaluator(tables())
        for kwargs in (dict(max_nodes=0), dict(max_seconds=-1), dict(cell_failure_bits=True), dict(checkpoint=1)):
            with self.assertRaises(ValueError):
                mode.cover(value, **kwargs)


if __name__ == '__main__':
    unittest.main()

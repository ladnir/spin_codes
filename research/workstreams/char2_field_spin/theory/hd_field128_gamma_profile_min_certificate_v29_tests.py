"""Tiny exact directed controls only: no actual bank, search or disk output."""
from contextlib import nullcontext
from copy import deepcopy
from fractions import Fraction as F
import gc
from hashlib import sha256
from types import SimpleNamespace
import unittest
from unittest.mock import MagicMock, patch
import weakref

from flint import arb, ctx

import hd_field128_gamma_profile_min_certificate_v29 as producer
import hd_field128_gamma_midpoint_certificate as cert


def multiply(a, b):
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])


def power(a, n):
    out = (F(1), F(0), F(0), F(1))
    for unused in range(n):
        out = multiply(out, a)
    return out


def fixture(pairs, weight, bank):
    return {f: tuple(F(x, 16) for x in
            ((3+f, 1+bank, 2+bank+f, 4) if bank == 0 else (2, 5+f, 1, 3+f)))
            for f in cert.tail.full_counts(pairs, weight)}


def as_arb(rows):
    return {f: tuple(cert.base.rational(x) for x in row) for f, row in rows.items()}


def exact_union(banks, prefixes, counts, repetitions):
    total, winners = F(0), [0]*len(banks)
    singles = [F(0)]*len(banks)
    for a in counts.full:
        for b in counts.full:
            count = counts.count(a, b)
            if not count:
                continue
            values = []
            for index, bank in enumerate(banks):
                matrix = power(multiply(bank[a], bank[b]), repetitions//2)
                value = count*prefixes[index]*(matrix[0]+matrix[1])
                values.append(value)
                singles[index] += value
            chosen = min(range(len(values)), key=values.__getitem__)
            total += values[chosen]
            winners[chosen] += 1
    return total, singles, winners


class ProfileMinimumDirectedTests(unittest.TestCase):
    def test_fraction_oracle_own_prefix_order_and_complete_partition(self):
        with ctx.workprec(256):
            for pairs, weight, repetitions in ((3, 2, 4), (4, 4, 8), (3, 4, 8)):
                counts = cert.midpoint.MidpointCounts(pairs, weight)
                banks, prefixes = [fixture(pairs, weight, i) for i in range(2)], [F(3, 8), F(5, 4)]
                oracle, singles, winners = exact_union(banks, prefixes, counts, repetitions)
                value, report = producer.directed_rows(list(map(as_arb, banks)),
                    list(map(cert.base.rational, prefixes)), pairs, weight, counts, cert,
                    repetitions=repetitions)
                self.assertTrue(value >= cert.base.rational(oracle))
                self.assertTrue(value < cert.base.rational(oracle)*(1+arb(2)**-180))
                self.assertEqual(report['witness_winning_profile_counts'], winners)
                self.assertEqual(report['profile_pairs_evaluated'], len(counts.full)**2)
                self.assertLessEqual(oracle, min(singles))
                self.assertTrue(all(row['outer_prefactors_already_included'] for row in report['rows']))
                fresh = producer.aggregate_rows(report['rows'], counts.full, cert)
                self.assertFalse(fresh['zero'])
                self.assertTrue(arb(2)**cert.base.rational(F(fresh['log2_upper_rational'])) >= cert.base.rational(oracle))

    def test_minimum_must_follow_complete_transfer_and_own_prefix(self):
        pairs, weight, reps = 2, 2, 4
        counts = cert.midpoint.MidpointCounts(pairs, weight)
        # Complementary diagonal caps make entry-wise mixing spuriously small.
        banks = [{f: tuple(map(F, (1, 0, 0, 4))) for f in counts.full},
                 {f: tuple(map(F, (4, 0, 0, 1))) for f in counts.full}]
        prefixes = [F(8), F(1, 8)]
        oracle, unused, wins = exact_union(banks, prefixes, counts, reps)
        wrong = {f: tuple(min(a, b) for a, b in zip(banks[0][f], banks[1][f])) for f in counts.full}
        splice, unused, unused_wins = exact_union([wrong], [min(prefixes)], counts, reps)
        self.assertGreater(oracle, splice)
        with ctx.workprec(256):
            value, report = producer.directed_rows(list(map(as_arb, banks)),
                list(map(cert.base.rational, prefixes)), pairs, weight, counts, cert, repetitions=reps)
            self.assertEqual(value, cert.base.rational(oracle))
            self.assertEqual(report['witness_winning_profile_counts'], wins)
        # Chronology is independently visible for the noncommuting toy matrices.
        toy = fixture(3, 2, 0)
        self.assertNotEqual(multiply(toy[0], toy[1]), multiply(toy[1], toy[0]))

    def test_single_witness_replays_frozen_transfer_256(self):
        with ctx.workprec(256):
            counts = cert.midpoint.MidpointCounts(2, 2)
            rows, prefix = as_arb(fixture(2, 2, 0)), arb(8)
            actual, report = producer.directed_rows([rows], [prefix], 2, 2, counts, cert)
            control, old = cert.transfer_sum(2, 2, rows, counts=counts)
            self.assertTrue(actual >= prefix*control*(1-arb(2)**-180))
            self.assertTrue(actual <= prefix*control*(1+arb(2)**-180))
            self.assertEqual(report['profile_pairs_evaluated'], old['profile_pairs_evaluated'])

    def test_zeros_and_invalid_counts_are_fail_closed(self):
        counts = cert.midpoint.MidpointCounts(3, 2)
        zeros = {f: (arb(0),)*4 for f in counts.full}
        value, report = producer.directed_rows([zeros], [arb(1)], 3, 2, counts, cert)
        self.assertTrue(value.is_zero())
        self.assertTrue(producer.aggregate_rows(report['rows'], counts.full, cert)['zero'])
        self.assertEqual(report['profile_pairs_evaluated'], 4)
        with patch.object(counts, 'count', return_value=-1):
            with self.assertRaises(cert.IncompleteTransfer) as caught:
                producer.directed_rows([zeros], [arb(1)], 3, 2, counts, cert)
        self.assertEqual(caught.exception.partial['completed_even_rows'], 0)
        for bad in (arb(0), arb(-1), arb('inf')):
            with self.assertRaises((ValueError, ArithmeticError)):
                producer.directed_rows([zeros], [bad], 3, 2, counts, cert)

    def test_stop_preserves_only_completed_rows_no_aggregate(self):
        class Budget:
            expired = False
            def check(self):
                if self.expired:
                    raise TimeoutError('tiny stop')
        budget, saved = Budget(), []
        def stop(row):
            saved.append(row)
            budget.expired = True
        counts = cert.midpoint.MidpointCounts(4, 4)
        rows = as_arb(fixture(4, 4, 0))
        with self.assertRaises(cert.IncompleteTransfer) as caught:
            producer.directed_rows([rows], [arb(1)], 4, 4, counts, cert,
                                   repetitions=4, budget=budget, checkpoint=stop)
        partial = caught.exception.partial
        self.assertEqual(partial['completed_even_rows'], 1)
        self.assertEqual(partial['pending_even_full_range'], [1, 2])
        self.assertEqual(len(saved), 1)
        self.assertFalse(saved[0]['complete_selected_weight'])
        with self.assertRaises(ValueError):
            producer.aggregate_rows(partial['rows'], counts.full, cert)

    def test_saved_row_partition_guards_and_no_double_prefix(self):
        counts = cert.midpoint.MidpointCounts(3, 2)
        rows = as_arb(fixture(3, 2, 0))
        with ctx.workprec(256):
            value, report = producer.directed_rows([rows], [arb(32)], 3, 2, counts, cert, repetitions=4)
        complete = report['rows']
        for altered in (complete[::-1], complete+[complete[0]], complete[:-1]):
            with self.assertRaises(ValueError):
                producer.aggregate_rows(altered, counts.full, cert)
        altered = deepcopy(complete)
        altered[0]['outer_prefactors_already_included'] = False
        with self.assertRaises(ValueError):
            producer.aggregate_rows(altered, counts.full, cert)
        aggregate = producer.aggregate_rows(complete, counts.full, cert)
        reconstructed = arb(2)**cert.base.rational(F(aggregate['log2_upper_rational']))
        self.assertTrue(reconstructed >= value)
        self.assertTrue(reconstructed < value*arb('1.00000001'))

    def test_literal_only_scope_frozen_pins_and_limits(self):
        self.assertEqual(sha256((producer.HERE/producer.literals.NAME).read_bytes()).hexdigest(), producer.LITERALS_SHA)
        self.assertEqual(sha256((producer.HERE/(producer.CERT_NAME+'.py')).read_bytes()).hexdigest(), producer.CERT_SHA)
        self.assertEqual(producer.HARD_SECONDS, 600)
        self.assertEqual(producer.TABLE_BYTES, 64 << 20)
        self.assertEqual(producer.portable.MAX_OUTPUT_BYTES, 1 << 20)
        self.assertEqual(producer.portable.RESERVE, 20 << 30)
        with patch.object(producer.portable, 'check_pins', side_effect=RuntimeError('changed')):
            with self.assertRaises(RuntimeError):
                producer.check_pins({'f.py': 'a'*64}, [{'source': {'path': 'input.json', 'sha256': 'b'*64}}])
        with patch.object(producer, 'remote_limits', return_value=nullcontext()), \
                patch.object(producer, 'helpers', side_effect=AssertionError('Invalid target must stop first')):
            with self.assertRaises(ValueError):
                producer.run([], 8192, 'unused')

    def test_sequential_banks_release_frozen_cache_cycles(self):
        references = []
        class Bank:
            def __init__(self, *args, **kwargs):
                if any(ref() is not None for ref in references):
                    raise AssertionError('A previous bank remains live')
                references.append(weakref.ref(self))
                self.constant, self.slope = arb(0), arb(0)
            def check_pins(self):
                pass
            def metadata(self):
                return {}
        def normalize(bank, weight, budget):
            cycle = SimpleNamespace(bank=bank)
            cycle.self = cycle
            return as_arb(fixture(3, 2, 0))
        proposal = dict(witness=dict(z='1/2', alpha='1/2', markers=[['1', '1']]),
                        source={'path': 'tiny.json', 'sha256': '0'*64}, source_marker_key='markers')
        enabled = gc.isenabled()
        gc.disable()
        try:
            with patch.object(cert.provider, 'AliasBank', Bank), \
                    patch.object(cert.frozen, 'regional_matrices', new=normalize):
                for unused in range(2):
                    rows, factor, metadata = producer.build_regions(proposal, 2, cert, cert.base.Budget())
                    self.assertEqual(factor, arb(1))
                    self.assertTrue(metadata['fresh_bank_rebuilt'])
                    self.assertTrue(all(ref() is None for ref in references))
                    self.assertEqual(set(rows), {0, 1})
        finally:
            if enabled:
                gc.enable()

    def test_output_reserve_and_cumulative_guards_before_open(self):
        directory, path = MagicMock(), MagicMock()
        directory.__truediv__.return_value = path
        path.exists.return_value = False
        output = producer.portable.Outputs(directory, lambda obj, unused: obj)
        with patch.object(producer.portable, 'reserve_guard', side_effect=OSError('reserve')):
            with self.assertRaises(OSError):
                output.write('row.json', {'row': 1})
        path.open.assert_not_called()
        self.assertEqual(output.used, 0)
        output.used = producer.portable.MAX_OUTPUT_BYTES
        with patch.object(producer.portable, 'reserve_guard') as guard:
            with self.assertRaises(OSError):
                output.write('row.json', {'row': 1})
            guard.assert_not_called()
        path.open.assert_not_called()


if __name__ == '__main__':
    unittest.main()

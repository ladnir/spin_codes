"""Independent small controls for the high-total-arity directed adapter."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_two_sum_gf16_high_arity_certificate as target
from hd_shifted_pairs_count import exact_profiles


def multiply(left, right):
    return tuple(left[2*i]*right[j]+left[2*i+1]*right[2+j]
                 for i in range(2) for j in range(2))


def power(value, exponent):
    result = (1, 0, 0, 1)
    for _ in range(exponent):
        result = multiply(result, value)
    return result


def fake_regions(mark, alpha, maximum, bits):
    assert alpha == F(1, 5) and bits == 60
    return {(maximum, f): tuple(arb(x)/arb(2)**128 for x in (2, 1, 3, 2))
            for f in range(maximum//2+1)}


def fake_receipt(precision):
    with patch.object(target, 'checked_dependencies'), \
            patch.object(target.uniform, 'regional_directed', side_effect=fake_regions):
        result = target.receipt(precision)
    result['source_sha256'] = sha256(Path(target.__file__).read_bytes()).hexdigest()
    return result


class HighArityTests(unittest.TestCase):
    def test_reviewed_dependencies_and_scope(self):
        target.checked_dependencies()
        self.assertEqual([len(target.profile_counts(r)) for r, _, _ in target.ANCHORS], [99, 168, 288])
        self.assertFalse(target.EXPECTED_SCOPE['whole_code_certificate'])
        self.assertFalse(target.EXPECTED_SCOPE['entire_active_group_weight_classes_covered'])
        self.assertEqual(target.EXPECTED_SCOPE['rowset_count_factor'], 1)
        self.assertEqual(target.EXPECTED_SCOPE['regional_amplitude_measure'], 'uniform-field-including-zero')

    def test_direct_density_prefix_equals_legacy_plus_exact_deficit(self):
        for precision in (256, 384):
            with ctx.workprec(precision):
                for r, a0, mark in target.ANCHORS:
                    alpha = target.uniform.rational(target.ALPHA)
                    legacy = alpha*(r*(2048-128*target.log2(arb(255)))
                        -target.log2(arb(255))-734003*target.log2(target.uniform.rational(mark)))
                    delta = 8*alpha*(16*r-a0)*target.log2(arb(255))
                    self.assertTrue((target.high_prefix(r, a0, mark)-legacy-delta).contains(0))
                    self.assertTrue((target.high_prefix(r, 16*r, mark)-legacy).contains(0))
        for args in ((True, 16, F(1, 2)), (1, 17, F(1, 2)), (2, 1, F(1, 2)), (1, 1, 1)):
            with self.assertRaises(ValueError):
                target.high_prefix(*args)

    def test_ordered_profile_sum_dominates_exact_indexed_cycles(self):
        values = {(3, 0): (2, 1, 3, 2), (3, 1): (1, 4, 2, 1)}
        self.assertNotEqual(multiply(values[3, 0], values[3, 1]), multiply(values[3, 1], values[3, 0]))
        exact = exact_profiles(3)
        expected_upper = 0
        for (a, b), count in target.profile_counts(3, 3).items():
            self.assertGreaterEqual(count, exact.get((3, a, b), 0))
            run = power(multiply(values[3, a], values[3, b]), 2)
            expected_upper += count*(run[0]+run[1])
        actual = sum(count*sum(power(multiply(values[3, a], values[3, b]), 2)[:2])
                     for (r, a, b), count in exact.items() if r == 3)
        with ctx.workprec(256):
            terms = target.source_profile_terms(3, {key: tuple(map(arb, row)) for key, row in values.items()},
                                                pairs=3, cycles=2)
            bounded = sum((row[3] for row in terms), arb(0))
            self.assertTrue(bounded.contains(expected_upper))
            self.assertGreaterEqual(expected_upper, actual)
        with self.assertRaises(ValueError):
            target.source_profile_terms(3, {}, pairs=3, cycles=2)

    def test_unique_row_set_comparison_over_gf4(self):
        # Two repetition RS[2,1] rows. H=[[1,1],[1,2]] is invertible over
        # GF4. Condition on each output diagonal, and enumerate all input
        # diagonals/messages. No source-row-set count multiplies the cap.
        def mul(a, b):
            out = 0
            while b:
                if b & 1:
                    out ^= a
                b >>= 1
                a <<= 1
                if a & 4:
                    a ^= 7
            return out
        prescaled = Counter()
        for diagonal in product(range(1, 4), repeat=4):
            for a, b in product(range(4), repeat=2):
                prescaled[tuple(mul(x, y) for x, y in zip(diagonal, (a, a, b, b)))] += 1
        averaged = {1: Counter(), 2: Counter()}
        for output_diagonal in product(range(1, 4), repeat=4):
            determined_rows = {}
            conditional = {1: Counter(), 2: Counter()}
            for v, count in prescaled.items():
                rows = (any(v[:2]), any(v[2:]))
                mixed = (v[0]^v[2], v[1]^v[3], v[0]^mul(2, v[2]), v[1]^mul(2, v[3]))
                y = tuple(mul(x, d) for x, d in zip(mixed, output_diagonal))
                self.assertNotIn(y, determined_rows)
                determined_rows[y] = rows
                for a0 in (1, 2):
                    if sum(rows) >= a0:
                        conditional[a0][y] += count
                        averaged[a0][y] += count
            for a0 in (1, 2):
                self.assertLessEqual(max(conditional[a0].values()), 81//3**a0)
        for a0 in (1, 2):
            self.assertLessEqual(max(averaged[a0].values()), 81**2//3**a0)

    def test_fresh_banks_once_per_anchor_at_requested_precision(self):
        calls = []
        def bank(*args):
            calls.append((ctx.prec, args))
            return fake_regions(*args)
        with patch.object(target, 'checked_dependencies'), \
                patch.object(target.uniform, 'regional_directed', side_effect=bank), \
                patch.object(target.uniform, 'weight_receipt', side_effect=AssertionError('Saved legacy receipt path')):
            result = target.receipt(256)
        self.assertEqual(calls, [(256, (z, target.ALPHA, r, 60)) for r, _, z in target.ANCHORS])
        self.assertEqual(len(result['entries']), 3)
        self.assertEqual(sum(len(row['profile_entries']) for row in result['entries']), 555)
        self.assertFalse(result['saved_receipts_used'])
        self.assertFalse(result['sharper_toeplitz_measure_used'])

    def test_scope_and_profile_validation_fail_closed(self):
        fixture = fake_receipt(256)
        target.checked_rows(fixture, 256)
        mutations = [lambda x: x.update(rowset_count_factor=16),
            lambda x: x.update(whole_code_certificate=True),
            lambda x: x['entries'].append(deepcopy(x['entries'][0])),
            lambda x: x['entries'][0].update(total_arity_minimum=246),
            lambda x: x['entries'][0]['profile_entries'].pop(),
            lambda x: x['entries'][0]['profile_entries'][0].update(indexed_source_count_upper=1),
            lambda x: x.update(source_sha256='unreviewed')]
        for mutate in mutations:
            bad = deepcopy(fixture)
            mutate(bad)
            with self.assertRaises(ValueError):
                target.checked_rows(bad, 256)

    def test_dual_precision_verifier_and_independent_rounded_sum(self):
        fixtures = {bits: fake_receipt(bits) for bits in (256, 384)}
        with tempfile.TemporaryDirectory() as directory:
            for bits, data in fixtures.items():
                (Path(directory)/f'hd_two_sum_gf16_high_arity_receipt{bits}.json').write_text(json.dumps(data))
            verified = target.verify(directory)
            exact = target.rounded_union([F(row['failure_log2_upper']) for row in fixtures[256]['entries']])
            self.assertEqual(verified['rounded_anchor_union_log2_upper'], str(exact))
            self.assertEqual(verified['ordered_profile_ceilings_checked'], 555)
            self.assertEqual(verified['independent_sum_precision_bits'], 512)
            bad = fixtures[384]
            bad['entries'][0]['profile_entries'][0]['failure_log2_upper'] = '123'
            (Path(directory)/'hd_two_sum_gf16_high_arity_receipt384.json').write_text(json.dumps(bad))
            with self.assertRaises(ValueError):
                target.verify(directory)

    def test_guards_and_cooperative_timeout(self):
        for kwargs in ({'precision': True}, {'precision': 512}, {'max_seconds_per_anchor': True},
                       {'max_seconds_per_anchor': 61}, {'max_seconds_per_anchor': float('nan')}):
            with patch.object(target.uniform, 'regional_directed') as bank, self.assertRaises(ValueError):
                target.receipt(**kwargs)
            bank.assert_not_called()
        with patch.object(target, 'checked_dependencies'), \
                patch.object(target, 'perf_counter', side_effect=(0., 0., 61.)), \
                patch.object(target.uniform, 'regional_directed') as bank, self.assertRaises(TimeoutError):
            target.receipt(256)
        bank.assert_not_called()


if __name__ == '__main__':
    unittest.main()

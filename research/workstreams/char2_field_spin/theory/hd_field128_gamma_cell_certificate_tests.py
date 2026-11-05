"""Arithmetic tests of the cell checker; synthetic banks prove no code claim."""

from fractions import Fraction as F
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_gamma_cell_certificate as c


class CellArithmeticTests(unittest.TestCase):
    def bank(self):
        markers = ((F(1, 2), F(1, 3)), (F(2, 3), F(3, 5)))
        logs = tuple(tuple(c.tail.log2arb(c.base.rational(x)) for x in row) for row in markers)
        bank = SimpleNamespace(markers=markers, log_markers=logs, constant=arb(-100), slope=arb(3))
        bank.scalar = lambda point, indices, count_matching=None: c.tail.PositiveBank.scalar(
            bank, point, indices, count_matching=count_matching)
        bank.numerator = lambda index, rect: tuple(arb(v) for v in ((1, 2, 3, 4), (2, 1, 4, 3))[index])
        return bank

    def provider_stub(self):
        """Test-only instance: bypass coefficient construction, with synthetic caps."""
        toy = self.bank()
        bank = object.__new__(c.provider.AliasBank)
        bank.__dict__.update(toy.__dict__)
        bank.__dict__.update(construction='synthetic-review-only', state=30, q=1 << 128,
                             distance=F(2, 5), forbidden=838860, precision=256,
                             z=F(1, 2), alpha=F(1, 4), grid=64, bits=80, pins={})
        bank.check_pins = lambda: None
        return bank

    def test_scalar_none_equals_independent_phi(self):
        with ctx.workprec(256):
            bank = self.bank()
            for point in ((F(17, 2), F(3, 2), F(5, 2)), (20, 10, 0), (5000, 3000, 2000)):
                for indices in ((0, 0), (0, 1), (1, 0)):
                    r, a, b = point
                    expected = bank.constant+bank.slope*c.base.rational(r)
                    for full, index in zip((a, b), indices):
                        e, s = 4096-r+full, r-2*full
                        logh = (arb(4097).lgamma()-sum((c.base.rational(x)+1).lgamma()
                                                                     for x in (e, s, full)))/arb(2).log()
                        expected -= 128*(logh+c.base.rational(s)*bank.log_markers[index][0]
                                         +c.base.rational(full)*bank.log_markers[index][1])
                    self.assertTrue(bank.scalar(point, indices, count_matching=None).overlaps(expected))

    def test_single_point_matches_directed_formula(self):
        with ctx.workprec(256):
            bank = self.bank()
            cell = c.Cell('interior', ((9, 9), (2, 2), (3, 3)))
            result = c.fixed_bound(bank, cell, (0, 1), 37, c.base.Budget())
            left, right = bank.numerator(0, None), bank.numerator(1, None)
            direct = (c.tail.log2arb(arb(37))+c.tail.log2arb(c.profiles.pair_row_sum128(left, right))
                      +bank.scalar((9, 2, 3), (0, 1), count_matching=None))
            self.assertGreaterEqual(F(result['log2_upper_rational']), c.base.rational_ceiling(direct))

    def test_vertex_envelope_bounds_rational_interior(self):
        with ctx.workprec(256):
            bank = self.bank()
            for cell in (c.Cell('interior', ((7, 10), (0, 4), (1, 4))),
                         c.Cell('even_face', ((4, 7), (1, 4))), c.Cell('odd_face', ((4, 7), (1, 4)))):
                result = c.fixed_bound(bank, cell, (0, 1), 1, c.base.Budget())
                ceiling = c.base.rational(F(result['scalar_log2_upper_rational']))
                vertices = cell.vertices()
                for first in vertices:
                    for second in vertices:
                        point = tuple((x+2*y)/3 for x, y in zip(first, second))
                        self.assertLessEqual(bank.scalar(cell.profile(point), (0, 1)).upper(), ceiling)

    def test_numerator_rectangles_cover_all_integer_profiles(self):
        with ctx.workprec(256):
            bank = self.bank()
            cell = c.Cell('interior', ((7, 10), (0, 4), (1, 4)))
            seen = []
            bank.numerator = lambda index, rect: seen.append((index, rect)) or (arb(1),)*4
            c.fixed_bound(bank, cell, (0, 1), 1, c.base.Budget())
            for r in range(7, 11):
                for a in range(5):
                    for b in range(1, 5):
                        if cell.contains((r, a, b)):
                            for parity, f in enumerate((a, b)):
                                self.assertTrue(all(lo <= x <= hi for x, (lo, hi) in
                                                    zip((r-2*f, f), seen[parity][1])))

    def test_scope_timeout_and_nonfinite_fail_closed(self):
        bank, cell = self.bank(), c.Cell('interior', ((9, 9), (2, 2), (3, 3)))
        with self.assertRaises(ValueError):
            c.Evaluator(bank)
        with self.assertRaises(ValueError):
            c.fixed_bound(bank, c.Cell('interior', ((3, 3), (0, 0), (1, 1)), pairs=4),
                          (0, 1), 1, c.base.Budget())
        expired = SimpleNamespace(check=lambda: (_ for _ in ()).throw(TimeoutError()))
        with self.assertRaises(TimeoutError):
            c.fixed_bound(bank, cell, (0, 1), 1, expired)
        bank.scalar = lambda *args, **kwargs: arb('nan')
        with self.assertRaises(ArithmeticError):
            c.fixed_bound(bank, cell, (0, 1), 1, c.base.Budget())

    def test_representatives(self):
        for p in range(1, 8):
            for kind in ('interior', 'even_face', 'odd_face'):
                ranges = ((-1, 2*p), (-1, p), (-1, p)) if kind == 'interior' else ((-1, p), (-1, p))
                cell = c.Cell(kind, ranges, pairs=p)
                point = c.representative(cell)
                self.assertEqual(point is None, cell.integer_empty())
                if point is not None:
                    self.assertTrue(cell.contains(point))

    def test_public_receipt_rebounds_every_shortlisted_marker_pair(self):
        with ctx.workprec(256):
            bank = self.provider_stub()
            evaluator = c.Evaluator(bank)
            cell = c.Cell('interior', ((9, 9), (2, 2), (3, 3)))
            selected = evaluator.candidate_indices(cell, c.base.Budget())
            self.assertEqual(len(selected), 2)
            self.assertTrue(all(1 <= len(row) <= 4 for row in selected))
            expected_pairs = {(a, b) for a in selected[0] for b in selected[1]}
            with patch.object(c, 'fixed_bound', wraps=c.fixed_bound) as bound:
                result = evaluator.receipt(cell, max_seconds=10)
            self.assertEqual({tuple(call.args[2]) for call in bound.call_args_list}, expected_pairs)
            self.assertEqual(result['marker_pairs_evaluated'], len(expected_pairs))
            self.assertTrue(result['complete_cell'])
            self.assertTrue(result['directed'])
            self.assertFalse(result['whole_code_certificate'])
            controls = [c.fixed_bound(bank, cell, pair, result['count']['support_count_upper'], c.base.Budget())
                        for pair in expected_pairs]
            self.assertEqual(F(result['log2_upper_rational']),
                             min(F(row['log2_upper_rational']) for row in controls))
            self.assertEqual(result['witness']['selected_markers'],
                             [[str(value) for value in bank.markers[index]]
                              for index in result['marker_indices']])

    def test_public_empty_receipt_skips_markers_and_transfer(self):
        with ctx.workprec(256):
            evaluator = c.Evaluator(self.provider_stub())
            cell = c.Cell('even_face', ((1, 1), (2, 2)))
            with patch.object(evaluator, 'candidate_indices') as candidates, \
                 patch.object(c, 'fixed_bound') as bound:
                result = evaluator.receipt(cell, max_seconds=10)
            candidates.assert_not_called()
            bound.assert_not_called()
            self.assertTrue(result['zero'])
            self.assertTrue(result['integer_empty'])
            self.assertTrue(result['complete_cell'])
            self.assertTrue(result['directed'])
            self.assertEqual(result['marker_pairs_evaluated'], 0)
            self.assertIsNone(result['log2_upper_rational'])

    def test_public_timeout_is_pending_and_pin_change_is_fatal(self):
        with ctx.workprec(256):
            evaluator = c.Evaluator(self.provider_stub())
            cell = c.Cell('interior', ((9, 9), (2, 2), (3, 3)))
            with patch.object(c.counts, 'count_diagnostics', side_effect=TimeoutError('mock before count')):
                before_count = evaluator.receipt(cell, max_seconds=10)
            calls, actual_bound = [], c.fixed_bound

            def interrupted(*args, **kwargs):
                calls.append(args[2])
                if len(calls) == 2:
                    raise TimeoutError('mock after completed candidate')
                return actual_bound(*args, **kwargs)

            with patch.object(c, 'fixed_bound', side_effect=interrupted):
                after_candidate = evaluator.receipt(cell, marker_pairs=((0, 0), (0, 1)), max_seconds=10)
            self.assertEqual(calls, [(0, 0), (0, 1)])
            for result in (before_count, after_candidate):
                self.assertEqual(result['status'], 'pending')
                self.assertFalse(result['directed'])
                self.assertFalse(result['complete_cell'])
                self.assertFalse(result['meets_cell_budget'])
                self.assertIsNone(result['log2_upper_rational'])
                self.assertEqual(result['failure']['error_type'], 'TimeoutError')
            with patch.object(c.profiles, 'check_source_pins', side_effect=RuntimeError('mock changed source')):
                with self.assertRaises(RuntimeError):
                    evaluator.receipt(cell, marker_pairs=((0, 0),), max_seconds=10)


if __name__ == '__main__':
    unittest.main()

"""Mocked bounded orchestration; no regional bank or quantitative gate."""
from contextlib import ExitStack
import unittest
from unittest.mock import patch

import hd_two_sum_gf16_group_seen_control as target


def fixture():
    return dict(banks=[dict(active_groups=r, log2_slopes=[0., 1., 2., 3.],
        source_strata=[dict(full_max=1, count=100+r)]) for r, _, _ in target.reference.low.ANCHORS],
        witnesses=[dict(active_groups=r, source_full_pair_max=1,
            template=dict(name=name), failure_log2_screen=float(r+i))
                   for r, _, _ in target.reference.low.ANCHORS
                   for i, name in enumerate(target.reference.TEMPLATES)])


def templates(bank):
    return [(dict(name=name, literal='unchanged', ordinal=i), dict(r=bank['active_groups'], i=i))
            for i, name in enumerate(target.reference.TEMPLATES)]


class GroupSeenControlTests(unittest.TestCase):
    def mocks(self, stack, source=None, evaluate=None):
        stack.enter_context(patch.object(target, 'checked_source', return_value=source or fixture()))
        stack.enter_context(patch.object(target.reference, 'restricted_curves', return_value='curves'))
        stack.enter_context(patch.object(target.reference, 'fixed_templates', side_effect=templates))
        stack.enter_context(patch.object(target.reference.cover.emptycut, 'prepared_envelope',
            side_effect=lambda curves, slopes, weights: weights))
        stack.enter_context(patch.object(target.reference.cover.emptycut, 'evaluate_prepared',
            side_effect=lambda prepared, r, mark, alpha: dict(failure_log2_screen=float(r+prepared['i']))))
        stack.enter_context(patch.object(target.reference, 'restrict_count', side_effect=lambda row, r, count: row))
        if evaluate is None:
            def evaluate(prepared, weights, r, mark, alpha, source_count, check_time):
                check_time()
                self.assertEqual(prepared, weights)
                self.assertEqual(source_count, 100+r)
                self.assertEqual(alpha, target.reference.low.ALPHA)
                self.assertEqual(mark, next(z for n, _, z in target.reference.low.ANCHORS if n == r))
                value = r-2 if weights['i'] == 0 else r+2
                return dict(active_groups=r, seen_all_failure_log2_screen=float(value),
                    aggregate_9r_failure_log2_screen=float(value+1), failure_log2_screen=float(value))
        return stack.enter_context(patch.object(target.seen, 'evaluate', side_effect=evaluate))

    def test_checked_source_loads_only_pinned_files(self):
        data = target.checked_source()
        self.assertEqual(len(data['banks']), 3)
        self.assertEqual(target.reference.FULL_MAX, 1)
        self.assertEqual(target.reference.TEMPLATES, ('old-all-zero', 'joint-empty-floor1'))

    def test_six_controls_preserve_scopes_counts_and_old_minimum(self):
        with ExitStack() as stack:
            evaluate = self.mocks(stack)
            result = target.run()
        self.assertEqual(evaluate.call_count, 6)
        self.assertEqual(len(result['witnesses']), 6)
        self.assertTrue(result['complete_requested_controls'])
        self.assertEqual(result['regional_banks_rebuilt'], 0)
        self.assertFalse(result['all_source_supports'])
        self.assertFalse(result['constraints_imposed_jointly'])
        self.assertFalse(result['individual_nine_column_floors_enforced'])
        for row in result['witnesses']:
            self.assertEqual(row['source_full_pair_max'], 1)
            self.assertEqual(row['template']['literal'], 'unchanged')
            self.assertEqual(row['best_with_old_failure_log2_screen'],
                min(row['previous_emptycut_failure_log2_screen'], row['failure_log2_screen']))
        self.assertTrue(all(row['template']['ordinal'] == 0 for row in result['selected']))

    def test_failed_baseline_blocks_new_control(self):
        source = fixture()
        source['witnesses'][0]['failure_log2_screen'] += 1
        with ExitStack() as stack:
            evaluate = self.mocks(stack, source=source)
            with self.assertRaises(ArithmeticError):
                target.run()
            evaluate.assert_not_called()

    def test_missing_baseline_blocks_new_control(self):
        source = fixture()
        source['witnesses'].pop(0)
        with ExitStack() as stack:
            evaluate = self.mocks(stack, source=source)
            with self.assertRaises(ArithmeticError):
                target.run()
            evaluate.assert_not_called()

    def test_timeout_retains_partial_rows_and_in_progress_scope(self):
        calls = 0
        def evaluate(prepared, weights, r, mark, alpha, source_count, check_time):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise TimeoutError('bounded')
            return dict(active_groups=r, seen_all_failure_log2_screen=-1.,
                        aggregate_9r_failure_log2_screen=-2., failure_log2_screen=-2.)
        with ExitStack() as stack:
            self.mocks(stack, evaluate=evaluate)
            result = target.run()
        self.assertFalse(result['complete_requested_controls'])
        self.assertEqual(result['status'], 'compute-budget-exceeded')
        self.assertEqual(len(result['witnesses']), 1)
        self.assertEqual(result['in_progress']['template'], 'joint-empty-floor1')
        self.assertNotIn('selected', result)

    def test_invalid_budget_and_initial_timeout_do_not_load(self):
        for budget in (31, True, float('nan')):
            with patch.object(target, 'checked_source') as load, self.assertRaises(ValueError):
                target.run(budget)
            load.assert_not_called()
        with patch.object(target, 'checked_source') as load:
            result = target.run(1e-12)
        load.assert_not_called()
        self.assertFalse(result['complete_requested_controls'])


if __name__ == '__main__':
    unittest.main()

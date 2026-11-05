"""Exact small source-geometry tests and mocked orchestration; no gate."""
from itertools import combinations
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_source_strata as target


def matching_counts(size, support):
    support = set(support)
    return tuple(sum(i in support and (i+1) % size in support for i in range(offset, size, 2))
                 for offset in (0, 1))


class SourceStrataTests(unittest.TestCase):
    def test_exact_indexed_partition_and_every_subset_constraint(self):
        for size, groups in ((6, 2), (8, 3), (8, 4)):
            counts = {f: 0 for f in range(groups//2+1)}
            for support in combinations(range(size), groups):
                full_max = max(matching_counts(size, support))
                counts[full_max] += 1
                for weight in range(groups+1):
                    for subset in combinations(support, weight):
                        for full in matching_counts(size, subset):
                            self.assertLessEqual(full, min(full_max, weight//2))
            self.assertEqual(target.checked_source_counts(size, groups), counts)
            self.assertEqual(sum(counts.values()), comb(size, groups))

    def test_restricted_curves_use_an_upper_bound_not_exact_f(self):
        maximum = 6
        regions = {(m, f): np.array([[m-2*f, f-m], [3*f-m, -m-f]], dtype=float)
                   for m in range(maximum+1) for f in range(m//2+1)}
        previous = None
        for full_max in range(maximum//2+1):
            curves = target.restricted_curves(regions, maximum, full_max)
            expected = np.asarray([np.max([regions[m, f] for f in range(min(full_max, m//2)+1)], axis=0).ravel()
                                   for m in range(maximum+1)]).T
            np.testing.assert_array_equal(curves, expected)
            if previous is not None:
                self.assertTrue(np.all(previous <= curves))
            previous = curves
        np.testing.assert_array_equal(previous, target.state.regional_curves(regions, maximum))
        self.assertEqual(target.restricted_curves(regions, maximum, 3)[0, 6], 6.)

    def test_source_count_is_outside_alpha_and_charged_once(self):
        prepared = dict(regions=16, typed=np.zeros((1, 16)), classes=((0, 0, (16, 0, 0, 0), 1),),
            path_a=np.array([0.]), slopes=np.zeros(4), intercepts=np.zeros(4), slack=np.zeros((4, 3)),
            log_c=-2., allowed_counts=np.array([[16, 0, 0, 0]]))
        identity = np.array([[0., -np.inf], [-np.inf, 0.]])
        with patch.object(target.state, 'kernel_from_types', return_value=identity):
            first = target.evaluate_stratum(prepared, 18, 1, 3, target.F(691, 700))
            second = target.evaluate_stratum(prepared, 18, 1, 21, target.F(691, 700))
        self.assertAlmostEqual(second['failure_log2_screen']-first['failure_log2_screen'], log2(7))
        expected = .2*(-target.low.fractional.LIMIT*log2(691/700)-log2(255))+log2(3)
        self.assertAlmostEqual(first['failure_log2_screen'], expected)
        self.assertTrue(first['source_count_outside_alpha'])
        self.assertTrue(first['fixed_source_across_all_columns'])
        self.assertFalse(first['exact_total_arity_restriction'])
        self.assertFalse(first['directed'])

    def test_frozen_reference_loads_without_provider(self):
        with patch.object(target.paired, 'paired_nz_regions') as provider:
            data = target.load_reference()
        self.assertEqual(len(data['banks']), 3)
        self.assertEqual(len(data['witnesses']), 45)
        provider.assert_not_called()

    def test_three_banks_forty_strata_and_eight_unchanged_templates(self):
        reference = target.load_reference()
        prepared_calls, evaluated = [], []
        def provider(mark, alpha, maximum):
            self.assertEqual((alpha, maximum), (target.low.ALPHA, 32))
            index = [x[2] for x in target.low.ANCHORS].index(mark)
            curves = np.asarray([[(-np.inf if x is None else x) for x in row]
                                 for row in reference['banks'][index]['curves']])
            return {(m, f): curves[:, m].reshape(2, 2) for m in range(33) for f in range(m//2+1)}
        def counts(size, groups):
            self.assertEqual(size, 8192)
            return {f: comb(size, groups)-groups//2 if f == 0 else 1 for f in range(groups//2+1)}
        def templates(groups, mark):
            return [(dict(name=f'template-{i}'), (groups, i)) for i in range(8)]
        def prepare(curves, slopes, weights):
            anchor_index = [x[0] for x in target.low.ANCHORS].index(weights[0])
            self.assertEqual(slopes.tolist(), reference['witnesses'][15*anchor_index]['log2_slopes'])
            prepared_calls.append(weights)
            return weights
        def evaluate(prepared, groups, full_max, count, mark):
            self.assertEqual(prepared[0], groups)
            evaluated.append((groups, full_max, count, prepared[1]))
            return dict(active_groups=groups, source_full_pair_max=full_max, source_supports=count,
                        failure_log2_screen=log2(count)-prepared[1])
        with patch.object(target.paired, 'paired_nz_regions', side_effect=provider) as build, \
                patch.object(target, 'checked_source_counts', side_effect=counts), \
                patch.object(target.control, 'type_templates', side_effect=templates) as get_templates, \
                patch.object(target.cut, 'prepared_envelope', side_effect=prepare), \
                patch.object(target, 'evaluate_stratum', side_effect=evaluate):
            result = target.run()
        self.assertEqual(build.call_count, 3)
        self.assertEqual(get_templates.call_count, 3)
        self.assertEqual(len(prepared_calls), 320)
        self.assertEqual(len(evaluated), 320)
        self.assertEqual(len(result['strata']), 40)
        self.assertTrue(result['complete_requested_anchors'])
        self.assertEqual(result['fitted_parameters'], 0)
        self.assertTrue(all(row['template_index'] == 7 for row in result['strata']))
        for row, (groups, _, _) in zip(result['selected'], target.low.ANCHORS):
            self.assertAlmostEqual(row['failure_log2_screen'], log2(comb(8192, groups))-7)
            self.assertEqual(row['exact_source_support_count'], comb(8192, groups))
            self.assertTrue(row['complete_source_partition'])
            self.assertFalse(row['extra_template_factor'])
        self.assertFalse(result['directed'])
        self.assertFalse(result['continuous_weight_band_covered'])

    def test_incomplete_partition_fails_closed(self):
        with patch.object(target.geometry, 'matching_max_counts', return_value={0: 1}), self.assertRaises(ArithmeticError):
            target.checked_source_counts(8, 4)

    def test_invalid_budget_and_timeout_do_not_build_banks(self):
        for budget in (31, True, float('nan')):
            with patch.object(target.paired, 'paired_nz_regions') as provider, self.assertRaises(ValueError):
                target.run(budget)
            provider.assert_not_called()
        with patch.object(target, 'load_reference') as load:
            result = target.run(1e-12)
        load.assert_not_called()
        self.assertEqual(result['status'], 'compute-budget-exceeded')
        self.assertFalse(result['complete_requested_anchors'])


if __name__ == '__main__':
    unittest.main()

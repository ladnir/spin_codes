"""Tiny structural and Arb checks; no state census or large proof computation."""
from fractions import Fraction as Q
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx
import packet_return_ablation as diagnostic


def family():
    return [arb_mat([[1, 0], [0, arb(1)/2]]),
            arb_mat([[arb(1)/8, arb(1)/4], [arb(1)/16, arb(1)/2]]),
            arb_mat([[arb(1)/32, arb(1)/8], [arb(1)/8, arb(1)/4]])]


class ReturnAblationTests(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 192

    def tearDown(self):
        ctx.prec = self.precision

    def test_deletes_only_declared_physical_entries_without_mutation(self):
        original = family()
        copies = [arb_mat(m) for m in original]
        for variant in diagnostic.DEFAULT_VARIANTS:
            changed = diagnostic.alter_physical(original, variant)
            for j, (old, new) in enumerate(zip(original, changed)):
                if variant == 'zero_only':
                    self.assertEqual(new, arb_mat([[old[0, 0]]]))
                    continue
                for i in range(2):
                    for k in range(2):
                        deleted = ((i == k == 0 and
                                    ((variant == 'no_occupied_zero_stay' and j > 0) or
                                     (variant == 'no_pair_zero_stay' and j == 2))) or
                                   (variant == 'no_live_return' and i > 0 and k == 0))
                        self.assertEqual(new[i, k], 0 if deleted else old[i, k])
        self.assertEqual(original, copies)

    def test_macro_is_rebuilt_after_deletion_and_keeps_order(self):
        ops = family()
        changed = diagnostic.alter_physical(ops, 'no_live_return')
        actual = diagnostic.kernel.convolve(changed)
        # The occupancy-one split has equal probability in either half.
        expected = (changed[0]*changed[1]+changed[1]*changed[0])/2
        self.assertEqual(actual[1], expected)
        old_macro = diagnostic.kernel.convolve(ops)
        wrong = diagnostic.alter_physical(old_macro, 'no_live_return')
        self.assertNotEqual(actual[2], wrong[2])

    def test_baseline_matches_retained_composition(self):
        ops = family()
        result = diagnostic.evaluate_physical(ops, group_count=8, regions=3,
            occupancy=3, tilt=Q(1, 10), beta=Q(3, 2), threshold=2)
        macro = diagnostic.kernel.convolve(ops)
        regional = diagnostic.polynomial.placement_power(macro, epochs=2,
            windows=4, rounding=diagnostic.q1.rounded, maximum_groups=3)
        complete = diagnostic.uniform.regional_uniform(regional, 3)**3
        expected = diagnostic.kernel.up(56*diagnostic.kernel.aq(Q(3, 2))**3 *
            (diagnostic.kernel.aq(Q(1, 10))*2).exp() * sum(complete[0, j] for j in range(2)))
        actual = diagnostic._arb_endpoint(result['baseline']['diagnostic_mass_dyadic'])
        self.assertLess(abs(actual-expected), arb(2)**-150)
        for value in result.values():
            self.assertFalse(value['is_probability_bound'])
            self.assertFalse(value['whole_code_certificate'])
            self.assertLessEqual(diagnostic._arb_endpoint(value['diagnostic_mass_dyadic']), actual)

    def test_counting_beta_and_cutoff_factors_are_preserved(self):
        ops = [arb_mat([[1]]) for _ in range(3)]
        result = diagnostic.evaluate_physical(ops, group_count=8, regions=3,
            occupancy=3, tilt=Q(1, 5), beta=Q(7, 3), threshold=2, variants=('baseline',))
        expected = 56*diagnostic.kernel.aq(Q(7, 3))**3*(diagnostic.kernel.aq(Q(2, 5))).exp()
        actual = diagnostic._arb_endpoint(result['baseline']['diagnostic_mass_dyadic'])
        self.assertLess(abs(actual-expected), arb(2)**-140)

    def test_zero_projection_agrees_with_full_dimension_projection(self):
        ops = family()
        scalar = diagnostic.alter_physical(ops, 'zero_only')
        full = [arb_mat([[m[0, 0], 0], [0, 0]]) for m in ops]
        a = diagnostic.evaluate_physical(scalar, group_count=4, regions=2,
            occupancy=2, tilt='1/8', beta=1, threshold=0, variants=('baseline',))
        b = diagnostic.evaluate_physical(full, group_count=4, regions=2,
            occupancy=2, tilt='1/8', beta=1, threshold=0, variants=('baseline',))
        self.assertEqual(a['baseline']['diagnostic_mass_dyadic'], b['baseline']['diagnostic_mass_dyadic'])

    def test_arb_does_not_underflow_tiny_positive_path(self):
        tiny = arb(2)**-2000
        ops = [arb_mat([[tiny]]) for _ in range(2)]
        result = diagnostic.evaluate_physical(ops, group_count=4, regions=2,
            occupancy=1, tilt='1/8', beta=1, threshold=0, variants=('baseline',))['baseline']
        value = diagnostic._arb_endpoint(result['diagnostic_mass_dyadic'])
        self.assertTrue(value > 0)
        self.assertFalse(result['exact_zero_in_modified_transfer'])
        self.assertLess(abs(value/(arb(2)**-15998)-1), arb(2)**-140)

    def test_exact_zero_has_no_fake_finite_log_or_certificate(self):
        ops = [arb_mat([[0]]) for _ in range(2)]
        result = diagnostic.evaluate_physical(ops, group_count=2, regions=2,
            occupancy=1, tilt='1/8', beta=1, threshold=0, variants=('baseline',))['baseline']
        self.assertEqual(result['diagnostic_mass_dyadic'], [0, 0])
        self.assertTrue(result['exact_zero_in_modified_transfer'])
        self.assertIsNone(result['negative_log2_diagnostic_mass'])
        self.assertIsNone(result['ratio_to_baseline_endpoint'])

    def test_run_builds_and_selects_once_and_emits_only_diagnostic_schema(self):
        data, record, sources = {'toy': True}, {'s': 22}, {'toy-source': 'toy-hash'}
        adapter = diagnostic.ladder.adapter_for(22)
        ops = [arb_mat([[arb(1)/2, arb(1)/4], [arb(1)/16, arb(1)/2]]) for _ in range(17)]
        with patch.object(diagnostic.ladder, '_prepare', return_value=(data, record, sources)), \
             patch.object(adapter, 'authenticate'), \
             patch.object(diagnostic.kernel, 'authenticate', return_value={}), \
             patch.object(diagnostic.kernel.sparse_kernel, 'outward_at_z', return_value=ops) as build, \
             patch.object(diagnostic.kernel.kernel_birth_density, 'refine_local', return_value=ops) as select, \
             patch.object(diagnostic.q1, 'source_snapshot', return_value=sources):
            result = diagnostic.run(data, record, K=4096, occupancy=3, tilt='1/8')
        self.assertEqual(build.call_count, 1)
        self.assertEqual(select.call_count, 1)
        self.assertEqual(result['schema'], diagnostic.SCHEMA)
        self.assertTrue(result['diagnostic_only'])
        self.assertFalse(result['is_probability_bound'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['all_occupancies_covered'])
        self.assertNotIn('occupancy_uppers', result)
        self.assertEqual(result['threshold'], 819)
        self.assertEqual(result['map_record'], record)

    def test_rejects_bad_variant_geometry_and_unprepared_run(self):
        for variants in ((), ('no_live_return',), ('baseline', 'baseline'), 'baseline'):
            with self.assertRaises(ValueError):
                diagnostic.evaluate_physical(family(), group_count=4, regions=1,
                    occupancy=1, tilt='1/8', beta=1, threshold=0, variants=variants)
        for updates in (dict(group_count=3), dict(regions=True), dict(occupancy=5),
                        dict(threshold=-1), dict(beta=0), dict(tilt=0)):
            settings = dict(group_count=4, regions=1, occupancy=1, tilt='1/8', beta=1, threshold=0)
            settings.update(updates)
            with self.assertRaises(ValueError):
                diagnostic.evaluate_physical(family(), **settings)
        with self.assertRaises(ValueError):
            diagnostic.alter_physical([arb_mat([[-1]])]*2, 'baseline')
        with self.assertRaises(ValueError):
            diagnostic.run(None, None, K=4096, occupancy=3, tilt='1/8')


if __name__ == '__main__':
    unittest.main()

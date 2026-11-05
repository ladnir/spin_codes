"""Lightweight bounded-search, shared-bank and complete-partition controls."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import numpy as np

import hd_xor_sparse_cohort_refine as gate
from hd_xor_geometry_strata import category_resources


class SparseCohortRefineTests(unittest.TestCase):
    def test_fit_uses_one_bank_and_two_caps(self):
        calls = []
        class Caps:
            def __init__(self, *args):
                calls.append('bank')
            def __call__(self, mu):
                calls.append(mu)
                return {'mu': mu}
        def row(n1, bounds, geometries, count, caps, label, z, p):
            self.assertEqual(caps['mu'], bounds[1])
            return dict(cohort_counts=[n1, 17-n1, 0], failure_log2_screen=-n1*10)
        with patch.object(gate, 'pair_regions_screen', return_value=None) as raw, \
                patch.object(gate, 'ordered_region_tables', return_value=None), \
                patch.object(gate, 'MultiplicityCapsFloat', Caps), \
                patch.object(gate, 'class_row', side_effect=row):
            result = gate.evaluate_joint(*gate.INITIAL)
        self.assertEqual(raw.call_count, 1)
        self.assertEqual(calls, ['bank', 14, 15])
        self.assertEqual([row['cohort_counts'][0] for row in result['entries']], [2, 3])
        self.assertAlmostEqual(result['objective_log2_screen'], float(np.logaddexp2(-20, -30)))

    def test_fit_hard_budget_and_gain_gate(self):
        def fake(z, p):
            value = ((float(z)-.994)/.0001)**2+(float(p)+65)**2
            return dict(mark=str(z), log2_potential=str(p), objective_log2_screen=value, entries=[])
        with patch.object(gate, 'evaluate_joint', side_effect=fake):
            fit = gate.fit_joint()
        self.assertLessEqual(fit['evaluations'], 32)
        self.assertEqual(fit['evaluations'], len(fit['proposals']))
        self.assertGreaterEqual(fit['gain_bits'], 0)
        self.assertEqual(fit['full_screen_authorized_by_gain'], fit['gain_bits'] >= 10)
        self.assertEqual(F(fit['proposals'][0]['mark']), gate.INITIAL[0])
        self.assertEqual(F(fit['proposals'][0]['log2_potential']), gate.INITIAL[1])

    def test_canonical_cap_identity(self):
        for n2, bins in gate.partitions().items():
            for _, high, _ in bins:
                for double in range(n2+1):
                    for single in range(18-double):
                        for cancelled in range(double+1):
                            self.assertEqual(category_resources(single, double, cancelled, max(1, high)),
                                category_resources(single, double, cancelled, gate.canonical_mu(high)))

    def test_summary_requires_complete_cells_and_separates_directed_subterm(self):
        rows = [dict(source_matching_counts=list(source), cohort_counts=[n1, 17-n1, 0],
            mu_range=list(bounds), failure_log2_screen=-200-n1) for source, n1, bounds in gate.expected_keys()]
        result = gate.summarize_rows(rows)
        self.assertEqual(len(rows), 2673)
        self.assertEqual(len(result['cohorts']), 18)
        self.assertAlmostEqual(result['failure_log2_screen'], float(np.logaddexp2.reduce([
            row['failure_log2_screen'] for row in rows])))
        expected = float(np.logaddexp2(float(gate.HOMOGENEOUS_DIRECTED_CEILING), result['nonhomogeneous_log2_screen']))
        self.assertAlmostEqual(result['preserved_homogeneous_mixed_log2_screen'], expected)
        for invalid in (rows[:-1], rows+[rows[0]]):
            with self.assertRaises(ValueError):
                gate.summarize_rows(invalid)

    def test_fixed_screen_rejects_incomplete_controls_before_work(self):
        with patch.object(gate, 'pair_regions_screen') as raw:
            with self.assertRaises(ValueError):
                gate.fixed_bank_screen(dict(mark='249/250', log2_potential='-64'), {})
            raw.assert_not_called()

    def test_current_pinned_controls_and_scope(self):
        from pathlib import Path
        directory = Path(gate.__file__).parent
        if not all((directory/name).exists() for name in gate.INPUT_PINS):
            self.skipTest('Frozen proposal/receipt artifacts are optional local fixtures')
        controls = gate.load_controls(directory)
        self.assertEqual(set(controls), gate.expected_keys())
        evidence = {row['input_evidence'] for row in controls.values()}
        self.assertIn('directed-homogeneous-cell', evidence)
        self.assertIn('sibling-FLOAT-schedule', evidence)
        self.assertIn('frozen-two-bank-FLOAT', evidence)

    def test_saved_proposal_recombination_without_operators(self):
        import json
        from pathlib import Path
        directory = Path(gate.__file__).parent
        fit_path = directory/'hd_xor_sparse_cohort_refine32_0_7_fit.json'
        screen_path = directory/'hd_xor_sparse_cohort_refine32_all_sources_float.json'
        if not fit_path.exists() or not screen_path.exists():
            self.skipTest('Bounded run artifacts are optional local fixtures')
        fit = json.loads(fit_path.read_text())['fit']
        screen = json.loads(screen_path.read_text())
        self.assertEqual(fit['evaluations'], 32)
        self.assertEqual(fit['best'], min(fit['proposals'], key=lambda row: row['objective_log2_screen']))
        self.assertGreaterEqual(fit['gain_bits'], 10)
        self.assertFalse(screen['directed'])
        self.assertEqual(screen['combined_evidence'], 'MIXED-PRELIMINARY')
        controls = gate.load_controls(directory)
        new = {gate.row_key(row): row for row in screen['raw_new_bank']['entries']}
        selected = {gate.row_key(row): row for row in screen['combined']['entries']}
        self.assertEqual(set(new), gate.expected_keys())
        self.assertEqual(set(selected), gate.expected_keys())
        for key, row in new.items():
            self.assertEqual((row['mark'], row['log2_potential']), (fit['best']['mark'], fit['best']['log2_potential']))
            self.assertEqual(selected[key]['failure_log2_screen'], min(row['failure_log2_screen'], controls[key]['failure_log2_screen']))
            self.assertEqual(row['source_supports'], gate.source_counts()[key[0]])
            self.assertEqual(row['sparse_support_floor'], 17*(34-key[1]))
            self.assertEqual(row['sparse_support_normalizer_exponent'], 16*(34-key[1]))
        total = float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in selected.values()]))
        self.assertAlmostEqual(total, screen['combined']['failure_log2_screen'], places=10)
        self.assertLess(screen['combined']['preserved_homogeneous_mixed_log2_screen'], -50)


if __name__ == '__main__':
    unittest.main()

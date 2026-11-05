from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
from unittest.mock import patch
import unittest

from flint import arb, ctx

import hd_lambda_renyi_zero_control as control


class ZeroControlTests(unittest.TestCase):
    def test_bank_pins_scope_and_precision_separation(self):
        for precision in (256, 384):
            path = control.DIRECTORY/f'hd_lambda_renyi_zero_reference_final_{precision}.json'
            data, (envelopes, count, _) = control.load_bank(path, precision)
            self.assertEqual(len(envelopes), 19)
            self.assertEqual(count, sum(control.A.profile_counts(18).values()))
            with self.assertRaises(ValueError):
                control.load_bank(path, 640-precision)
            for key, replacement in (('alpha', '1/4'), ('source_count', str(count+1)),
                                     ('precision_bits', 640-precision), ('whole_column_flags', False)):
                bad = deepcopy(data)
                bad[key] = replacement
                with self.assertRaises(ValueError):
                    control.validate_bank(bad, precision)

    def test_reference_selection_never_uses_uniform_for_nonzero_eta(self):
        with ctx.workprec(256), patch.object(control.reference, 'reference_log', return_value=arb(7)) as evaluate:
            nonzero = control.choose_reference('bank', control.ETAS[0], 123, F(-100))
            zero = control.choose_reference('bank', F(0), 123, F(-100))
        self.assertEqual(nonzero['selected_reference']['reference_source_log2_upper'], '7')
        self.assertEqual(zero['selected_reference']['reference_source_log2_upper'], '-100')
        self.assertEqual(len(nonzero['all_reference_controls']), 9)
        self.assertEqual(len(zero['all_reference_controls']), 10)
        self.assertEqual(evaluate.call_count, 18)
        for call in evaluate.call_args_list:
            self.assertEqual(call.kwargs['minimum_flags'], 162)
            self.assertEqual(call.kwargs['source_count'], 123)

    def test_weights_match_direct_new_reference_atoms(self):
        prepared = {}
        for a in control.ARITIES:
            old_z3 = control.reference.activity.zero_atom(a, 3)
            old_z4 = control.reference.activity.zero_atom(a, 4)
            # Valid synthetic moment upper bounds with two actual column classes.
            one = (max(F(1), old_z3)+1-old_z3, max(F(1), old_z4)+2-old_z4)
            two = (one[0]+1, one[1]+3)
            prepared[a] = Counter({(one,)*7+(two,)*9: 2})
        with ctx.workprec(256):
            for p in (3, 4):
                for eta in (F(0), control.ETAS[0]):
                    weights, _ = control.weights_from_signatures(prepared, p, eta)
                    for a, weight in zip(control.ARITIES, weights):
                        signature = next(iter(prepared[a]))
                        old_z = control.reference.activity.zero_atom(a, p)
                        atoms = [control.reference.reference_moment_atoms(a, p, pair[p-3]+old_z, eta)
                                 for pair in signature]
                        tail = control.reference.activity.activity_tail(atoms[0][0], tuple(n for _, n in atoms), 9)
                        expected = 2*tail**control.A.uniform.rational(F(1, p))
                        self.assertTrue(weight.overlaps(expected))

    def test_normalization_beta_projective_and_reference_once(self):
        with ctx.workprec(256):
            total, p, reference = 216, 3, F(-777)
            coefficient = arb(255)*arb(255)**(8*total)/arb(256)**(16*total)
            got = control.total_bound(coefficient, total, p, reference)
            threshold = -control.A.uniform.rational(control.ALPHA)*734003*control.A.log2(
                control.A.uniform.rational(control.MARK))
            self.assertTrue(got.overlaps(threshold+control.A.uniform.rational(reference)))

    def test_completed_controls_keep_density_and_reference_eta_together(self):
        calls = []
        def choose(_prepared, eta, _source, _uniform):
            return dict(eta=str(eta), all_reference_controls=[], selected_reference=dict(
                flag_mark='1', kind='test', reference_source_log2_upper=str(-10 if eta else -20)))
        def weights(_prepared, p, eta):
            calls.append((p, eta))
            return [arb(1)]*6, []
        with ctx.workprec(256), patch.object(control, 'choose_reference', side_effect=choose), \
                patch.object(control, 'weights_from_signatures', side_effect=weights):
            data = control.combined_controls('moments', 'bank', 123, F(-20))
        self.assertEqual(len(data['controls']), 14)
        self.assertEqual(len(data['selected_totals']), 48)
        for row, (p, eta) in zip(data['controls'], calls):
            self.assertEqual(row['eta'], str(eta))
            self.assertEqual(row['moment_order'], p)
            self.assertEqual(row['reference_source_log2_upper'], str(-10 if eta else -20))

    def test_receipt_cached_only_and_immutable_inputs(self):
        precision = 256
        path = control.DIRECTORY/'hd_lambda_renyi_zero_reference_final_256.json'
        fake = dict(controls=[], reference_controls=[], selected_totals=[], union_log2_upper='1')
        with patch.object(control.interpolate, 'checked_inputs', return_value=({}, {})), \
                patch.object(control.interpolate, 'prepare_signatures', return_value='signatures'), \
                patch.object(control.interpolate, 'inputs_unchanged', return_value=True), \
                patch.object(control, 'combined_controls', return_value=fake), \
                patch.object(control.A.uniform, 'regional_directed', side_effect=AssertionError('No provider')):
            data = control.receipt(path, precision=precision)
        self.assertEqual(data['fresh_regional_banks'], 0)
        self.assertEqual(data['total_arity_range'], [198, 245])
        self.assertEqual(data['cached_bank_precision_bits'], precision)
        self.assertFalse(data['whole_code_certificate'])

    def test_guards_and_sources(self):
        self.assertEqual(control.dependencies()['hd_lambda_renyi_zero_reference.py'],
                         control.PINS['hd_lambda_renyi_zero_reference.py'])
        with self.assertRaises(ValueError):
            control.weights_from_signatures({}, 2, control.ETAS[0])
        with self.assertRaises(ValueError):
            control.choose_reference({}, F(1, 2), 1, F(0))
        with self.assertRaises(ValueError):
            control.weights_from_signatures({}, 3, control.ETAS[0])


if __name__ == '__main__':
    unittest.main()

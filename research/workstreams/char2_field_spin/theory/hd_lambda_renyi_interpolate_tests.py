"""Tiny interpolation/interface tests; no regional provider or inventory run."""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from unittest.mock import patch
import unittest

from flint import arb, ctx

import hd_lambda_renyi_interpolate as target


class InterpolationTests(unittest.TestCase):
    def test_pins_and_prescribed_orders(self):
        self.assertEqual(target.dependencies()['hd_lambda_renyi_activity.py'],
                         target.PINS['hd_lambda_renyi_activity.py'])
        self.assertEqual(target.ORDERS, tuple(F(k, 8) for k in range(25, 32)))
        self.assertEqual(target._orders((4, F(7, 2), 3)), (F(3), F(7, 2), F(4)))

    def test_interpolation_endpoints_and_zero_atom(self):
        n3, n4 = F(43), F(171)
        with ctx.workprec(256):
            for p, exact in ((3, n3), (4, n4)):
                self.assertTrue(target.nonzero_interpolation(n3, n4, p).contains(
                    target.A.uniform.rational(exact)))
                for arity in (11, 12, 16):
                    self.assertTrue(target.zero_contribution(arity, p).contains(
                        target.A.uniform.rational(target.activity.zero_atom(arity, p))))

    def test_nonzero_submeasure_logconvexity(self):
        # A subprobability measure with masses1/5,2/5 at likelihoods1,4.
        n3, n4, exact_half_order = F(129, 5), F(513, 5), F(257, 5)
        with ctx.workprec(256):
            upper = target.nonzero_interpolation(n3, n4, F(7, 2))
            self.assertTrue(upper > target.A.uniform.rational(exact_half_order))
            self.assertTrue(upper**2 == target.A.uniform.rational(n3*n4)
                            or (upper**2).contains(target.A.uniform.rational(n3*n4)))

    def test_positive_tail_matches_frozen_endpoint(self):
        z, atoms = F(2, 3), (F(1, 7), F(3, 5), F(4, 3), F(7, 2))
        with ctx.workprec(256):
            for floor in range(5):
                actual = target.positive_tail(target.A.uniform.rational(z),
                    tuple(target.A.uniform.rational(x) for x in atoms), floor)
                expected = target.activity.activity_tail(z, atoms, floor)
                self.assertTrue(actual.overlaps(expected))

    def test_original_label_multiplicities_and_paired_signatures(self):
        moments = {}
        for b in range(3):
            a = 16-b
            for complement in combinations(range(16), b):
                # Distinct deterministic labels exercise pairing, not just counts.
                n3 = F(1+sum(complement))
                n4 = F(2+sum(i*i for i in complement))
                moments[complement] = {3: target.activity.zero_atom(a, 3)+n3,
                                       4: target.activity.zero_atom(a, 4)+n4}
        prepared = target.prepare_signatures(moments, (14, 15, 16))
        schedules = [target.binding.column_permutations('input-inverse-output-xor', j)[0]
                     for j in range(16)]
        for a, signatures in prepared.items():
            self.assertEqual(sum(signatures.values()), comb(16, a))
            for complement in combinations(range(16), 16-a):
                pairs = []
                for before in schedules:
                    mapped = tuple(sorted(before[i] for i in complement))
                    pairs.append((F(1+sum(mapped)), F(2+sum(i*i for i in mapped))))
                self.assertIn(tuple(sorted(pairs)), signatures)
        with ctx.workprec(256):
            weights, diagnostic = target.scheduled_weights(prepared, 3, minimum=0)
            for a, weight, row in zip((14, 15, 16), weights, diagnostic):
                reference, _ = target.activity.scheduled_weights(moments, 3, minimum=0, arities=(a,))
                self.assertTrue(weight.overlaps(reference[0]))
                self.assertEqual(row['labelled_rowsets'], comb(16, a))

    def test_receipt_uses_one_common_order_and_normalization(self):
        # Isolated controller test: neither inventory validation nor provider runs.
        profile_terms = [(i, 0, 1, arb(1)) for i in range(99)]
        with patch.object(target, 'dependencies', return_value={'frozen': 'yes'}), \
                patch.object(target, 'checked_inputs', return_value=({}, {})), \
                patch.object(target, 'inputs_unchanged', return_value=True), \
                patch.object(target, 'prepare_signatures', return_value={}), \
                patch.object(target, 'scheduled_weights', return_value=([arb(1)]*6, [])), \
                patch.object(target.A.uniform, 'regional_directed', return_value={}) as provider, \
                patch.object(target.A, 'source_profile_terms', return_value=profile_terms):
            result = target.receipt(precision=256, orders=(F(7, 2),))
        provider.assert_called_once_with(target.MARK, target.ALPHA, 18, 60)
        self.assertEqual(result['controls'][0]['outer_theta'], '7/25')
        self.assertEqual(result['total_arity_range'], [198, 245])
        self.assertEqual(len(result['selected_totals']), 48)
        self.assertEqual(result['projective_normalizations'], 1)
        self.assertEqual(result['inventory_sha256'], target.INVENTORY_PINS)
        self.assertFalse(result['whole_code_certificate'])
        with ctx.workprec(256):
            # At total198 the coefficient is1; all beta powers remain outside it.
            expected = (target.A.uniform.rational(F(7, 25))*(198*(128-8*target.A.log2(arb(255)))
                        - target.A.log2(arb(255)))
                        - target.A.uniform.rational(target.ALPHA)*734003
                        * target.A.log2(target.A.uniform.rational(target.MARK))
                        + target.A.log2(arb(99)))
            self.assertEqual(result['selected_totals'][0]['failure_log2_upper'],
                             str(target.A.uniform.rational_ceiling(expected)))

    def test_guards_and_inventory_pin_failure(self):
        for orders in ((), (3, 3), (F(10, 3),), (True,), (3.5,)):
            with self.assertRaises(ValueError):
                target._orders(orders)
        for args in ((F(-1), F(1), F(7, 2)), (F(1), F(-1), 3), (F(1), F(1), 5)):
            with self.assertRaises(ValueError):
                target.nonzero_interpolation(*args)
        with patch.object(target.Path, 'read_bytes', return_value=b'{}'):
            with self.assertRaisesRegex(ValueError, 'Unreviewed exact inventory'):
                target.checked_inputs('not-used-small', 'not-used-middle')
        with self.assertRaises(ValueError):
            target.prepare_signatures({}, (11, 11))


if __name__ == '__main__':
    unittest.main()

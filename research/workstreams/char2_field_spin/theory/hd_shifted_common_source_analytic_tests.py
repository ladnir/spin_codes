"""Exact/light tests; no production event calculation or provider rebuild."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json
from math import log2
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import hd_shifted_common_source_analytic as gate
from hd_shifted_common_source_partition import good_probability, matrix_product, IDENTITY


class AnalyticGoodRouteTests(unittest.TestCase):
    def test_refresh_pointwise_atom_with_toy_superregular_matrix(self):
        # Every square minor of [[1,1],[1,2]] over F5 is nonzero.
        q, dimension = 5, 2
        diagonal_choices = list(product(range(1, q), repeat=dimension))
        denominator = len(diagonal_choices)**2
        for x in product(range(q), repeat=dimension):
            if not any(x):
                continue
            counts = {}
            for d1 in diagonal_choices:
                v = [d1[i]*x[i] % q for i in range(dimension)]
                w = ((v[0]+v[1]) % q, (v[0]+2*v[1]) % q)
                for d2 in diagonal_choices:
                    y = tuple(d2[i]*w[i] % q for i in range(dimension))
                    counts[y] = counts.get(y, 0)+1
            self.assertNotIn((0, 0), counts)
            self.assertEqual(sum(counts.values()), denominator)
            self.assertLessEqual(F(max(counts.values()), denominator), F(1, (q-1)**dimension))

    def test_return_uses_flat_nonzero_reference_not_a_stochastic_row(self):
        states = [x for x in product(range(5), repeat=2) if any(x)]
        for target in product(range(5), repeat=2):
            hits = sum(all((q+c) % 5 == 0 for q, c in zip(state, target)) for state in states)
            self.assertEqual(F(hits, len(states)), F(1, 24) if any(target) else 0)
        cap = gate.good_region_cap(F(1, 2), 4, q=5, state=2, emission_distance=3)
        self.assertEqual(cap['return_entry'], F(1, 24))
        self.assertGreater(sum(cap['packet'][2:]), 1)  # Deliberately not row-stochastic.

    def test_positive_powers_and_birth_return(self):
        cap = gate.good_region_cap(F(1, 2), 6, q=5, state=2, emission_distance=3)
        packet = cap['packet']
        beta, delta = cap['beta'], F(1, 16)
        recurrence = [F(1), F(0)]
        for k in range(2, 7):
            recurrence.append(beta*recurrence[-1]+delta*recurrence[-2])
        self.assertEqual(cap['power00'], recurrence[1:])
        for k in range(1, 7):
            exact = F(0)
            for middle in product((0, 1), repeat=k-1):
                path = (0,)+middle+(0,)
                weight = F(1)
                for before, after in zip(path, path[1:]):
                    weight *= packet[2*before+after]
                exact += weight
            self.assertEqual(exact, cap['power00'][k-1])
        self.assertEqual(cap['power00'][0], 0)
        self.assertEqual(cap['power00'][1], delta)
        self.assertEqual(cap['good00'], max(cap['power00']))

    def test_empty_contraction_preserves_order(self):
        cap = gate.good_region_cap(F(1, 2), 3, q=5, state=2, emission_distance=3)
        a = cap['packet']
        b = tuple(x*F(i+1, 4) for i, x in enumerate(a))
        c = tuple(x*F(4-i, 5) for i, x in enumerate(a))
        empty = (F(1), F(0), F(0), cap['empty_live_cap'])
        with_empty, without_empty = IDENTITY, IDENTITY
        for packet in (a, empty, b, empty, c):
            with_empty = matrix_product(with_empty, packet)
        for packet in (a, b, c):
            without_empty = matrix_product(without_empty, packet)
        self.assertNotEqual(matrix_product(a, b), matrix_product(b, a))
        self.assertTrue(all(x <= y for x, y in zip(with_empty, without_empty)))
        self.assertLessEqual(with_empty[0], cap['good00'])
        with self.assertRaises(ValueError):
            gate.good_region_cap(F(1), 3)

    def test_exact_tiny_probability_and_conditioning_power(self):
        p = good_probability(2, 0, block_length=2, packets=2, cap=1)
        self.assertEqual(p, F(13, 14))
        row = gate.conditional_comparison(F(3, 4), F(1, 16), -6., F(1, 2))
        # p*sqrt(1/16) + (1-p)*sqrt((1/64)/(1-p)) = 1/4.
        self.assertAlmostEqual(row['good_plus_bad_log2'], -2.)
        self.assertAlmostEqual(row['unconditional_jensen_log2'], -3.)
        self.assertEqual(row['minimum_with_unconditional_log2'], -3.)
        self.assertIsNone(gate.conditional_comparison(0, F(1, 16), -6., F(1, 2))['good_component_log2'])
        self.assertIsNone(gate.conditional_comparison(1, F(1, 16), -6., F(1, 2))['bad_component_log2'])
        self.assertAlmostEqual(gate.conditional_comparison(1, F(1, 16), -6., F(1, 2))['good_plus_bad_log2'], -2.)

    def test_saved_scope_and_byte_pins(self):
        value = json.loads((gate.DIRECTORY/gate.ARTIFACT).read_bytes())
        curve, line = gate.checked_artifact(value)
        self.assertEqual(curve, value['regional_envelope']['log2_curves'][9])
        self.assertEqual(line, 9*value['regional_envelope']['log2_slope'])
        self.assertEqual(gate.load_saved_curve(), (curve, line))
        for key, replacement in [('fresh_regional_alpha', .2), ('mark', '99/100'),
                                 ('active_groups', 17), ('whole_code_certificate', True)]:
            changed = deepcopy(value)
            changed[key] = replacement
            with self.assertRaises(ValueError):
                gate.checked_artifact(changed)
        changed = deepcopy(value)
        changed['regional_envelope']['log2_curves'][9] = None
        with self.assertRaises(ValueError):
            gate.checked_artifact(changed)
        with TemporaryDirectory() as directory:
            path = Path(directory)/gate.ARTIFACT
            path.write_bytes(b'{}')
            with self.assertRaises(ValueError):
                gate.load_saved_curve(path)

    def test_mocked_run_scope_and_exact_production_contraction(self):
        # Only mock probabilities are used here; no real r18 calculation.
        with patch.object(gate, 'good_probability', return_value=F(13, 14)) as probability:
            answer = gate.run()
        self.assertEqual([call.args for call in probability.call_args_list], [(18, 0), (16, 1)])
        self.assertEqual(answer['provider_rebuilds'], 0)
        self.assertFalse(answer['directed'])
        self.assertFalse(answer['whole_column_bound'])
        self.assertFalse(answer['whole_code_certificate'])
        self.assertFalse(answer['regional00_zero_claim'])
        self.assertTrue(answer['probability_outside_fractional_power'])
        self.assertLessEqual(F(answer['empty_live_cap']), 1)
        self.assertEqual(F(answer['beta'])/int(answer['flat_reference_dimension']), F(1, 255**16))
        self.assertTrue(all(row['actual_full_pairs_not_fixed'] for row in answer['cases']))
        json.dumps(answer, allow_nan=False)


if __name__ == '__main__':
    unittest.main()

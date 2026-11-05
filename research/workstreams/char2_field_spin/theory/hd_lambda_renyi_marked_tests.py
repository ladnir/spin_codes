from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from unittest.mock import patch
import unittest

import hd_lambda_renyi_marked as marked
import hd_lambda_arity5_relations as relations
from hd_lambda_outer import matrix
from hd_quaternary_outer_relations import multiply_gf4


def direct_inputs(rows):
    """Enumerate only tiny GF4 input spaces, independently of relation kernels."""
    arity = len(rows[0])
    result = []
    for x in product(range(4), repeat=arity):
        y = []
        for row in rows:
            value = 0
            for c, v in zip(row, x):
                value ^= multiply_gf4(c, v)
            y.append(value)
        result.append((sum(v != 0 for v in x), sum(bool(v) << i for i, v in enumerate(y))))
    return tuple(result)


def direct_replica_polynomial(inputs, arity, width, p):
    # Each ordered tuple contributes only when all its outputs have SAME
    # support. This does not use a grouped polynomial or its convolution.
    integer_terms = Counter()
    for replicas in product(inputs, repeat=p):
        first_mask = replicas[0][1]
        if all(mask == first_mask for _, mask in replicas):
            integer_terms[sum(weight for weight, _ in replicas), first_mask.bit_count()] += 1
    result = [F(0)]*(arity*p+1)
    for (degree, support_weight), count in integer_terms.items():
        result[degree] += F(4**(width*(p-1))*count, 3**((p-1)*support_weight))
    return tuple(result)


class MarkedRenyiTests(unittest.TestCase):
    def test_same_support_power_matches_direct_input_replicas(self):
        # Full output support has BOTH input weights1 and2. Masks01 and10
        # have equal weight but must NOT be combined before the p-th power.
        inputs = direct_inputs(((1, 1), (1, 2)))
        joint = dict(Counter(inputs))
        self.assertGreater(joint[1, 3], 0)
        self.assertGreater(joint[2, 3], 0)
        actual = marked.marked_polynomials_from_joint(joint, 2, width=2, q=4)
        for p in (2, 3):
            self.assertEqual(actual[p], direct_replica_polynomial(inputs, 2, 2, p))
        # Odd degrees require mixed input weights in the SAME support.
        self.assertGreater(actual[2][3], 0)
        self.assertGreater(actual[3][4], 0)
        one = direct_inputs(((1,), (2,)))
        one_polys = marked.marked_polynomials_from_joint(dict(Counter(one)), 1, width=2, q=4)
        self.assertEqual(one_polys[4], direct_replica_polynomial(one, 1, 2, 4))

    def test_raw_mass_identity_against_uniform_target_enumeration(self):
        for rows in (((1,), (2,)), ((1, 1), (1, 2))):
            arity, width = len(rows[0]), len(rows)
            inputs = direct_inputs(rows)
            joint = dict(Counter(inputs))
            counts = Counter(mask for _, mask in inputs)
            polynomials = marked.marked_polynomials_from_joint(joint, arity, width=width, q=4)
            moments = {p: F(0) for p in (2, 3, 4)}
            for y in product(range(4), repeat=width):
                mask = sum(bool(v) << i for i, v in enumerate(y))
                likelihood = F(4**(width-arity)*counts[mask], 3**mask.bit_count())
                for p in moments:
                    moments[p] += likelihood**p/4**width
            for p, coefficients in polynomials.items():
                self.assertEqual(sum(coefficients), 4**(arity*p)*moments[p])
                self.assertEqual(coefficients[0], 4**(width*(p-1)))
                self.assertEqual(len(coefficients), arity*p+1)
                self.assertTrue(all(isinstance(c, F) and c >= 0 for c in coefficients))

    def test_real_lambda_relation_provider_matches_tiny_field_exhaustion(self):
        h = matrix(lambda_value=2, q=4)
        for source in ((1,), (0, 1), (1, 3)):
            rows = tuple(tuple(row[i] for i in source) for row in h)
            joint = dict(Counter(direct_inputs(rows)))
            self.assertEqual(joint, relations.joint_counts(source, 'static', 0,
                                                         lambda_value=0xbc, q=4))
            expected = marked.marked_polynomials_from_joint(joint, len(source), q=4)
            actual = marked.mapped_source_polynomials(source, q=4)
            self.assertEqual(actual['polynomials'], expected)
            self.assertEqual(actual['raw_allowed_input_mass'], 4**len(source))

    def test_actual_gf256_mass_and_raw_zero_scaling(self):
        data = marked.mapped_source_polynomials((0, 1), orders=(2, 3, 4))
        joint = relations.joint_counts((0, 1), 'static', 0, lambda_value=0xbc)
        counts = Counter()
        for (_, mask), count in joint.items():
            counts[mask] += count
        for p, coefficients in data['polynomials'].items():
            normalized = F(256**(p*14), 256**16)*sum(
                (F(count**p, 255**((p-1)*mask.bit_count())) for mask, count in counts.items()), F())
            self.assertEqual(sum(coefficients), 256**(2*p)*normalized)
            self.assertEqual(coefficients[0], 256**(16*(p-1)))

    def test_local_symmetry_witness_and_no_mutable_cache_leak(self):
        first = marked.mapped_source_polynomials((1, 3), q=4)
        image = marked.orbits.orthogonal_maps()[19]
        moved = tuple(image[i] for i in (1, 3))
        second = marked.mapped_source_polynomials(moved, q=4)
        self.assertEqual(first['polynomials'], second['polynomials'])
        for row in (first, second):
            image = row['representative_to_source_permutation']
            self.assertEqual(tuple(sorted(image[i] for i in row['representative'])), row['mapped_source'])
            self.assertIn(image, marked.orbits.orthogonal_maps())
        first['polynomials'].clear()
        self.assertEqual(marked.mapped_source_polynomials((1, 3), q=4)['polynomials'], second['polynomials'])

    def test_all_label_orbits_partition_without_arithmetic_inventory(self):
        for arity in range(1, 6):
            seen = set()
            for orbit in marked.source_orbits(arity):
                rep = orbit['representative']
                for source, image in orbit['members']:
                    self.assertNotIn(source, seen)
                    seen.add(source)
                    self.assertEqual(tuple(sorted(image[i] for i in rep)), source)
            self.assertEqual(seen, set(combinations(range(16), arity)))

    def test_mock_inventory_retains_every_label_and_exact_coefficient_strings(self):
        calls = []
        def fake(source, *, orders):
            calls.append(source)
            a = len(source)
            return dict(polynomials={p: (F(256**(16*(p-1))),)+(F(1, 3),)*(a*p) for p in orders},
                        raw_allowed_input_mass=256**a, joint_entries=7, output_masks=5)
        with patch.object(marked, 'mapped_source_polynomials', side_effect=fake):
            data = marked.inventory(arities=(1, 2, 3, 4, 5), max_seconds=10)
        self.assertTrue(data['complete'])
        self.assertEqual(data['expected_labelled_sources'], 6884)
        self.assertEqual(len(data['rows']), 6884)
        self.assertEqual(len(calls), len(data['representatives']))
        self.assertEqual(len(set(calls)), len(calls))
        for row in data['rows']:
            representative = data['representatives'][row['orbit_index']]
            image = row['representative_to_source_permutation']
            self.assertEqual(sorted(image[i] for i in representative['representative']), row['mapped_source'])
            self.assertEqual(representative['coefficients']['3'][1], '1/3')
        self.assertEqual(data['completed_labelled_sources'], sum(comb(16, a) for a in range(1, 6)))
        self.assertFalse(data['whole_word_chronology_quotiented'])
        self.assertFalse(data['rs_normalization_included'])
        self.assertFalse(data['whole_code_certificate'])
        self.assertEqual(data['regional_provider_runs'], 0)

    def test_budget_and_source_change_are_not_complete_evidence(self):
        with patch.object(marked, 'perf_counter', side_effect=(0., 2., 2.)), \
                patch.object(marked, 'mapped_source_polynomials') as calculate:
            result = marked.inventory(arities=(1,), max_seconds=1)
        calculate.assert_not_called()
        self.assertFalse(result['complete'])
        self.assertEqual(result['rows'], [])
        with patch.object(marked, 'perf_counter', side_effect=(0., 2.)), \
                patch.object(marked, 'checked_dependencies', side_effect=({'a': 'one'}, {'a': 'two'})):
            with self.assertRaisesRegex(ValueError, 'Source changed'):
                marked.inventory(arities=(1,), max_seconds=1)

    def test_guards(self):
        good = {(0, 0): 1, (1, 1): 3}
        for orders in ((), (2, 2), (1,), (5,), (True,)):
            with self.assertRaises(ValueError):
                marked.marked_polynomials_from_joint(good, 1, width=1, q=4, orders=orders)
        for bad in ({(0, 0): 1, (1, 1): 2}, {(0, 0): 1, (1, 0): 3},
                    {(0, 1): 1, (1, 1): 3}, {(0, 0): 1, (1, 2): 3},
                    {(0, 0): 1, (True, 1): 3}, {(0, 0): 1, (1, 1): F(3)}):
            with self.assertRaises(ValueError):
                marked.marked_polynomials_from_joint(bad, 1, width=1, q=4)
        for source in ((), (0, 0), (16,), (True,), tuple(range(6))):
            with self.assertRaises(ValueError):
                marked.mapped_source_polynomials(source)
        for arities in ((), (1, 1), (6,), (True,)):
            with self.assertRaises(ValueError):
                marked.inventory(arities=arities, max_seconds=1)
        for q in (2, 16, True):
            with self.assertRaises(ValueError):
                marked.marked_polynomials_from_joint(good, 1, width=1, q=q)
        for budget in (0, 601, float('nan'), True):
            with self.assertRaises(ValueError):
                marked.inventory(arities=(1,), max_seconds=budget)


if __name__ == '__main__':
    unittest.main()

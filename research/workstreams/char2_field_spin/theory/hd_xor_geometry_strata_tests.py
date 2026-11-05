"""Small exact checks for fixed constituent geometry and ordered caps."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb
import unittest

import numpy as np
from flint import arb, ctx

from hd_xor_geometry_strata import (MultiplicityCapsArb, MultiplicityCapsFloat,
    category_resources, common_pair_caps_arb, common_pair_caps_float,
    geometry_stratum_counts, pair_resource_arb, pair_resource_float,
    scalar_route_receipt)
from hd_xor_joint_columns import resource_envelope_float


class GeometryStrataTests(unittest.TestCase):
    def test_integer_counts_equal_labelled_assignment_enumeration(self):
        for choices in range(1, 5):
            for groups in range(1, 7):
                exact = Counter(max(Counter(word).values()) for word in product(range(choices), repeat=groups))
                expected = {mu: exact[mu] for mu in range(1, groups+1)}
                self.assertEqual(geometry_stratum_counts(groups, choices), expected)

    def test_production_counts_partition_once_and_common_pair_has_28_choices(self):
        for groups in range(1, 33):
            counts = geometry_stratum_counts(groups)
            self.assertEqual(sum(counts.values()), 28**groups)
            self.assertEqual(counts[groups], 28)
            if groups <= 28:
                expected = 1
                for value in range(28-groups+1, 29):
                    expected *= value
                self.assertEqual(counts[1], expected)
            else:
                self.assertEqual(counts[1], 0)

    def test_multiplicity_constraints_hold_for_actual_cancelled_multigraphs(self):
        edges = tuple(combinations(range(4), 2))
        # Every small multigraph and cancellation subgraph, including repeated
        # edges. Singleton endpoints range over all eight labelled regions.
        for geometry in product(edges, repeat=3):
            mu = max(Counter(geometry).values())
            for single in range(3):
                double = 3-single
                for singleton_rows in product(range(8), repeat=single):
                    for cancelled_flags in product(range(2), repeat=double):
                        k = sum(cancelled_flags)
                        occupancy = [single+double-k-singleton_rows.count(row) for row in range(8)]
                        for edge, flag in zip(geometry[single:], cancelled_flags):
                            if flag:
                                for row in edge:
                                    occupancy[row] += 1
                        lower, upper, pair_limit, target = category_resources(single, double, k, mu)
                        self.assertGreaterEqual(min(occupancy), lower)
                        self.assertLessEqual(max(occupancy), upper)
                        self.assertEqual(sum(occupancy)-8*lower, target)
                        self.assertTrue(all(occupancy[i]+occupancy[j] <= pair_limit
                                            for i, j in combinations(range(8), 2)))
        # Degree concentration at one vertex gives the additional 7*mu bound.
        self.assertEqual(category_resources(0, 9, 9, 1)[1], 7)

    def test_ordered_resource_caps_dominate_all_words_and_agree_across_backends(self):
        common = np.array([[[1., 0.], [.125, .25]], [[.125, .5], [.0625, .25]],
                           [[.25, .125], [.03125, .5]], [[.125, .25], [.25, .125]]])
        # Distinct labelled tables catch left/right order reversals.
        labelled = np.stack([common, common[:, ::-1, ::-1], common*.5])
        with np.errstate(divide='ignore'):
            logs = np.log2(labelled)
        directed = [[tuple(arb(float(x)) for x in matrix.flatten()) for matrix in table] for table in labelled]
        for lower, upper, limit in ((0, 3, 3), (1, 3, 4), (0, 2, 4), (1, 1, 1)):
            floating = pair_resource_float(logs, -2., lower, upper, limit, 3)
            bounded = pair_resource_arb(directed, arb(.25), lower, upper, limit, 3, labelled_regions=True)
            for word in product(range(lower, upper+1), repeat=3):
                if any(a+b > limit for a, b in combinations(word, 2)):
                    continue
                vector = np.array([1., .25])
                for position in reversed(range(3)):
                    vector = labelled[position, word[position]] @ vector
                degree = sum(word)-3*lower
                self.assertTrue(all(arb(float(vector[row])) <= bounded[degree][row] for row in range(2)))
            np.testing.assert_allclose(np.exp2(floating),
                [[float(x) for x in vector] for vector in bounded], rtol=1e-14)

    def test_one_letter_has_no_pair_constraint_and_inactive_constraint_preserves_old_path(self):
        matrices = np.array([[[1., 0.], [0., .25]], [[.25, .5], [.125, .25]],
                             [[.125, .25], [.03125, .5]]])
        with np.errstate(divide='ignore'):
            logs = np.log2(matrices)
        actual = pair_resource_float(logs, -2., 0, 2, 0, 1)
        np.testing.assert_allclose(np.exp2(actual), [matrix @ [1., .25] for matrix in matrices])
        np.testing.assert_array_equal(pair_resource_float(logs, -2., 0, 2, 4, 8),
                                      resource_envelope_float(logs, -2., 0, 2, 8))

    def test_cancellation_mixture_precedes_potential_quotient(self):
        tables = np.ones((2, 2, 2))
        floating = MultiplicityCapsFloat(tables, -1., 1)
        directed = MultiplicityCapsArb([(arb(1),)*4]*2, arb(.5), 1)
        # d=1: k=0 reaches total0, k=1 reaches total2. The two branches
        # attain their maxima in different vector coordinates.
        def vectors(lower, upper, limit):
            return [(arb(254), arb(0))]*3 if lower else [(arb(0), arb(255))]*3
        directed._envelope = vectors
        floating._envelope = lambda *args: np.array([[-np.inf if x == 0 else float(x.log()/arb(2).log())
                                                      for x in vector] for vector in vectors(*args)])
        cap = directed(1)[0, 1]
        expected = arb(254*254)/255
        self.assertTrue(cap.overlaps(expected))
        self.assertLess(cap, arb(254*254)/255+2)
        self.assertAlmostEqual(floating(1)[0, 1], float(expected.log()/arb(2).log()), places=12)

    def test_invalid_tables_and_scope_are_rejected(self):
        with self.assertRaises(ValueError):
            geometry_stratum_counts(33)
        with self.assertRaises(ValueError):
            pair_resource_float(np.zeros((2, 2, 2)), 0., 0, 2, 4)
        with self.assertRaises(ValueError):
            MultiplicityCapsFloat(np.full((2, 2, 2), np.nan), 0., 1)
        with self.assertRaises(ValueError):
            MultiplicityCapsFloat(np.zeros((2, 2, 2)), 0., 1)(0)

    def test_directed_cache_rejects_cross_precision_reuse(self):
        with ctx.workprec(256):
            evaluator = MultiplicityCapsArb([(arb(1),)*4]*2, arb(.5), 1)
        with ctx.workprec(384):
            with self.assertRaisesRegex(ValueError, 'precision context'):
                evaluator(1)

    def test_common_pair_caps_equal_exact_labelled_source_and_coin_enumeration(self):
        groups = 2
        raw = [((1, 0), (0, Fraction(1, 4))),
               ((Fraction(1, 8), Fraction(1, 2)), (Fraction(1, 16), Fraction(1, 4))),
               ((Fraction(1, 4), Fraction(1, 8)), (Fraction(1, 32), Fraction(1, 2)))]
        tables = [raw if position % 3 else [tuple(reversed(matrix)) for matrix in raw]
                  for position in range(8)]
        expected = {(s, d): Fraction(0) for s in range(3) for d in range(3-s)}
        for left, right in combinations(range(8), 2):
            for sources in product(range(4), repeat=groups):
                single = sum(source in (1, 2) for source in sources)
                double = sources.count(3)
                mixed = [Fraction(0), Fraction(0)]
                for flags in product(range(2), repeat=double):
                    flags_iter = iter(flags)
                    masks = []
                    for source in sources:
                        masks.append(0 if source == 0 else 255 ^ (1 << left) if source == 1
                                     else 255 ^ (1 << right) if source == 2
                                     else (1 << left) | (1 << right) if next(flags_iter) else 255)
                    occupancies = [sum(mask >> position & 1 for mask in masks) for position in range(8)]
                    vector = [Fraction(1), Fraction(1, 4)]
                    for position in reversed(range(8)):
                        matrix = tables[position][occupancies[position]]
                        vector = [sum(matrix[row][col]*vector[col] for col in range(2)) for row in range(2)]
                    probability = Fraction(254**(double-sum(flags)), 255**double)
                    mixed = [mixed[row]+probability*vector[row] for row in range(2)]
                expected[single, double] = max(expected[single, double], mixed[0], 4*mixed[1])
        with np.errstate(divide='ignore'):
            logs = np.log2(np.array(tables, dtype=float))
        floating = common_pair_caps_float(logs, -2., groups)
        with ctx.workprec(256):
            directed_tables = [[tuple(arb(Fraction(x).numerator)/Fraction(x).denominator
                                      for row in matrix for x in row) for matrix in table] for table in tables]
            directed = common_pair_caps_arb(directed_tables, arb(1)/4, groups, labelled_regions=True)
            for category, exact in expected.items():
                value = arb(exact.numerator)/exact.denominator
                self.assertTrue(directed[category] >= value or directed[category].overlaps(value))
                self.assertAlmostEqual(2**floating[category], float(exact), places=14)
                self.assertAlmostEqual(float(directed[category]), float(exact), places=14)

    def test_small_directed_receipt_has_explicit_scalar_scope(self):
        receipt = scalar_route_receipt(1, precision=128)
        self.assertTrue(receipt['directed'])
        self.assertTrue(receipt['full_homogeneous_group_class'])
        self.assertFalse(receipt['paired_route_certificate'])
        self.assertFalse(receipt['whole_code_certificate'])
        self.assertEqual(receipt['source_pair_assignments'], 28)
        self.assertEqual(receipt['total_support_floor'], 34)
        self.assertEqual(receipt['projective_normalizations'], 1)
        Fraction(receipt['failure_log2_upper'])


if __name__ == '__main__':
    unittest.main()

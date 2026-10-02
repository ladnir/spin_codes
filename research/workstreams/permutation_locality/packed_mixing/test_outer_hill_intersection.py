import itertools
from math import comb
import unittest
from unittest.mock import patch

import canonical_counts as counts
import outer_hill_incidence as incidence
import outer_hill_intersection as intersection
from test_canonical_counts import binary_rank
from test_outer_hill_incidence import toy_code_and_dimensions


def synthetic_production_premises():
    # These fixtures exercise composition only; they are not authenticated
    # production premises and never produce a certificate.
    dimensions = [min(u, 2) if u <= 64 else min(u, 65) if u < 160
                  else min(128, 66+u-160) for u in range(257)]
    proof = dict(schema='canonical-octet-h5-exact-1', length=256, dimension=128,
                 block_width=8, minimum_distance=38, rank_cumulative_caps_at5=[0]*4)
    return dimensions, proof


class IntersectionTests(unittest.TestCase):
    def test_exact_weighted_formula_for_every_rank(self):
        actual = intersection.ekr_incidence_numerators(32, 128, 20, 66, 2)
        for rank, value in enumerate(actual, 1):
            self.assertEqual(value, (comb(32, 20)-comb(31, 11))*counts.rank_total(65, 4, rank)
                             + comb(31, 11)*counts.rank_total(66, 4, rank))
            self.assertLess(value, comb(32, 20)*counts.rank_total(66, 4, rank))
        # Summing ranks independently recovers the complete nonzero tuple count.
        self.assertEqual(sum(actual), comb(32, 20)*((1 << 260)-1)
                         + comb(31, 11)*((1 << 264)-(1 << 260)))

    def test_small_intersecting_families_exhaustively(self):
        for groups, size in ((4, 2), (5, 2)):
            subsets = [sum(1 << j for j in selected)
                       for selected in itertools.combinations(range(groups), size)]
            maximum = 0
            for bits in range(1 << len(subsets)):
                family = [value for i, value in enumerate(subsets) if bits >> i & 1]
                if all(a & b for a, b in itertools.combinations(family, 2)):
                    maximum = max(maximum, len(family))
            self.assertEqual(maximum, comb(groups-1, size-1))

    def test_toy_shortened_spaces_and_tuple_incidence(self):
        # Even-weight code on four of six coordinates, dimension3.
        code, dimensions = toy_code_and_dimensions([0b000011, 0b000101, 0b001001], 6)
        groups, size, high = 6, 4, dimensions[4]
        self.assertEqual((high, dimensions[2]), (3, 1))
        numerators = intersection.ekr_incidence_numerators(groups, 3, size, high,
                                                         dimensions[2], rows=2)
        actual = [0, 0]
        high_complements = []
        for selected in itertools.combinations(range(groups), size):
            mask = sum(1 << j for j in selected)
            subcode = [word for word in code if word & ~mask == 0]
            if len(subcode) == 1 << high:
                high_complements.append(((1 << groups)-1) ^ mask)
            for words in itertools.product(subcode, repeat=2):
                rank = binary_rank(words)
                if rank:
                    actual[rank-1] += 1
        self.assertTrue(all(a & b for a, b in itertools.combinations(high_complements, 2)))
        self.assertLessEqual(len(high_complements), comb(groups-1, groups-size-1))
        self.assertTrue(all(a <= b for a, b in zip(actual, numerators)))
        # Check the same incidence divisor against every actual support prefix.
        for h in range(size+1):
            exact = [0, 0]
            for words in itertools.product(code, repeat=2):
                rank = binary_rank(words)
                if rank and (words[0] | words[1]).bit_count() <= h:
                    exact[rank-1] += 1
            self.assertTrue(all(a <= b//comb(groups-h, size-h)
                                for a, b in zip(exact, numerators)))

    def test_strict_contradiction_and_ekr_domain_required(self):
        for args in ((6, 3, 4, 2, 1), (6, 3, 2, 3, 0),
                     (6, 3, 6, 3, 0), (True, 3, 4, 3, 1),
                     (6, 3, 4, 4, 1)):
            with self.assertRaises(ValueError):
                intersection.ekr_incidence_numerators(*args)

    def test_production_caps_only_tighten_and_preserve_h5_and_total(self):
        dimensions, proof = synthetic_production_premises()
        before = (dimensions[:], dict(proof))
        old = incidence.production_rank_caps(dimensions, proof)
        new = intersection.production_rank_caps(dimensions, proof)
        self.assertEqual((dimensions, proof), before)
        for previous, row in zip(old, new):
            self.assertEqual(row[:6], (0,)*6)
            self.assertEqual(row[-1], previous[-1])
            self.assertTrue(all(a <= b for a, b in zip(row, previous)))
            self.assertTrue(all(a <= b for a, b in zip(row, row[1:])))
        self.assertEqual(counts.total_caps(new)[-1], (1 << 512)-1)

    def test_production_rejects_missing_or_weaker_premises(self):
        dimensions, proof = synthetic_production_premises()
        bad64 = dimensions[:]
        bad64[64] = 3
        bad160 = dimensions[:]
        bad160[160] = 67
        bad_dimension = [min(d, 127) for d in dimensions]
        for bad in (dimensions[:-1], bad64, bad160, bad_dimension):
            with self.assertRaises(ValueError):
                intersection.production_rank_caps(bad, proof)
        with self.assertRaises(ValueError):
            intersection.production_rank_caps(dimensions, dict(proof, dimension=127))

    def test_authentication_is_fresh_and_new_schema_is_explicit(self):
        dimensions, proof = synthetic_production_premises()
        with patch.object(intersection, 'authenticated_bch_dimensions', return_value=dimensions) as auth, \
                patch.object(intersection, 'check_production', return_value=proof) as exact, \
                patch.object(intersection, 'transport_cdf', return_value=('synthetic',)) as transport:
            result, premises = intersection.authenticated_bch_cdf(17, False)
        auth.assert_called_once_with(17, False)
        exact.assert_called_once_with()
        transport.assert_called_once()
        self.assertEqual(result, ('synthetic',))
        self.assertEqual(premises['schema'], 'packed-canonical-full32-h5-incidence-ekr-expected-cdf-1')
        self.assertEqual(premises['shortened_dimensions'], dimensions)
        self.assertEqual(premises['count_refinements']['exact_h5'], proof)
        ekr = premises['count_refinements']['ekr']
        self.assertEqual(ekr['maximum_high_dimension_subsets'], comb(31, 11))
        self.assertEqual(ekr['rank_incidence_numerators'],
                         list(intersection.ekr_incidence_numerators(32, 128, 20, 66, 2)))
        self.assertFalse(ekr['assumes_lower_dimension_family_intersecting'])
        self.assertNotIn('verified', premises)


if __name__ == '__main__':
    unittest.main()

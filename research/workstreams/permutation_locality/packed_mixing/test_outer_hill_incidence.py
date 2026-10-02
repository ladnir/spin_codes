import itertools
import unittest
from unittest.mock import patch

import canonical_counts as counts
import outer_hill_incidence as hill
from local_models import full_block
from test_canonical_counts import binary_rank


def toy_code_and_dimensions(generators, length):
    code = [0]
    for row in generators:
        code += [word ^ row for word in code]
    dimensions = []
    for size in range(length+1):
        dimension = 0
        for positions in itertools.combinations(range(length), size):
            mask = sum(1 << i for i in positions)
            mass = sum(word & ~mask == 0 for word in code)
            dimension = max(dimension, mass.bit_length()-1)
        dimensions.append(dimension)
    return code, dimensions


class Incidence(unittest.TestCase):
    def test_exhaustive_toy_rank_cdfs(self):
        for generators in ([0b11000011, 0b00111100, 0b10101010],
                           [0b11110000, 0b11001100, 0b10101010, 0b11111111]):
            code, dimensions = toy_code_and_dimensions(generators, 8)
            for width in (1, 2, 4, 8):
                groups = 8//width
                actual = [[0]*(groups+1) for _ in range(2)]
                mask = (1 << width)-1
                for words in itertools.product(code, repeat=2):
                    rank = binary_rank(words)
                    if not rank:
                        continue
                    union = words[0] | words[1]
                    support = sum(bool(union & (mask << (width*b))) for b in range(groups))
                    for h in range(support, groups+1):
                        actual[rank-1][h] += 1
                old = counts.canonical_rank_caps(dimensions, width, 2)
                new = hill.incidence_rank_caps(dimensions, width, 2)
                for exact, cap, previous in zip(actual, new, old):
                    self.assertTrue(all(a <= b <= c for a, b, c in zip(exact, cap, previous)))
                    self.assertEqual(exact[-1], cap[-1])
                    self.assertTrue(all(a <= b for a, b in zip(cap, cap[1:])))
                local = full_block(width, 2)
                self.assertTrue(all(a <= b <= c for a, b, c in zip(
                    counts.transport_cdf(counts.total_caps(actual), local),
                    counts.transport_cdf(counts.total_caps(new), local),
                    counts.transport_cdf(counts.total_caps(old), local))))

    def test_strict_improvement_hamming_code(self):
        _, dimensions = toy_code_and_dimensions(
            [0b11110000, 0b11001100, 0b10101010, 0b11111111], 8)
        self.assertEqual(counts.canonical_rank_caps(dimensions, 1, 2)[0][4], 45)
        self.assertEqual(hill.incidence_rank_caps(dimensions, 1, 2)[0][4], 42)

    def test_optional_original_caps_and_no_mutation(self):
        _, dimensions = toy_code_and_dimensions([0b11110000, 0b11001100], 8)
        original = [list(row) for row in counts.canonical_rank_caps(dimensions, 1, 2)]
        before = (dimensions[:], [row[:] for row in original])
        new = hill.incidence_rank_caps(dimensions, 2, 2, original)
        self.assertTrue(all(cap[h] <= original[r][2*h]
                            for r, cap in enumerate(new) for h in range(5)))
        self.assertEqual((dimensions, original), before)

    def test_invalid_dimensions(self):
        for dimensions, width in (([0, 0, 1], 3), ([0, 1, 0], 1), ([0, True], 1)):
            with self.assertRaises(ValueError):
                hill.incidence_rank_caps(dimensions, width)

    def test_production_composition_keeps_dimensions(self):
        dimensions = [min(u, 128) for u in range(257)]
        proof = dict(schema='canonical-octet-h5-exact-1', length=256, dimension=128,
                     block_width=8, minimum_distance=38, rank_cumulative_caps_at5=[0]*4)
        before = dimensions[:]
        ranks = hill.production_rank_caps(dimensions, proof)
        self.assertEqual(dimensions, before)
        self.assertTrue(all(row[:6] == (0,)*6 for row in ranks))
        self.assertEqual(counts.total_caps(ranks)[-1], (1 << 512)-1)
        with self.assertRaises(ValueError):
            hill.production_rank_caps(dimensions, dict(proof, dimension=127))
        with self.assertRaises(ValueError):
            hill.production_rank_caps(dimensions, dict(proof, rank_cumulative_caps_at5=[-1]*4))

    def test_authenticated_entry_freshly_calls_both_checkers(self):
        dimensions = [min(u, 128) for u in range(257)]
        proof = dict(schema='canonical-octet-h5-exact-1', length=256, dimension=128,
                     block_width=8, minimum_distance=38, rank_cumulative_caps_at5=[0]*4)
        with patch.object(hill, 'authenticated_bch_dimensions', return_value=dimensions) as auth, \
                patch.object(hill, 'check_production', return_value=proof) as exact, \
                patch.object(hill, 'transport_cdf', return_value=('synthetic',)) as transport:
            result, premises = hill.authenticated_bch_cdf(17, False)
        auth.assert_called_once_with(17, False)
        exact.assert_called_once_with()
        transport.assert_called_once()
        self.assertEqual(result, ('synthetic',))
        self.assertEqual(premises['shortened_dimensions'], dimensions)
        self.assertEqual(premises['count_refinements']['exact_h5'], proof)


if __name__ == '__main__':
    unittest.main()

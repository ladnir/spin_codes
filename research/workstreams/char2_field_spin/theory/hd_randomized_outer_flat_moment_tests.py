"""Exact geometry and positive-DP checks for the new flat-state moment."""

from fractions import Fraction as F
from itertools import combinations
from math import comb
from random import Random
import unittest

from hd_randomized_outer_flat_moment import (
    flat_moments, profile_rank_lower, compare_selected,
)
from hd_blockwise_shared_transfer import rank as small_rank
from hd_blockwise_shared_maps import matrices, rank


def restricted_ranks(feedthrough, state, inputs, zeros, rank_function=rank):
    base = tuple(tuple(feedthrough[i][j] for j in inputs) for i in zeros)
    joint = tuple(row+tuple(state[i]) for row, i in zip(base, zeros))
    output_rank = rank_function(base)
    return output_rank, rank_function(joint)-output_rank


class FlatMomentTests(unittest.TestCase):
    def test_single_tile_cauchy_vandermonde_rank(self):
        feedthrough, state, _, _, _ = matrices(64)
        for tile in range(8):
            for a in range(9):
                for v in range(9):
                    inputs = tuple(8*tile+i for i in range(a))
                    zeros = tuple(8*tile+7-i for i in range(v))
                    output, residual = restricted_ranks(feedthrough, state, inputs, zeros)
                    self.assertEqual(output, min(a, v))
                    self.assertEqual(residual, min(4, max(v-a, 0)))

    def test_distinct_one_and_three_dimensional_tile_constraints(self):
        feedthrough, state, _, _, _ = matrices(64)
        generator = Random(0xB364)
        for first, second in combinations(range(8), 2):
            for dimension in (1, 3):
                a, v = (3, 4) if dimension == 1 else (1, 4)
                inputs = tuple(8*t+i for t in (first, second)
                               for i in generator.sample(range(8), a))
                zeros = tuple(8*t+i for t in (first, second)
                              for i in generator.sample(range(8), v))
                _, actual = restricted_ranks(feedthrough, state, inputs, zeros)
                self.assertGreaterEqual(actual, 2 if dimension == 1 else 4)

    def test_packet_rank_lower_on_bounded_actual_masks(self):
        feedthrough, state, _, _, _ = matrices(64)
        generator = Random(0xF1A7)
        for count in (0, 1, 8, 16, 24, 32, 40, 48, 56, 60, 63, 64):
            for zero_count in (8, 32, 61, 64):
                inputs = tuple(sorted(generator.sample(range(64), count)))
                zeros = tuple(sorted(generator.sample(range(64), zero_count)))
                a = tuple(sum(x//8 == t for x in inputs) for t in range(8))
                v = tuple(sum(x//8 == t for x in zeros) for t in range(8))
                base, actual = restricted_ranks(feedthrough, state, inputs, zeros)
                self.assertEqual(base, sum(min(x, y) for x, y in zip(a, v)))
                self.assertLessEqual(profile_rank_lower(a, v), actual)

    def test_exhaustive_toy_packet_and_exact_nonzero_state_removal(self):
        # GF8, two width2 tiles shifted by4, one constant state coordinate.
        feedthrough = ((5, 4, 0, 0), (4, 5, 0, 0),
                       (0, 0, 5, 4), (0, 0, 4, 5))
        state = ((1,),)*4
        mark = F(1, 3)
        exact = [F(0) for unused in range(5)]
        for input_mask in range(16):
            inputs = tuple(i for i in range(4) if input_mask >> i & 1)
            a = tuple(sum(i//2 == t for i in inputs) for t in range(2))
            for zero_mask in range(16):
                zeros = tuple(i for i in range(4) if zero_mask >> i & 1)
                v = tuple(sum(i//2 == t for i in zeros) for t in range(2))
                base, residual = restricted_ranks(feedthrough, state, inputs, zeros,
                                                  lambda matrix:small_rank(matrix, 8))
                lower = profile_rank_lower(a, v, width=2, state=1)
                self.assertEqual(lower, residual)
                probability = F(8**(1-residual)-1, 7*8**base)
                exact[len(inputs)] += mark**(4-len(zeros))*(1-mark)**len(zeros)*probability
        exact = tuple(value/comb(4, j) for j, value in enumerate(exact))
        from flint import arb, ctx
        previous = ctx.prec
        ctx.prec = 160
        try:
            result = flat_moments(mark, 8, width=2, tiles=2, state=1, directed=True)
            for value, reference in zip(result, exact):
                self.assertTrue((value-arb(reference.numerator)/reference.denominator).contains(0))
        finally:
            ctx.prec = previous

    def test_directed_moments_and_strict_middle_improvement(self):
        from flint import arb, ctx
        previous = ctx.prec
        ctx.prec = 160
        try:
            values = flat_moments(F(1, 256), directed=True)
            self.assertTrue(values[32].log()/arb(2).log() < -510)
            self.assertTrue(values[48].log()/arb(2).log() < -511)
            exact_last = ((1+arb(65535)/256)/65536)**64
            self.assertTrue((values[64]-exact_last).contains(0))
            algebraic = flat_moments(arb(2)**(-arb(1)/2), directed=True)
            self.assertTrue(all(0 < value <= 1 for value in algebraic))
        finally:
            ctx.prec = previous
        comparison = compare_selected(F(1, 256), (16, 24, 32, 40, 48))
        self.assertLessEqual(comparison['peak_states'], 1000)
        self.assertTrue(all(item['gain_bits'] > 60 for item in comparison['positions']))


if __name__ == '__main__':
    unittest.main()

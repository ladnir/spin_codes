"""Exact, bounded checks of the zero/one/two-input actual-mask rank kernels."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from hashlib import sha256
from math import comb
from random import Random
import unittest

import numpy as np

from cauchy_diagonal_packet_reference import LABELS, NODES
from cauchy_joint_return import matrix_rank, packet_forms, MUL, rank_return_profiles
from mask_state_rank import (joint_ranks_feasible, local_mask_transfer, _one_placement_terms,
                             _placement_terms, input_profiles, one_input_profiles,
                             placement_orbits, translate_mask)


class MaskRankTests(unittest.TestCase):
    def test_batched_rank_and_feasibility_against_scalar(self):
        rng = Random(901731)
        for width in range(1, 7):
            rows = max(1, width-1)
            matrices = [[[rng.randrange(256) for _ in range(width)] for _ in range(rows)]
                        for _ in range(24)]
            # Include coordinate-killing, dependent, and zero equation rows.
            matrices[0] = [[int(i == j) for j in range(width)] for i in range(rows)]
            matrices[1] = [[0]*width for _ in range(rows)]
            following = [[[rng.randrange(256) for _ in range(width)] for _ in range(4)]
                         for _ in matrices]
            expected, masks = [], []
            coordinate = [[int(i == j) for i in range(width)] for j in range(width)]
            for index, (equations, future) in enumerate(zip(matrices, following)):
                rank = matrix_rank(equations, width)
                mask = sum((matrix_rank(equations+[form], width) > rank) << j
                           for j, form in enumerate(future))
                target = mask if index % 2 else mask ^ 1
                alive = (all(matrix_rank(equations+[form], width) > rank for form in coordinate)
                         and target == mask)
                masks.append(target)
                expected.append((rank, alive))
            ranks, alive = joint_ranks_feasible(matrices, following, masks)
            self.assertEqual(list(zip(ranks.tolist(), alive.tolist())), expected)

    def test_zero_packet_preserves_actual_mask(self):
        y, tau = F(9, 10), F(7, 2)
        value = local_mask_transfer(0, y, tau)
        for h in range(16):
            for k in range(16):
                if h != k:
                    self.assertEqual(value[h][k], 0)
            expected = (F(1) if not h else sum(comb(16, z)*y**(16-z)*(1-y)**z/tau**z
                                               for z in range(h.bit_count())))
            self.assertEqual(value[h][h], expected)
        self.assertEqual(value[1][1], y**16)

    def test_one_packet_birth_and_joint_return(self):
        y, tau = F(9997, 10000), F(35)
        value = local_mask_transfer(1, y, tau, phase=0)
        self.assertEqual(value[0][15], y)
        self.assertEqual(sum(value[0]), y)
        self.assertEqual(value[15][0], y**16/tau**4)
        self.assertTrue(all(value[h][0] == 0 for h in range(1, 15)))
        # Outside the incoming support, a single full-support injection
        # cannot be cancelled by the refreshed incoming state.
        for h in range(1, 16):
            for k in range(16):
                if k & (15 ^ h) != 15 ^ h:
                    self.assertEqual(value[h][k], 0)
        # Independent phase/placement checks without building further tables.
        for phase in (1, 2, 3):
            for index in (0, 9):
                terms = _one_placement_terms(phase, 15, index)
                self.assertEqual({(rank, zeros): count for (k, rank, zeros), count in terms.items() if k == 0},
                                 {(4, 0): 1})

    def test_translation_orbits_and_original_one_input_profiles(self):
        for phase in range(4):
            for offset in range(4):
                for i in range(16):
                    for j in range(4):
                        self.assertEqual(LABELS[phase][i ^ offset] ^ NODES[j ^ offset],
                                         LABELS[phase][i] ^ NODES[j])
        for weight, number in ((1, 4), (2, 36)):
            orbits = placement_orbits(weight)
            self.assertEqual(len(orbits), number)
            images = [tuple(sorted(i ^ offset for i in inputs))
                      for inputs, offsets in orbits for offset in offsets]
            self.assertEqual(sorted(images), list(combinations(range(16), weight)))
        frozen = ('1d4dba3eb9b5c830a693bfeaf47de699fc7e9b50a2246b4ae69853fa1632c04e',
                  '4965a748888402e57508a014a6dd09141db5131604ec70e92ecab387b477436e',
                  'b87c1fc279d07470c368674f045e7bf7275a53a7673185d6b79025ee508d13c1',
                  'e08ef9453df4692fd8cd1bf16b514069e6a2c8ea2db451906e31da188174d8ca')
        for phase in range(4):
            self.assertEqual(sha256(repr(one_input_profiles(phase)[1:]).encode()).hexdigest(), frozen[phase])

    def test_two_input_fixed_placement_translation(self):
        for phase, incoming, inputs in ((0, 3, (0, 1)), (1, 7, (0, 6)), (2, 15, (0, 9)), (3, 5, (1, 14))):
            original = _placement_terms(phase, incoming, inputs)
            for offset in (1, 2, 3):
                actual = _placement_terms(phase, translate_mask(incoming, offset),
                                          tuple(sorted(i ^ offset for i in inputs)))
                expected = Counter({(translate_mask(k, offset), rank, zeros): count
                                    for (k, rank, zeros), count in original.items()})
                self.assertEqual(actual, expected)

    def test_two_input_all_phase_births_and_returns(self):
        # Return oracle uses a projective arrangement, not this joint matrix
        # enumeration. Birth oracle exhausts all 255 projective input ratios.
        y, tau = F(9, 10), F(35)
        for phase in range(4):
            profiles = input_profiles(2, phase)
            returns = Counter()
            for incoming in range(1, 16):
                for rank, zeros, count in profiles[incoming][0]:
                    returns[rank, zeros] += count
            oracle, = rank_return_profiles(2, phase=phase, average_inputs=True)
            self.assertEqual(returns, Counter({(rank, zeros): count*120 for rank, zeros, count in oracle}))
            feasible = Counter()
            for inputs in combinations(range(16), 2):
                injection, _ = packet_forms(phase, inputs)
                supports = set()
                for ratio in range(1, 256):
                    supports.add(sum((row[0] ^ MUL[row[1]][ratio] != 0) << j
                                     for j, row in enumerate(injection)))
                self.assertNotIn(0, supports)
                for mask in supports:
                    feasible[mask] += 1
            value = local_mask_transfer(2, y, tau, phase)
            for mask in range(16):
                expected = F(feasible[mask], 120)*y**2/tau**int(mask != 15)
                self.assertEqual(value[0][mask], expected)
            self.assertEqual(value[0][0], 0)

    def test_complete_profiles_use_one_mark_convention(self):
        y, tau = F(4, 5), F(19)
        for weight in (0, 1, 2):
            matrix = local_mask_transfer(weight, y, tau, 0)
            for incoming, row in enumerate(input_profiles(weight, 0)):
                for target, profile in enumerate(row):
                    expected = sum(F(count, comb(16, weight))*y**(16-z)*(1-y)**z/tau**rank
                                   for rank, z, count in profile)
                    self.assertEqual(matrix[incoming][target], expected)

    def test_joint_diagonal_rank_charge_small_field(self):
        def multiply(a, b):
            result = 0
            while b:
                if b & 1:
                    result ^= a
                b >>= 1
                a <<= 1
                if a & 4:
                    a ^= 7
            return result

        table = [[multiply(a, b) for b in range(4)] for a in range(4)]

        def rank(rows, width):
            work = [list(row) for row in rows]
            r = 0
            for column in range(width):
                pivot = next((j for j in range(r, len(work)) if work[j][column]), None)
                if pivot is None:
                    continue
                work[r], work[pivot] = work[pivot], work[r]
                inverse = next(a for a in range(1, 4) if table[a][work[r][column]] == 1)
                work[r] = [table[inverse][a] for a in work[r]]
                for j in range(r+1, len(work)):
                    coefficient = work[j][column]
                    work[j] = [a ^ table[coefficient][b] for a, b in zip(work[j], work[r])]
                r += 1
            return r

        rng = Random(602117)
        forms = [v for v in product(range(4), repeat=3) if any(v)]
        comparisons = 0
        for n in range(1, 5):
            for _ in range(8):
                coefficients = [[rng.randrange(4) for _ in range(n)] for _ in range(3)]
                source = [rng.choice(forms) for _ in range(n)]
                histogram = Counter()
                for diagonal in product(range(1, 4), repeat=n):
                    equations = []
                    for row in coefficients:
                        equation = [0]*3
                        for coefficient, d, form in zip(row, diagonal, source):
                            scalar = table[coefficient][d]
                            equation = [a ^ table[scalar][b] for a, b in zip(equation, form)]
                        equations.append(equation)
                    r = rank(equations, 3)
                    if all(rank(equations+[list(form)], 3) > r for form in source):
                        histogram[r] += 1
                charge = rank(coefficients, n)
                for tau in (1, 2, 3):
                    actual = sum(F(count, 3**n*tau**r) for r, count in histogram.items())
                    self.assertLessEqual(actual, F(1, tau**charge))
                    comparisons += 1
        self.assertEqual(comparisons, 96)

    def test_no_unproved_high_weight_fallback(self):
        for weight in (3, 4, 16):
            with self.assertRaises(NotImplementedError):
                local_mask_transfer(weight, F(9, 10), 35)


if __name__ == '__main__':
    unittest.main()

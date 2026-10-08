"""Independent finite checks for the new design gates; no native timings."""

from fractions import Fraction
import itertools
import math
import random
import unittest

import numpy as np

import rate_design_collision as collision
import rate_design_coset as coset
import rate_design_region_refresh as refresh


class CollisionTests(unittest.TestCase):
    def test_gap_count_against_positive_enumeration(self):
        for copies in range(2, 7):
            for weight in range(1, 5):
                counts = [0]*(weight*(copies-1)+1)
                for gaps in itertools.product(range(1, copies), repeat=weight):
                    counts[sum(gaps)] += math.prod(2*(copies-d) for d in gaps)
                running = 0
                for maximum, count in enumerate(counts):
                    running += count
                    self.assertEqual(collision.distance_sum_count(copies, weight, maximum), running)

    def test_pair_triple_and_two_pair_labels(self):
        for period, copies in ((2, 3), (3, 2), (2, 4)):
            rows = period*copies
            pair = Fraction(sum(a % period == b % period for a, b in itertools.permutations(range(rows), 2)), rows*(rows-1))
            triple = Fraction(sum(a % period == b % period == c % period for a, b, c in itertools.permutations(range(rows), 3)), rows*(rows-1)*(rows-2))
            double = Fraction(sum(a % period == b % period and c % period == d % period
                                  for a, b, c, d in itertools.permutations(range(rows), 4)), rows*(rows-1)*(rows-2)*(rows-3))
            self.assertEqual(pair, Fraction(copies-1, rows-1))
            self.assertEqual(triple, Fraction((copies-1)*(copies-2), (rows-1)*(rows-2)))
            self.assertEqual(double, pair*Fraction((copies-2)*(copies-3)+(period-1)*copies*(copies-1), (rows-2)*(rows-3)))

    def test_diagonal_ratio_constraints(self):
        # Ratios of two independent nonzero GF4 diagonals are independent
        # uniform symbols in its three-element multiplicative group.
        for weight in (2, 3):
            length = 2*weight
            first = tuple(range(weight))
            for overlap in range(weight+1):
                second = tuple(range(overlap))+tuple(range(weight, 2*weight-overlap))
                count = 0
                for ratios in itertools.product(range(3), repeat=length):
                    count += len({ratios[j] for j in first}) == 1 and len({ratios[j] for j in second}) == 1
                rank = 2*(weight-1) if overlap == 0 else 2*weight-overlap-1
                self.assertEqual(Fraction(count, 3**length), Fraction(1, 3**rank))

    def test_finite_no_go_and_scope(self):
        old = collision.collision_lower(20, 16)
        self.assertGreater(old['failure_lower'], Fraction(1, 1000))
        self.assertGreater(collision.collision_lower(20, 1024)['failure_lower'], Fraction(1, 1 << 40))
        new = collision.collision_lower(20, 16, 32, 24)
        self.assertLess(new['first_bonferroni'], Fraction(1, 1 << 40))
        self.assertLessEqual(new['failure_lower'], new['first_bonferroni'])
        with self.assertRaises(ValueError):
            collision.collision_lower(20, 15)

    def test_joint_match_tail_against_enumeration(self):
        marginal, joint = Fraction(1, 3), Fraction(1, 7)
        weights = [1-2*marginal+joint, marginal-joint, marginal-joint, joint]
        for length in range(1, 6):
            for minimum in range(length+1):
                direct = Fraction(0)
                for outcomes in itertools.product(range(4), repeat=length):
                    if sum(value & 1 for value in outcomes) >= minimum and sum(value >> 1 for value in outcomes) >= minimum:
                        direct += math.prod(weights[value] for value in outcomes)
                self.assertEqual(collision.joint_match_tail(length, minimum, marginal, joint), direct)

    def test_rank_pair_gate(self):
        for k_log in (16, 18, 20):
            self.assertGreater(collision.kernel_pair_lower(k_log)['failure_lower'], Fraction(1, 1 << 40))
        self.assertGreater(collision.kernel_pair_lower(20, 32)['failure_lower'], Fraction(1, 1 << 40))
        self.assertLess(collision.kernel_pair_lower(20, 64)['first_bonferroni'], Fraction(1, 1 << 40))


class CosetTests(unittest.TestCase):
    def test_subfield_partition_and_basis(self):
        points, betas = coset.geometry()
        self.assertEqual(len(points), 16)
        self.assertEqual(len(betas), 16)
        self.assertEqual({u ^ beta for beta in betas for u in points}, set(range(256)))
        for beta in betas:
            for state in ([0, 0, 0], [1, 2, 3], [255, 17, 128]):
                self.assertEqual(coset.basis(coset.basis(state, beta), beta), state)
                self.assertEqual(coset.basis(coset.basis(state, beta, True), beta, True), state)

    def test_global_local_equivalence(self):
        generator = random.Random(719)
        for beta in coset.geometry()[1]:
            for _ in range(4):
                x = [generator.randrange(256) for _ in range(16)]
                q = [generator.randrange(256) for _ in range(3)]
                global_y, global_q = coset.global_step(x, q, beta)
                local_y, local_s = coset.local_step(x, coset.basis(q, beta))
                self.assertEqual(global_y, local_y)
                self.assertEqual(coset.basis(global_q, beta), local_s)

    def test_adjoint_with_region_changes(self):
        generator = random.Random(613)
        betas = [coset.geometry()[1][j] for j in (0, 1, 9, 15)]
        x = [[[generator.randrange(256) for _ in range(16)] for _ in range(3)] for _ in betas]
        y = [[[generator.randrange(256) for _ in range(16)] for _ in range(3)] for _ in betas]
        ex, aty = coset.encode(x, betas), coset.adjoint(y, betas)
        left = coset.xor_sum(coset.mul(a, b) for rr, ss in zip(ex, y) for r, s in zip(rr, ss) for a, b in zip(r, s))
        right = coset.xor_sum(coset.mul(a, b) for rr, ss in zip(x, aty) for r, s in zip(rr, ss) for a, b in zip(r, s))
        self.assertEqual(left, right)

    def test_partial_first_last_regions(self):
        points, betas = coset.geometry()
        steps = 4
        regions = [[[0]*16 for _ in range(steps)] for _ in betas]
        # Four distinct-region Vandermonde columns have this nonzero kernel
        # vector. Place first/last injections at the worst packet endpoints.
        alphas = betas[:4]
        values = []
        for i, alpha in enumerate(alphas):
            denominator = 1
            for j, other in enumerate(alphas):
                if i != j:
                    denominator = coset.mul(denominator, alpha ^ other)
            values.append(coset.power(denominator, 254))
        for region, packet, value in zip(range(4), (steps-1, 1, 2, 0), values):
            regions[region][packet][0] = value
        regions[15][steps-1][0] = 1
        output = coset.encode(regions, betas)
        weight = sum(value != 0 for region in output for packet in region for value in packet)
        self.assertGreaterEqual(weight, coset.one_row_lower(steps))
        self.assertTrue(all(value == 0 for region in output[4:15] for packet in region for value in packet))
        for steps in (171, 683, 2731):
            self.assertGreater(coset.one_row_lower(steps)*10, 16*16*steps)


class RegionalRefreshTests(unittest.TestCase):
    def test_nonadjacent_return_count(self):
        for weight in range(1, 13):
            for returns in range(weight//2+1):
                valid = sum(all(b-a > 1 for a, b in zip(selected, selected[1:]))
                            for selected in itertools.combinations(range(1, weight), returns))
                self.assertEqual(valid, math.comb(weight-returns, returns))

    def test_exact_one_row_bound(self):
        for k_log, period in itertools.product((16, 18, 20), (16, 64)):
            result = refresh.return_count_upper(k_log, period=period)
            self.assertEqual(result['permitted_intermediate_regions'], 3)
            self.assertEqual(result['rows'] % period, 0)
            self.assertLess(result['failure_upper'], Fraction(1, 1 << 40))
            self.assertEqual(result['failure_upper'], sum(entry['contribution'] for entry in result['entries']))
            self.assertEqual(result['entries'][1]['required_returns'], 3)

    def test_window_firstmoment_obstruction(self):
        for period in (16, 64):
            result = refresh.window_firstmoment_obstruction(period=period)
            self.assertEqual(result['active'], 12)
            self.assertEqual(result['kernel_dimension'], 164)
            self.assertLessEqual(result['total_window'], result['bad_weight_threshold'])
            self.assertGreater(result['expected_bad_words_lower'], 1 << 27)
            self.assertFalse(result['is_failure_probability'])
            self.assertLess(result['route_probability'], Fraction(1, 1 << 1200))

    def test_regional_matrices_against_exact_small_field(self):
        from large_state_screen import mul
        q, steps, state, z = 16, 2, 2, Fraction(3, 4)
        exact_zero, exact_return, exact_continue = Fraction(0), Fraction(0), Fraction(0)
        for a, b in itertools.product(range(q), repeat=2):
            if not (a or b):
                continue
            before = [a ^ mul(j, b) for j in range(16)]
            weight_before = sum(v != 0 for v in before)
            exact_zero += z**(steps*weight_before)
            for lane in range(16):
                updated = [a ^ 1, b ^ lane]
                after = [updated[0] ^ mul(j, updated[1]) for j in range(16)]
                weight_after = sum(v != 0 for v in after)
                weight_active = sum((value ^ int(j == lane)) != 0 for j, value in enumerate(before))
                for packet in range(steps):
                    marked = z**(packet*weight_before+weight_active+(steps-packet-1)*weight_after)
                    if any(updated):
                        exact_continue += marked
                    else:
                        exact_return += marked
        exact_zero /= q**state-1
        exact_return /= (q**state-1)*16*steps
        exact_continue /= (q**state-1)*16*steps
        r0, r1 = refresh.region_matrices(float(z), steps, state, q)
        self.assertAlmostEqual(r0[1, 1], float(exact_zero), places=14)
        self.assertAlmostEqual(r1[1, 0], float(exact_return), places=14)
        self.assertGreaterEqual(r1[1, 1]*(1+1e-14), float(exact_continue))


if __name__ == '__main__':
    unittest.main()

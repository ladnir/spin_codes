from fractions import Fraction as Rat
from itertools import product
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_post_refresh_packet as cap


def multiply(a, b):
    answer = 0
    while b:
        if b & 1:
            answer ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return answer


def branch(row, u, c, target):
    """Classify the actual positive branches, not overlapping cone supports."""
    if row == cap.Z:
        return cap.Z if not c else cap.L
    if row == cap.L:
        return cap.D if not c else cap.R if not target[1] else cap.B
    if row == cap.R:
        return cap.Z if u == c else cap.L
    if row == cap.B:
        return cap.D if u == c else cap.R if not target[1] else cap.B
    if row == cap.D:
        return cap.P if u == c else cap.R if not target[1] else cap.B
    if row == cap.P:
        return cap.P if u == c else cap.G if not target[1] else cap.F
    if row == cap.F:
        return cap.D if u == c else cap.R if not target[1] else cap.B
    if row == cap.G:
        return cap.Z if u == c else cap.L
    raise ValueError('Unknown input cone')


def propagate(row, measure, support, z, fixed_input=None):
    assignments = ([fixed_input] if fixed_input is not None else
        [tuple(values[support.index(i)] if i in support else 0 for i in range(2))
         for values in product(range(1, 4), repeat=len(support))])
    result = {state: {} for state in range(8)}
    for (u, v), mass in measure.items():
        for x in assignments:
            mark = z**sum(value != 0 for value in (x[0] ^ u, x[1] ^ u))
            c = x[0] ^ x[1]
            for refresh in range(1, 4):
                target = v, multiply(refresh, u ^ c) ^ v
                output = result[branch(row, u, c, target)]
                output[target] = output.get(target, Rat(0))+mass*mark/(3*len(assignments))
    return result


class PostRefreshPacketTests(TestCase):
    def check_measure(self, row, measure, support, z, input_bound=None):
        matrix = cap.packet_caps((len(support),), z, q=4, length=2, rank=1)
        if input_bound is None:
            input_bound = sum(measure.values())
        actual = propagate(row, measure, support, z)
        for target, distribution in actual.items():
            coefficient = input_bound*matrix[8*row+target]
            if target in (cap.P, cap.G):
                self.assertLessEqual(sum(distribution.values()), coefficient)
                self.assertTrue(all(u and (v == u if target == cap.P else v == 0)
                                    for u, v in distribution))
            elif target == cap.F:
                # Minimal mass of a dominating conditional-flat family:
                # choose first-coordinate mass 3*max_v actual(u,v).
                self.assertTrue(all(u and v for u, v in distribution))
                required = sum(3*max((distribution.get((u, v), Rat(0))
                    for v in range(1, 4)), default=Rat(0)) for u in range(1, 4))
                self.assertLessEqual(required, coefficient)
            else:
                for (u, v), mass in distribution.items():
                    if target == cap.Z:
                        self.assertEqual((u, v), (0, 0))
                        denominator = 1
                    elif target == cap.L:
                        self.assertTrue(not u and v)
                        denominator = 3
                    elif target == cap.R:
                        self.assertTrue(u and not v)
                        denominator = 3
                    elif target == cap.D:
                        self.assertTrue(u and u == v)
                        denominator = 3
                    else:
                        self.assertEqual(target, cap.B)
                        self.assertTrue(u and v)
                        denominator = 9
                    self.assertLessEqual(mass, coefficient/denominator,
                        (row, target, support, z, (u, v), mass, coefficient))

    def test_all_eight_q4_cones_and_every_input_support(self):
        measures = [(cap.Z, {(0, 0): Rat(1)}),
            (cap.L, {(0, v): Rat(1, 3) for v in range(1, 4)}),
            (cap.R, {(u, 0): Rat(1, 3) for u in range(1, 4)}),
            (cap.B, {(u, v): Rat(1, 9) for u, v in product(range(1, 4), repeat=2)}),
            (cap.D, {(u, u): Rat(1, 3) for u in range(1, 4)})]
        for u in range(1, 4):
            measures.extend(((cap.P, {(u, u): Rat(1)}),
                (cap.G, {(u, 0): Rat(1)}),
                (cap.F, {(u, v): Rat(1, 3) for v in range(1, 4)})))
        for mask in range(4):
            support = tuple(i for i in range(2) if mask >> i & 1)
            for z in (Rat(1, 3), Rat(3, 4), Rat(1)):
                for row, measure in measures:
                    self.check_measure(row, measure, support, z)

    def test_weighted_arbitrary_point_and_conditional_flat_histories(self):
        first = {1: Rat(2, 7), 2: Rat(3, 5), 3: Rat(4, 3)}
        mixtures = [(cap.P, {(u, u): mass for u, mass in first.items()}, sum(first.values())),
            (cap.G, {(u, 0): mass for u, mass in first.items()}, sum(first.values())),
            (cap.F, {(u, v): mass/3 for u, mass in first.items() for v in range(1, 4)},
             sum(first.values()))]
        # A weighted actual history may have different second-coordinate
        # conditionals for each u, provided each is dominated by mass(u)*mu.
        conditional = {(u, v): mass/3*Rat((u+2*v) % 3, 2)
            for u, mass in first.items() for v in range(1, 4)}
        mixtures.append((cap.F, conditional, sum(first.values())))
        for row, measure, bound in mixtures:
            for support in ((), (0,), (1,), (0, 1)):
                self.check_measure(row, measure, support, Rat(2, 3), bound)

    def test_diagonal_spike_persists_but_first_point_pays_refresh_hit(self):
        z = Rat(3, 4)
        diagonal = {(u, u): Rat(1, 3) for u in range(1, 4)}
        actual = propagate(cap.D, diagonal, (0,), z, fixed_input=(1, 0))
        self.assertEqual(actual[cap.P], {(1, 1): z/3})
        self.assertGreater(actual[cap.P][1, 1], (z/3)/3)
        same = propagate(cap.P, {(1, 1): Rat(1)}, (0,), z, fixed_input=(1, 0))
        self.assertEqual(same[cap.P], {(1, 1): z})
        different = propagate(cap.P, {(1, 1): Rat(1)}, (0,), z, fixed_input=(2, 0))
        self.assertEqual(different[cap.G], {(1, 0): z**2/3})
        # F is a conditional-flat comparison, not an unrestricted pair mass.
        self.assertEqual(sum(different[cap.F].values()), 2*z**2/3)
        self.assertEqual(3*max(different[cap.F].values()), z**2)

    def test_exact_three_packet_return_needs_two_refresh_hits(self):
        for feedback in product(range(4), repeat=3):
            count = 0
            for keys in product(range(1, 4), repeat=3):
                u = v = 0
                for c, key in zip(feedback, keys):
                    u, v = v, multiply(key, u ^ c) ^ v
                a, b = multiply(keys[0], feedback[0]), multiply(keys[1], feedback[1])
                self.assertEqual((u, v) == (0, 0), a == b == feedback[2])
                count += (u, v) == (0, 0)
            expected = Rat(1) if not any(feedback) else Rat(1, 9) if all(feedback) else Rat(0)
            self.assertEqual(Rat(count, 27), expected)

    def test_generic_components_and_selected_D_decomposition(self):
        z = Rat(3, 4)
        details = cap.packet_details((1, 0, 0, 0), z)
        matrix = cap.packet_caps((1, 0, 0, 0), z)
        size, beta = details['dimension'], details['beta']
        self.assertEqual(details['sup_point_return'], z**12/255)
        self.assertEqual(matrix[8*cap.P+cap.G], beta*details['sup_point_live']/size)
        self.assertEqual(matrix[8*cap.D+cap.P], details['old_return_cap'])
        self.assertEqual(matrix[8*cap.D+cap.F], 0)
        self.assertEqual(matrix[8*cap.D+cap.G], 0)
        self.assertEqual(matrix[8*cap.B+cap.R], matrix[8*cap.B+cap.B]/size)
        self.assertEqual(cap.metadata()['terminal_vector'], [1]*8)

    def test_mocked_production_reuses_old_row_only_at_proved_entries(self):
        # Wiring only: no production table or regional bank is constructed.
        occupancy = np.array([[1, 0, 0, 0]])
        old = np.array([[[1., 0.], [1e-100, 1e-90]]])
        with patch('hd_xor_multigroup.profile_arrays', return_value=(occupancy, np.ones(1))), \
             patch('hd_xor_multigroup.nonzero_matrices', return_value=(np.zeros(1), np.ones(1), old)):
            matrix = cap.raw_packet_table(.9)[1, 0, 0, 0]
        for i, j in ((cap.R, cap.Z), (cap.B, cap.D), (cap.D, cap.P)):
            self.assertEqual(matrix[i, j], 1e-100)
        for i, j in ((cap.R, cap.L), (cap.B, cap.B)):
            self.assertEqual(matrix[i, j], 1e-90)
        self.assertEqual(matrix[cap.B, cap.R], 1e-90/(256**16-1))
        self.assertGreater(matrix[cap.D, cap.B], 1e-90)
        self.assertEqual(matrix[cap.D, cap.F], 0.)
        self.assertTrue(cap.metadata()['inverse_refresh_empty_cap_used_for_B'])

    def test_guards_and_scope(self):
        for occupancy in ((-1,), (3,), ()):
            with self.assertRaises(ValueError):
                cap.packet_caps(occupancy, Rat(1, 2), q=4, length=2, rank=1)
        with self.assertRaises(ValueError):
            cap.packet_caps((0,), 0, q=4, length=2, rank=1)
        with self.assertRaises(ValueError):
            cap.raw_packet_table(True)
        self.assertEqual(cap.STATE_ORDER, ('Z', 'L', 'R', 'B', 'D', 'P', 'F', 'G'))
        self.assertFalse(cap.metadata()['whole_code_certificate'])
        self.assertFalse(cap.metadata()['state_scaling_applied'])


if __name__ == '__main__':
    main()

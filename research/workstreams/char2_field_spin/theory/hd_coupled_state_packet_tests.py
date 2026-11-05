from fractions import Fraction as F
from itertools import product
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_coupled_state_packet as cap


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


def propagate(measure, support, z, fixed_input=None):
    assignments = ([fixed_input] if fixed_input is not None else
                   [tuple(values[support.index(i)] if i in support else 0 for i in range(2))
                    for values in product(range(1, 4), repeat=len(support))])
    answer = {}
    for (u, v), mass in measure.items():
        for x in assignments:
            y = x[0] ^ u, x[1] ^ u
            marked = z**sum(value != 0 for value in y)
            c = x[0] ^ x[1]
            for refresh in range(1, 4):
                target = v, multiply(refresh, u ^ v ^ c)
                answer[target] = answer.get(target, F(0))+mass*marked/(3*len(assignments))
    return answer


class CoupledPacketTests(TestCase):
    def check_measure(self, row, measure, support, z, fixed_input=None):
        matrix = cap.packet_caps((len(support),), z, q=4, length=2, rank=1)
        total_mass = sum(measure.values())
        result = propagate(measure, support, z, fixed_input)
        point_mass = sum(value for (u, v), value in result.items() if u and not v)
        if row == cap.L:
            self.assertLessEqual(point_mass, total_mass*matrix[5*row+cap.P])
        for (u, v), value in result.items():
            if not u and not v:
                bound = matrix[5*row+cap.Z]
            elif not u:
                bound = matrix[5*row+cap.L]/3
            elif not v:
                if row == cap.L:
                    continue  # P is constrained by total mass, not flat density.
                bound = matrix[5*row+cap.R]/3
            else:
                bound = matrix[5*row+cap.B]/9
            self.assertLessEqual(value, total_mass*bound,
                                 (row, support, z, (u, v), value, bound))

    def test_all_q4_scalar_states_inputs_and_refreshes(self):
        measures = [(cap.Z, {(0, 0): F(1)}),
                    (cap.L, {(0, v): F(1, 3) for v in range(1, 4)}),
                    (cap.R, {(u, 0): F(1, 3) for u in range(1, 4)}),
                    (cap.B, {(u, v): F(1, 9) for u, v in product(range(1, 4), repeat=2)})]
        measures += [(cap.P, {(u, 0): F(1)}) for u in range(1, 4)]
        for mask in range(4):
            support = tuple(i for i in range(2) if mask >> i & 1)
            for z in (F(1, 3), F(3, 4), F(1)):
                for row, measure in measures:
                    self.check_measure(row, measure, support, z)

    def test_weighted_point_history_is_not_a_new_union(self):
        measure = {(1, 0): F(2, 7), (2, 0): F(3, 5), (3, 0): F(4, 3)}
        for support in ((), (0,), (0, 1)):
            self.check_measure(cap.P, measure, support, F(2, 3))

    def test_fixed_input_cancellation_leaves_an_atom(self):
        z = F(3, 4)
        measure = {(0, v): F(1, 3) for v in range(1, 4)}
        result = propagate(measure, (0,), z, fixed_input=(1, 0))
        self.assertEqual(result[1, 0], z/3)
        self.assertEqual(sum(value for (u, v), value in result.items() if u and not v), z/3)
        # A P mass z/3 is valid; treating that same number as flat density
        # would assert mass z/9 at this point and fail.
        self.assertGreater(result[1, 0], (z/3)/3)
        self.check_measure(cap.L, measure, (0,), z, fixed_input=(1, 0))

    def test_rank_and_return_floors(self):
        one = cap.packet_details((1, 0, 0, 0), F(3, 4))
        self.assertEqual(one['rank'], 1)
        self.assertEqual(one['return_weight_floor'], 12)
        self.assertEqual(one['sup_point_return'], F(3, 4)**12/255)
        self.assertEqual(cap.support_floors((4, 4, 4, 4), 16, 4), (21, 36))
        empty = cap.packet_caps((0,), F(3, 4), q=4, length=2, rank=1)
        self.assertEqual(empty[5*cap.P+cap.Z], 0)
        self.assertEqual(empty[5*cap.P+cap.L], F(3, 4)**2)

    def test_unsharpened_uniform_emission_is_positive_exact_control(self):
        z = F(2, 3)
        for weight in range(3):
            support = tuple(range(weight))
            assignments = [tuple(values[i] if i < weight else 0 for i in range(2))
                           for values in product(range(1, 4), repeat=weight)]
            direct = sum(z**sum(value != 0 for value in (x[0] ^ u, x[1] ^ u))
                         for u in range(1, 4) for x in assignments)/(3*len(assignments))
            value = cap.packet_details((weight,), z, q=4, length=2, rank=1)
            self.assertGreaterEqual(value['uniform_live_emission'], direct)
        self.assertEqual(cap.packet_details((0,), z, q=4, length=2, rank=1)['uniform_live_emission'], z**2)

    def test_production_old_row_reuse_does_not_change_B(self):
        # Mocked table only; this does NOT construct the production bank.
        occupancy = np.array([[0, 0, 0, 0]])
        old = np.array([[[1., 0.], [0., 1e-30]]])
        with patch('hd_xor_multigroup.profile_arrays', return_value=(occupancy, np.ones(1))), \
             patch('hd_xor_multigroup.nonzero_matrices', return_value=(np.zeros(1), np.ones(1), old)):
            matrix = cap.raw_packet_table(.9)[0, 0, 0, 0]
        self.assertEqual(matrix[cap.R, cap.L], 1e-30)
        self.assertGreater(matrix[cap.B, cap.B], 1e-30)
        self.assertTrue(cap.metadata()['point_cone_is_arbitrary_mass'])
        self.assertFalse(cap.metadata()['inverse_refresh_empty_cap_used_for_B'])

    def test_guards(self):
        for occupancy in ((-1,), (3,), ()):
            with self.assertRaises(ValueError):
                cap.packet_caps(occupancy, F(1, 2), q=4, length=2, rank=1)
        with self.assertRaises(ValueError):
            cap.packet_caps((0,), 0, q=4, length=2, rank=1)
        with self.assertRaises(ValueError):
            cap.raw_packet_table(True)


if __name__ == '__main__':
    main()

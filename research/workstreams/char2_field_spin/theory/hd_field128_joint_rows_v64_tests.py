"""Unexecuted tiny tests; root runs remotely."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest
from unittest.mock import patch

import hd_field128_joint_rows_v64 as new
from hd_field128_wide_state_v60_tests import mul, mv, rank


class Tests(unittest.TestCase):
    def test_independent_joint_zero_rows(self):
        # GF16, width3/state2: CA!=0, all principal blocks invertible.
        a = ((1, 0), (0, 1), (1, 1))
        c = tuple(zip(*a))
        p = tuple(zip(*(mv(a, col) for col in zip(*c))))
        l = tuple(tuple(int(i == k) ^ mul(2, p[i][k]) for k in range(3))
                  for i in range(3))
        for j in range(4):
            table = new._row_table(j, 16, 3, 2)
            for support in combinations(range(3), j):
                self.assertEqual(rank([[l[i][k] for k in support] for i in support]), j)
                outside = tuple(i for i in range(3) if i not in support)
                joint = [[0]*(4-j) for unused in range(j+1)]
                for values in product(range(16), repeat=j):
                    x = [0]*3
                    for i, value in zip(support, values):
                        x[i] = value
                    if not any(mv(c, x)):
                        continue
                    y = mv(l, x)
                    az = sum(y[i] == 0 for i in support)
                    bz = sum(y[i] == 0 for i in outside)
                    joint[az][bz] += 1
                self.assertEqual(table.masses,
                    tuple(F(sum(row), 16**j) for row in joint))
                for row in table.rows:
                    counts = joint[row.active_zeros]
                    for e, bound in enumerate(row.bounds):
                        actual = F(sum(comb(b, e)*count for b, count in enumerate(counts)
                                       if b >= e), 16**j)
                        self.assertLessEqual(actual, bound)
                for z in (F(1, 3), F(1, 2), F(3, 4)):
                    direct = sum((F(count, 16**j)*z**(3-az-bz)
                        for az, row in enumerate(joint)
                        for bz, count in enumerate(row)), F(0))
                    refined, old, row_cap = new._components(j, z, 16, 3, 2)
                    self.assertLessEqual(direct, refined[2])
                    self.assertLessEqual(refined[2], old)
                    self.assertLessEqual(direct, row_cap)

    def test_generic_inversion_and_v57_row_evaluation(self):
        for n, state in ((3, 2), (5, 3), (6, 2)):
            for j in range(n+1):
                table = new._row_table(j, 16, n, state)
                self.assertEqual(sum(table.masses), 1-F(1, 16**min(j, state)))
                self.assertEqual(table.masses[-1], 0)
                if j <= state:
                    self.assertEqual(table.masses[:-1], tuple(
                        F(comb(j, a)*15**(j-a), 16**j) for a in range(j)))
                for z in (F(1, 4), F(3, 4)):
                    separate = sum((z**(j-row.active_zeros)
                        *new.old.tail.evaluate_envelope(row.bounds, z)['refined']
                        for row in table.rows), F(0))
                    refined, old, got = new._components(j, z, 16, n, state)
                    self.assertEqual(got, separate)
                    self.assertEqual(refined[2], min(old, separate))
                    base = new.old._moments(j, z, 16, n, state, True)
                    self.assertEqual((refined[0], refined[1], refined[3]),
                                     (base[0], base[1], base[3]))

    def test_public_normalization_and_guards(self):
        local, z = (F(1, 3), F(1, 5), F(1, 7), F(1, 11)), F(1, 2)
        for state in new.STATES:
            d = new.FIELD_ORDER**state-1
            rho = F(d, (new.FIELD_ORDER-1)**state)
            with patch.object(new, 'exact_moments', return_value=local):
                self.assertEqual(new.exact_packet(1, z, state),
                    (local[0], rho*local[1], local[2]/d, rho*local[3]))
                self.assertEqual(new.exact_packet(0, z, state)[3], z**(65-state))
            self.assertEqual(new.principal_bad_gamma_upper(state),
                             new.old.principal_bad_gamma_upper(state))
        for state in (True, 31, 33, 64):
            with self.assertRaises(ValueError):
                new.exact_packet(0, z, state)
        for j, z in ((True, F(1, 2)), (65, F(1, 2)), (1, .5), (1, True), (1, 0), (1, 1)):
            with self.assertRaises(ValueError):
                new.exact_packet(j, z, 40)

    def test_cache_source_and_batch_boundaries(self):
        self.assertIs(new._row_table(2, 16, 3, 2), new._row_table(2, 16, 3, 2))
        pins = new.source_pins()
        new.check_source_pins(pins)
        with patch.dict(new.FROZEN, {'../HD_FIELD128_JOINT_ROWS.md': '0'*64}):
            with self.assertRaises(RuntimeError):
                new.verify_frozen()
        local = (F(1), F(0), F(0), F(1, 3))
        with patch.object(new, 'verify_frozen') as check, \
             patch.object(new, '_components', return_value=(local, F(0), F(0))) as core:
            rows = new.local_logs(F(1, 2), 40)
            self.assertEqual(rows.shape, (65, 2, 2))
            self.assertEqual(check.call_count, 2)
            self.assertEqual([call.args[0] for call in core.call_args_list], list(range(65)))
            self.assertEqual(rows[0, 0, 1], float('-inf'))
            self.assertFalse(rows.flags.writeable)


if __name__ == '__main__':
    unittest.main()

import unittest
from fractions import Fraction as Q
from itertools import combinations

from constraint_counts import Constraints, normalize, verify_dual, dual_residual
from dual_moments import exact_tuple_cdfs
from joint_support import span
from bch_joint_support import rank_total


def dimensions(words, n):
    result = []
    for t in range(n+1):
        best = 0
        for positions in combinations(range(n), t):
            mask = sum(1 << j for j in positions)
            size = sum(w & ~mask == 0 for w in words)
            best = max(best, size.bit_length()-1)
        result.append(best)
    return result


class ConstraintCountTests(unittest.TestCase):
    def test_exact_box_residual(self):
        constraints = [normalize({0: 1, 1: 2}, 1)]
        # A deliberately imperfect rounded dual remains a valid bound.
        row, rhs = constraints[0]
        for multiplier in (Q(0), Q(1), Q(12)-Q(1, 1024), Q(20)):
            upper, correction = verify_dual(constraints, [], {0: Q(3), 1: Q(1)},
                                           [multiplier], [], 2)
            self.assertGreaterEqual(upper, 3)
            self.assertGreaterEqual(correction, 0)
            for x, y in ((Q(0), Q(0)), (Q(1), Q(0)), (Q(0), Q(1, 2))):
                self.assertLessEqual(3*x+y, upper)
        with self.assertRaises(ValueError):
            verify_dual(constraints, [], {0: Q(1)}, [Q(-1)], [], 2)

    def test_complement_and_containment_rows_on_exact_codes(self):
        for basis, n in (([3, 5], 4), ([15, 51, 85], 7), ([0x97, 0x4b, 0x2d, 0x1e], 8)):
            code = span(basis)
            dual = [w for w in range(1 << n) if all((w & b).bit_count() % 2 == 0 for b in basis)]
            spectra = [[sum(w.bit_count() == u for w in words) for u in range(n+1)] for words in (code, dual)]
            exact = [exact_tuple_cdfs(words, n, 4) for words in (code, dual)]
            inflated = [[[2*x for x in row] for row in side] for side in exact]
            dims, ddims = dimensions(code, n), dimensions(dual, n)
            for last in ((n+1)//2, n):
                system = Constraints(*inflated, dims, ddims, len(basis), last, spectra=spectra,
                                     ones=tuple((1 << n)-1 in words for words in (code, dual)))
                values = []
                for side, h, u in system.keys:
                    count = exact[side][h-1][u]//rank_total(h, 4, h)
                    values.append(Q(count, system.caps[side][h-1][u]))
                for row, rhs in system.inequalities:
                    self.assertLessEqual(sum(a*values[i] for i, a in row.items()), rhs)
                for row, rhs in system.equalities:
                    self.assertEqual(sum(a*values[i] for i, a in row.items()), rhs)
                for h in range(1, min(4, len(basis))+1):
                    result = system.solve(h, last)
                    count = exact[0][h-1][last]//rank_total(h, 4, h)
                    self.assertGreaterEqual(result['cap'], count)
                    self.assertLessEqual(result['cap'], system.caps[0][h-1][last])
                with self.assertRaises(ValueError):
                    system.solve(4, last+1)

    def test_exact_equality_dual(self):
        equalities = [normalize({0: 1, 1: 1}, 1)]
        upper, correction = verify_dual([], equalities, {0: Q(1), 1: Q(1)}, [], [Q(2)], 2)
        self.assertEqual(upper, 1)
        self.assertEqual(correction, 0)

    def test_residual_dual_composition(self):
        inequalities = [normalize({0: 1, 1: 2}, 1)]
        objective = {0: Q(3), 1: Q(1)}
        budget, residual = dual_residual(inequalities, [], objective, [Q(1)], [], 2)
        before, _ = verify_dual(inequalities, [], objective, [Q(1)], [], 2)
        addition, _ = verify_dual(inequalities, [], dict(enumerate(residual)), [Q(11)], [], 2)
        after, _ = verify_dual(inequalities, [], objective, [Q(12)], [], 2)
        self.assertEqual(after, budget+addition)
        self.assertLess(after, before)
        self.assertEqual(after, 3)

    def test_repaired_solver_on_small_code(self):
        n, basis = 4, [3, 5]
        code = span(basis)
        dual = [w for w in range(1 << n) if all((w & b).bit_count() % 2 == 0 for b in basis)]
        exact = [exact_tuple_cdfs(words, n, 4) for words in (code, dual)]
        inflated = [[[2*x for x in row] for row in side] for side in exact]
        system = Constraints(*inflated, dimensions(code, n), dimensions(dual, n), len(basis), n)
        for h in (1, 2):
            normal = system.solve(h, 3)
            repaired = system.solve(h, 3, repair_rounds=2)
            actual = exact[0][h-1][3]//rank_total(h, 4, h)
            self.assertGreaterEqual(repaired['cap'], actual)
            self.assertLessEqual(repaired['cap'], normal['cap'])
            if 'dual_witness' in repaired:
                self.assertEqual(system.verify_witness(repaired['dual_witness']), repaired['cap'])
                invalid = dict(repaired['dual_witness'], y=[[0, '-1']])
                with self.assertRaises(ValueError):
                    system.verify_witness(invalid)
                invalid = dict(repaired['dual_witness'], equality_count=-1)
                with self.assertRaises(ValueError):
                    system.verify_witness(invalid)
        for invalid in (-1, 4, True):
            with self.assertRaises(ValueError):
                system.solve(1, 3, repair_rounds=invalid)

    def test_invalid_objective_index(self):
        with self.assertRaises(ValueError):
            verify_dual([], [], {2: Q(1)}, [], [], 2)


if __name__ == '__main__':
    unittest.main()

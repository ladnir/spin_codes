"""Independent map/interface checks for the paired t64/s16 proof candidate."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import re
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'proof'))
import disjoint_pair_maps as maps
from flint import arb, ctx


class DisjointPairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 256
        cls.wrapper, cls.record = maps.prepare()
        cls.physical = maps.screen.q1.kernel_t64.authenticate(cls.wrapper)
        cls.rows = [int(row, 16) for row in cls.record['expansion_rows_hex']]
        # Independent truth-table calculation from the declared ANF terms.
        cls.columns = []
        for x in range(64):
            bits = [1, *[(x >> i) & 1 for i in range(6)]]
            bits += [sum(((x >> i) & 1)*((x >> j) & 1) for i, j in group) % 2
                     for group in maps.GROUPS]
            cls.columns.append(sum(value << i for i, value in enumerate(bits)))
        cls.repo = HERE.parents[3]

    def test_declared_pairs_partition_quadratic_terms(self):
        terms = [tuple(pair) for group in maps.GROUPS for pair in group]
        self.assertEqual(sorted(terms), list(combinations(range(6), 2)))
        self.assertEqual(sorted(map(len, maps.GROUPS)), [1]*3+[2]*6)
        for group in maps.GROUPS:
            if len(group) == 2:
                self.assertFalse(set(group[0]) & set(group[1]))

    def test_truth_tables_and_generated_cpp_columns(self):
        self.assertEqual(self.columns, self.record['feedback_columns'])
        source = (self.repo / 'spin/experiments/k16_codesign_100us/kernel/T64PairedMap.h').read_text()
        literal = re.search(r'columns\[64\]=\{([^}]+)\}', source).group(1)
        self.assertEqual([int(value, 16) for value in literal.split(',')], self.columns)
        expected_terms = [(0,), *[(1 << i,) for i in range(6)],
                          *[tuple((1 << i) | (1 << j) for i, j in group) for group in maps.GROUPS]]
        statements = re.findall(r'out\[(\d+)\]=([^;]+);', source)
        self.assertEqual(len(statements), 16)
        for (index, expression), terms in zip(statements, expected_terms):
            self.assertEqual(tuple(map(int, re.findall(r'z\[(\d+)\]', expression))), terms)

    def test_full_state_images_and_CA_zero(self):
        weights = Counter()
        for state in range(65536):
            weights[sum((state & column).bit_count() % 2 for column in self.columns)] += 1
        self.assertEqual({str(w): count for w, count in sorted(weights.items())},
                         self.record['expansion_spectrum'])
        self.assertEqual(sum(weights.values()), 65536)
        self.assertEqual(weights[0], 1)
        self.assertTrue(all((a & b).bit_count() % 2 == 0 for a in self.rows for b in self.rows))

    def test_macro_is_two_chronological_physical_steps(self):
        self.assertEqual(self.wrapper['physical_steps'], 2)
        self.assertEqual(self.wrapper['physical_step_bits'], 64)
        self.assertEqual(self.wrapper['macro_step_bits'], 128)
        self.assertEqual(self.wrapper['physical_windows'], 16)
        self.assertEqual(self.wrapper['macro_windows'], 32)
        self.assertEqual(self.wrapper['state_continuity'], 'retained_between_halves')
        q1 = maps.screen.q1
        kernel = q1.kernel_t64
        z = (-kernel.aq(Fraction(1, 200))).exp()
        local = kernel.sparse_kernel.outward_at_z(self.physical, z)
        local = kernel.kernel_birth_density.refine_local(self.physical, local, z, Fraction(1, 2))
        direct = q1.placement(local, epochs=32, windows=16, maximum_groups=2, rounding=q1.rounded)
        macro = kernel.local_operators(self.wrapper, Fraction(1, 200))
        composed = q1.placement(macro, epochs=16, windows=32, maximum_groups=2, rounding=q1.rounded)
        for occupancy in range(3):
            for i in range(direct[occupancy].nrows()):
                for j in range(direct[occupancy].ncols()):
                    self.assertLess(abs(direct[occupancy][i, j]-composed[occupancy][i, j]), arb(2)**-220)


if __name__ == '__main__':
    unittest.main()

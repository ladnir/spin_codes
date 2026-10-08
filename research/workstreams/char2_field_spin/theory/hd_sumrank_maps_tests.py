"""Exact finite-field regressions for the proof-side sum-rank inner."""

from itertools import combinations
import random
import unittest

from hd_blockwise_shared_transfer import mm, power, rank
from hd_sumrank_maps import construction, mask_ranks, Q, transpose


class SumRankMapsTests(unittest.TestCase):
    def test_field_and_nilpotence(self):
        m = construction()
        f, a, c = (m[key] for key in ('F', 'A', 'C'))
        self.assertEqual(mm(f, f, Q), tuple(tuple(int(i == j) for j in range(64)) for i in range(64)))
        self.assertEqual(mm(c, a, Q), ((0,)*4,)*4)
        self.assertEqual(mm(c, mm(f, a, Q), Q), ((0,)*4,)*4)
        self.assertEqual(len({power(x, 4369, Q) for x in m['representatives']}), 8)

    def test_local_cauchy_minors(self):
        f = construction()['local_f']
        # All minors through order two, then random larger minors below.
        for size in (1, 2):
            for rows in combinations(range(8), size):
                for cols in combinations(range(8), size):
                    self.assertEqual(rank(tuple(tuple(f[i][j] for j in cols) for i in rows), Q), size)

    def test_sumrank_after_basefield_basis_changes(self):
        rng = random.Random(318921)
        m = construction()
        for _ in range(80):
            dims = [0]*8
            for __ in range(rng.randrange(1, 10)):
                i = rng.randrange(8)
                dims[i] = min(4, dims[i]+1)
            columns = []
            for block, dimension in zip(m['blocks'], dims):
                if not dimension:
                    continue
                while True:
                    local = tuple(tuple(rng.choice(m['subfield']) for j in range(dimension)) for i in range(4))
                    if rank(local, Q) == dimension:
                        break
                columns.extend(transpose(mm(block, local, Q)))
            self.assertEqual(rank(transpose(columns), Q), min(4, sum(dims)))

    def test_actual_mask_ranks(self):
        rng = random.Random(82653)
        for _ in range(90):
            # Include both small and dense masks and complementary cases.
            potential = sorted(rng.sample(range(64), rng.randrange(65)))
            zeros = sorted(rng.sample(range(64), rng.randrange(65)))
            receipt = mask_ranks(potential, zeros)
            for kind in ('quiet', 'flat', 'input'):
                self.assertEqual(receipt[kind+'_actual'], receipt[kind+'_expected'], (potential, zeros, receipt))

    def test_larger_state_masks_and_nilpotence(self):
        rng=random.Random(912371)
        for state in (5,6,8,12,16):
            m=construction(state)
            zero=((0,)*state,)*state
            self.assertEqual(mm(m['C'],m['A'],Q),zero)
            self.assertEqual(mm(m['C'],mm(m['F'],m['A'],Q),Q),zero)
            for _ in range(20):
                potential=sorted(rng.sample(range(64),rng.randrange(65)))
                zeros=sorted(rng.sample(range(64),rng.randrange(65)))
                receipt=mask_ranks(potential,zeros,state)
                for kind in ('quiet','flat','input'):
                    self.assertEqual(receipt[kind+'_actual'],receipt[kind+'_expected'],(state,receipt))


if __name__ == '__main__':
    unittest.main()

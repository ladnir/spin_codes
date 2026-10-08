"""Exact tests for zero-set expansion of a common-rank potential."""

from fractions import Fraction as F
from itertools import product
from math import comb
import unittest

from rank_packet_operator_tests import MUL16, rank16
from sparse_multirow_zero_set import local_rank_transfer, transfer


class ZeroSetRankTests(unittest.TestCase):
    def test_positive_zero_set_identity(self):
        for width in range(1, 9):
            for weight in range(width+1):
                for y in (F(1, 3), F(9, 10), F(1)):
                    value = sum(comb(width-weight, selected)*y**(width-selected)
                                *(1-y)**selected for selected in range(width-weight+1))
                    self.assertEqual(value, y**weight)

    def test_input_support_sum_is_exact(self):
        x, y, tau = F(2, 7), F(4, 5), F(13, 2)
        whole = transfer(x, y, tau)
        for r in range(2):
            for k in range(2):
                reference = sum(comb(16, w)*x**w*local_rank_transfer(w, y, tau)[r][k]
                                for w in range(17))
                self.assertEqual(whole[r][k], reference)

    def test_small_cauchy_packet_common_quotient(self):
        # Four GF4 labels in GF16; two residue nodes outside that subfield.
        # A[u,j]=u+node_j+1 is row-scaled Cauchy; C[j,u]=u+node_j.
        labels, nodes, eta = (0, 1, 6, 7), (2, 3), 2
        cases = 0
        for state in ([[0, 0], [0, 0]], [[1, 0], [2, 0]], [[1, 0], [0, 1]]):
            incoming = int(rank16(state) > 0)
            for weight in (0, 1):
                form = [1, 1]
                histogram = [[[0]*3 for _ in range(5)] for _ in range(2)]
                for d0, d1 in product(range(1, 16), repeat=2):
                    refreshed = [[MUL16[d][a] for a in row] for d, row in zip((d0, d1), state)]
                    for dx in (range(1, 16) if weight else (1,)):
                        x = [MUL16[dx][a] for a in form] if weight else [0, 0]
                        following = [[a ^ MUL16[node][b] for a, b in zip(row, x)]
                                     for node, row in zip(nodes, refreshed)]
                        outputs = []
                        for label in labels:
                            value = [MUL16[eta][MUL16[label ^ nodes[0] ^ 1][a]
                                                ^ MUL16[label ^ nodes[1] ^ 1][b]]
                                     for a, b in zip(*refreshed)]
                            if label == 0:
                                value = [a ^ b for a, b in zip(value, x)]
                            outputs.append(value)
                        for mask in range(16):
                            selected = mask.bit_count()
                            equations = [outputs[j] for j in range(4) if mask >> j & 1]
                            for returning in (False, True):
                                eqs = equations+(following if returning else [])
                                er = rank16(eqs)
                                if weight and rank16(eqs+[form]) == er:
                                    continue
                                if incoming and rank16(eqs+state) == er:
                                    continue
                                if not returning and rank16(eqs+following) == er:
                                    continue
                                histogram[0 if returning else 1][selected][er] += 1
                for y in (F(1, 2), F(9, 10)):
                    for tau in (F(2), F(15)):
                        upper = local_rank_transfer(weight, y, tau, q=16, width=4, state=2)
                        for destination in range(2):
                            actual = sum(F(count, 15**(2+weight))*y**(4-z)*(1-y)**z/tau**r
                                         for z, ranks in enumerate(histogram[destination])
                                         for r, count in enumerate(ranks))
                            self.assertLessEqual(actual, upper[incoming][destination],
                                                 (state, weight, tau, y, destination))
                        cases += 1
        self.assertEqual(cases, 24)


if __name__ == '__main__':
    unittest.main()

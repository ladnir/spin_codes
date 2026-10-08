"""Exact small-field checks of the Cauchy cap operator and its cone."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
import unittest

import hd_cauchy_rank as cover
from diagonal_refresh_rank import MUL,INV,matrix_rank as rank


def dot(left,right):
    value = 0
    for a,b in zip(left,right):
        value ^= MUL[a][b]
    return value


def matmul(left,right):
    return [[dot(row,column) for column in zip(*right)] for row in left]


def inverse(matrix):
    size = len(matrix)
    work = [list(row)+[int(i==j) for j in range(size)] for i,row in enumerate(matrix)]
    for col in range(size):
        pivot = next(i for i in range(col,size) if work[i][col])
        work[col],work[pivot] = work[pivot],work[col]
        work[col] = [MUL[INV[work[col][col]]][x] for x in work[col]]
        for row in range(size):
            if row != col:
                factor = work[row][col]
                work[row] = [a^MUL[factor][b] for a,b in zip(work[row],work[col])]
    return [row[size:] for row in work]


def maps():
    # Full 4x4 Cauchy over GF8, then normalize the TWO bottom rows.
    raw = [[INV[a^b] for b in range(4,8)] for a in range(4)]
    f = [row[:2] for row in raw[:2]]
    a = [row[2:] for row in raw[:2]]
    c = matmul(inverse([row[2:] for row in raw[2:]]),
               [row[:2] for row in raw[2:]])
    return f,a,c


def exact_histogram(state_forms,input_form=None,input_position=0,extra_history=()):
    """Enumerate current diagonals, equations and actual-zero/live covers.

    Each history item is a list of linear forms whose joint map must remain
    nonzero. Testing this necessary condition may overcount feasible live
    witnesses, which is safe for testing an upper cover.
    """
    f,a,c = maps()
    dimension = len(state_forms[0])
    zero = [0]*dimension
    incoming_rank = rank(state_forms)
    weight = int(input_form is not None)
    history = list(extra_history)
    if incoming_rank:
        history.append(state_forms)
    if weight:
        history.append([input_form])
    histogram = Counter()
    for ds in product(range(1,8),repeat=2):
        refreshed = [[MUL[d][x] for x in row] for d,row in zip(ds,state_forms)]
        for dx in (range(1,8) if weight else (1,)):
            x = [MUL[dx][v] for v in input_form] if weight else zero
            out = [[dot(a[j],[refreshed[t][v] for t in range(2)])
                    ^ MUL[f[j][input_position]][x[v]] for v in range(dimension)]
                   for j in range(2)]
            successor = [[refreshed[j][v]^MUL[c[j][input_position]][x[v]]
                          for v in range(dimension)] for j in range(2)]
            for mask in range(4):
                equations = [out[j] for j in range(2) if mask>>j&1]
                gained = rank(equations)
                live = all(rank(equations+list(item)) > gained for item in history)
                following = rank(equations+successor)-gained
                if live and following:
                    if not weight:
                        assert following == incoming_rank-gained
                    histogram[(mask.bit_count(),gained,following)] += 1
                returned = equations+successor
                return_gain = rank(returned)
                if all(rank(returned+list(item)) > return_gain for item in history):
                    histogram[(mask.bit_count(),return_gain,0)] += 1
    return incoming_rank,weight,histogram,7**(2+weight)


def weighted(histogram,denominator,y,tau):
    result = [F(0)]*3
    for (selected,gain,following),count in histogram.items():
        result[following] += (F(count,denominator)*y**(2-selected)
                              *(1-y)**selected/tau**gain)
    return result


class CauchyRankTests(unittest.TestCase):
    def test_small_full_cauchy_and_schur_profiles(self):
        f,a,c = maps()
        for matrix in (f,a,c):
            self.assertTrue(all(x for row in matrix for x in row))
            self.assertEqual(rank(matrix),2)
        returned = [[left^right for left,right in zip(l,r)]
                    for l,r in zip(f,matmul(a,c))]
        self.assertTrue(all(x for row in returned for x in row))
        self.assertEqual(rank(returned),2)
        for pivot in range(2):
            for column in range(2):
                other = 1-pivot
                ratio = MUL[f[other][column]][INV[f[pivot][column]]]
                residual = [l^MUL[ratio][r] for l,r in zip(a[other],a[pivot])]
                self.assertTrue(all(residual))

    def test_quiet_and_singleton_joint_successor_domination(self):
        configurations = (
            ([[0,0],[0,0]],[1,0]),
            ([[1,0],[0,0]],None),
            ([[1,0],[1,0]],None),
            ([[1,0],[0,1]],None),
            ([[1,0],[0,0]],[1,0]),
            ([[1,0],[0,0]],[0,1]),
            ([[1,0],[2,0]],[1,0]),
            ([[1,0],[2,0]],[0,1]),
            ([[1,0],[0,1]],[1,0]),
            ([[1,0],[0,1]],[1,1]),
            ([[1,0,0],[0,1,0]],[0,0,1]),
        )
        checked = 0
        for state,input_form in configurations:
            for position in (0,1):
                for extra in ((), ([[1]+[0]*(len(state[0])-1)],)):
                    incoming,weight,histogram,denominator = exact_histogram(
                        state,input_form,position,extra)
                    for y,tau in product((F(1,2),F(1)),(F(1),F(2),F(7))):
                        actual = weighted(histogram,denominator,y,tau)
                        matrix = cover.local_matrix(weight,y,tau,q=8,state=2,width=2)
                        for cap in (range(incoming,3) if incoming else (0,)):
                            for ray in cover.cone_generators(2):
                                self.assertLessEqual(sum(x*v for x,v in zip(actual,ray)),
                                    sum(x*v for x,v in zip(matrix[cap],ray)),
                                    (state,input_form,position,extra,y,tau,cap,ray))
                                checked += 1
        self.assertEqual(checked,1224)

    def test_actual_four_state_cap_cone_and_birth(self):
        for q,y,tau in product((256,65536),(F(1,3),F(1)),(F(1),F(2),F(7))):
            for weight in range(5):
                matrix = cover.local_matrix(weight,y,tau,q)
                self.assertTrue(all(value >= 0 for row in matrix for value in row))
                for ray in cover.cone_generators():
                    values = [sum(a*b for a,b in zip(row,ray)) for row in matrix]
                    self.assertTrue(all(values[j] <= values[j+1] for j in range(1,4)))
                if weight == 1:
                    self.assertEqual(matrix[0],(0,y**4,0,0,0))
            for cap in range(1,5):
                self.assertLessEqual(cover.singleton_return(cap,y,tau,q),y**4/tau**4)

    def test_log_matrix_polynomial_composition(self):
        import numpy as np
        y,tau = F(1,2),F(2)
        matrices = np.asarray([cover.local_matrix(w,y,tau) for w in range(5)],float)
        with np.errstate(divide='ignore'):
            logs = np.log(matrices)
        actual = np.exp(cover.log_polynomial_product(logs,logs,8))
        expected = np.zeros_like(actual)
        for i,left in enumerate(matrices):
            for j,right in enumerate(matrices):
                expected[i+j] += left@right
        self.assertTrue(np.allclose(actual,expected,rtol=1e-13,atol=1e-15))
        for coefficient in actual:
            for ray in cover.cone_generators():
                values = coefficient@np.asarray(ray,float)
                self.assertTrue(all(values[j] <= values[j+1]+1e-14 for j in range(1,4)))


if __name__ == '__main__':
    unittest.main()

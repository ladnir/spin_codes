"""Tiny exact oracles only; root runs remotely with Python -B/BLAS1."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest
from unittest.mock import patch

import hd_field128_wide_state_v60 as candidate


def mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 0x13
    return out


def mv(a, x):
    out = []
    for row in a:
        value = 0
        for b, c in zip(row, x):
            value ^= mul(b, c)
        out.append(value)
    return tuple(out)


def rank(a):
    a = [list(row) for row in a]
    r = 0
    for c in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = next(x for x in range(1, 16) if mul(x, a[r][c]) == 1)
        a[r] = [mul(inv, x) for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                factor = a[i][c]
                a[i] = [x ^ mul(factor, y) for x, y in zip(a[i], a[r])]
        r += 1
    return r


class Tests(unittest.TestCase):
    def test_nonzero_CA_tiny_principal_and_return_oracle(self):
        # Width3/state2 is beyond half width. A,C are MDS over GF16;
        # CA=[[0,1],[1,0]]!=0. gamma=2 has degree4 over their GF2 entries.
        a = ((1, 0), (0, 1), (1, 1))
        c = tuple(zip(*a))
        self.assertEqual(tuple(mv(c, col) for col in zip(*a)), ((0, 1), (1, 0)))
        p = tuple(tuple(mv(a, col)) for col in zip(*c))
        p = tuple(zip(*p))
        l = tuple(tuple(int(i == j) ^ mul(2, p[i][j]) for j in range(3))
                  for i in range(3))
        for size in range(4):
            for s in combinations(range(3), size):
                self.assertEqual(rank([[l[i][j] for j in s] for i in s]), size)
                rows = []
                for values in product(range(16), repeat=size):
                    x = [0]*3
                    for index, value in zip(s, values):
                        x[index] = value
                    syndrome, y = mv(c, x), mv(l, x)
                    rows.append((bool(any(syndrome)), tuple(x), y))
                for z in (F(1, 2), F(1, 3)):
                    quiet = sum((z**sum(v != 0 for v in x)
                                 for nonzero, x, y in rows if not nonzero), F(0))/16**size
                    birth = sum((z**sum(v != 0 for v in x)
                                 for nonzero, x, y in rows if nonzero), F(0))/16**size
                    ret = sum((z**sum(v != 0 for v in y)
                               for nonzero, x, y in rows if nonzero), F(0))/16**size
                    base = candidate._moments(size, z, 16, 3, 2, False)
                    refined = candidate._moments(size, z, 16, 3, 2, True)
                    self.assertEqual((quiet, birth), base[:2])
                    self.assertLessEqual(ret, refined[2])
                    self.assertLessEqual(refined[2], base[2])
                    # Active Y marginal and its syndrome-zero restriction
                    # agree with uniform X, despite CA != 0.
                    active = sum((z**sum(y[i] != 0 for i in s)
                                  for nonzero, x, y in rows if nonzero), F(0))/16**size
                    self.assertEqual(active, birth)
                for dc in range(size+1):
                    for d in combinations(s, dc):
                        outside = tuple(i for i in range(3) if i not in s)
                        for ec in range(len(outside)+1):
                            for e in combinations(outside, ec):
                                actual = F(sum(nonzero and all(y[i] == 0 for i in d+e)
                                               for nonzero, x, y in rows), 16**size)
                                self.assertLessEqual(actual,
                                    candidate.joint_probability(size, dc, ec, 16, 3, 2))

    def test_independent_uniform_nonzero_live(self):
        a = ((1, 0), (0, 1), (1, 1))
        for s in ((), (0,), (2,)):
            weights = [0]*4
            for values in product(range(16), repeat=len(s)):
                x = [0]*3
                for i, value in zip(s, values):
                    x[i] = value
                for q in product(range(16), repeat=2):
                    if not any(q):
                        continue
                    emission = mv(a, q)
                    weights[sum(v ^ w != 0 for v, w in zip(x, emission))] += 1
            for z in (F(1, 2), F(2, 3)):
                direct = sum((count*z**w for w, count in enumerate(weights)), F(0))
                direct /= 16**len(s)*(16**2-1)
                self.assertEqual(direct, candidate._moments(len(s), z, 16, 3, 2, False)[3])

    def test_generic32_control(self):
        for j in (0, 1, 2):
            old = candidate.maps.exact_moments(j, F(1, 2), 1 << 128, 64, 32, True)
            new = candidate.exact_moments(j, F(1, 2), 32, tail_refine=False)
            self.assertEqual(new, old)

    def test_public_packet_normalization_and_identity(self):
        local, z = (F(1, 3), F(1, 5), F(1, 7), F(1, 11)), F(1, 2)
        self.assertEqual(len({candidate.construction(s) for s in candidate.STATES}), 3)
        for state in candidate.STATES:
            d = candidate.FIELD_ORDER**state-1
            density = F(d, (candidate.FIELD_ORDER-1)**state)
            with patch.object(candidate, 'exact_moments', return_value=local):
                self.assertEqual(candidate.exact_packet(1, z, state),
                    (local[0], density*local[1], local[2]/d, density*local[3]))
                self.assertEqual(candidate.exact_packet(0, z, state)[3], z**(65-state))
            self.assertEqual(candidate.principal_bad_gamma_upper(state),
                F(state*(1 << 64), (1 << 128)-(1 << 64)))
            with patch.object(candidate, 'source_pins', return_value={'mock': 'pin'}):
                meta = candidate.metadata(state)
            for key in ('CA_zero_assumed', 'involution_assumed', 'old_numerical_bounds_reused',
                        'setup_risk_included_in_packet', 'whole_code_certificate'):
                self.assertFalse(meta[key])

    def test_guards_and_source_identity(self):
        pins = candidate.source_pins()
        candidate.check_source_pins(pins)
        self.assertEqual(pins[candidate.tail.NAME], candidate.FROZEN[candidate.tail.NAME])
        with patch.dict(candidate.FROZEN, {candidate.tail.NAME: '0'*64}):
            with self.assertRaises(RuntimeError):
                candidate.verify_frozen()
        for state in (True, 16, 33, 64):
            with self.assertRaises(ValueError):
                candidate.exact_moments(0, F(1, 2), state)
        for j, z in ((True, F(1, 2)), (65, F(1, 2)), (1, .5), (1, True), (1, 0), (1, 1)):
            with self.assertRaises(ValueError):
                candidate.exact_moments(j, z, 36)


if __name__ == '__main__':
    unittest.main()

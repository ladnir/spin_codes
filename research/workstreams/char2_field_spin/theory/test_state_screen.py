"""Exact small-field audits of the wider-state floating screening envelope.

No numerical screen is promoted to a certificate by these tests. GF(4)
enumeration checks every input support, input value, and refreshed nonzero
state for MDS maps of state dimension 1 through 4. Independent rational
matrix products check the regional conditioning and whole-code normalization.
"""

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import math
import unittest
from unittest.mock import patch

import state_screen as screen


def mul4(a, b):
    """Polynomial multiplication modulo X^2+X+1, independent of other models."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def dot4(a, b):
    result = 0
    for x, y in zip(a, b):
        result ^= mul4(x, y)
    return result


def mul16(a, b):
    """Polynomial multiplication modulo X^4+X+1, independent of the model."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 19
    return result


GF16_PRODUCT = tuple(tuple(mul16(a, b) for b in range(16)) for a in range(16))
GF16_A = tuple((1, j, mul16(j, j), mul16(mul16(j, j), j)) for j in range(16))


def dot16(a, b):
    result = 0
    for x, y in zip(a, b):
        result ^= GF16_PRODUCT[x][y]
    return result


@lru_cache(maxsize=None)
def enumerated(width, state):
    a = []
    for point in range(width):
        row = [1]
        for _ in range(1, state):
            row.append(mul4(row[-1], point))
        a.append(row)
    refreshed = [(r, tuple(dot4(row, r) for row in a))
                 for r in product(range(4), repeat=state) if any(r)]
    counts = {}
    for x in product(range(4), repeat=width):
        mask = tuple(i for i, value in enumerate(x) if value)
        marked = counts.setdefault(mask, [[[0]*(width+1) for _ in range(2)]
                                          for _ in range(2)])
        cx = tuple(dot4(tuple(row[j] for row in a), x) for j in range(state))
        marked[0][int(any(cx))][len(mask)] += 1
        for r, ar in refreshed:
            next_active = int(any(u ^ v for u, v in zip(r, cx)))
            weight = sum(u != v for u, v in zip(x, ar))
            marked[1][next_active][weight] += 1
    return counts


def marked_matrix(counts, weight, state, z):
    return [[sum(F(n)*z**d for d, n in enumerate(counts[i][j]))
             / (3**weight * (4**state-1 if i else 1))
             for j in range(2)] for i in range(2)]


def support_count(length, dimension, weight, q):
    # The alternate (q-1)*sum form of the MDS exact-support count.
    excess = weight-(length-dimension)-1
    return ((q-1)*sum((-1)**j * math.comb(weight-1, j) * q**(excess-j)
                     for j in range(excess+1))) if excess >= 0 else 0


def rational_local(z, width, state, q=256):
    output = []
    for w in range(width+1):
        count = 1 if w == 0 else support_count(w, max(0, w-state), w, q)
        p0 = F(count, (q-1)**w)
        emission = F(0)
        f = F(1, q-1) + F(q-2, q-1)*z
        for zero_set in range(1 << width):
            zeros = zero_set.bit_count()
            if zeros >= state:
                continue
            active_zeros = (zero_set & ((1 << w)-1)).bit_count()
            inactive_nonzeros = width-w-(zeros-active_zeros)
            multiplicity = support_count(width, state, width-zeros, q)
            emission += (F(multiplicity, q**state-1)
                         * z**(inactive_nonzeros+active_zeros)
                         * f**(w-active_zeros))
        ret = (1-p0)/(q**state-1)
        output.append([[z**w*p0, z**w*(1-p0)],
                       [ret*z**max(0, width-state+1-w), emission-ret*z**width]])
    return output


def add(a, b):
    return [[a[i][j]+b[i][j] for j in range(2)] for i in range(2)]


def scale(a, value):
    return [[a[i][j]*value for j in range(2)] for i in range(2)]


def matmul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


def rational_region(z, h, rows, width, state, q=256):
    # Literal conditioning on each h-subset, not polynomial convolution or
    # Bernoulli rescaling. Used only for tiny regions.
    matrices = rational_local(z, width, state, q)
    result = [[F(0), F(0)], [F(0), F(0)]]
    for support in combinations(range(rows), h):
        counts = [0]*(rows//width)
        for location in support:
            counts[location//width] += 1
        value = [[F(1), F(0)], [F(0), F(1)]]
        for w in counts:
            value = matmul(value, matrices[w])
        result = add(result, value)
    return scale(result, F(1, math.comb(rows, h)))


class StateScreenTests(unittest.TestCase):
    def assertMatrixClose(self, actual, expected):
        for i in range(2):
            for j in range(2):
                self.assertTrue(math.isclose(float(actual[i][j]), float(expected[i][j]),
                                             rel_tol=3e-12, abs_tol=1e-14),
                                (i, j, actual[i][j], expected[i][j]))

    def test_gl2_refresh_is_uniform_nonzero(self):
        counts = {}
        source = (2, 3)
        matrices = 0
        for a, b, c, d in product(range(4), repeat=4):
            if mul4(a, d) == mul4(b, c):
                continue
            image = (dot4((a, b), source), dot4((c, d), source))
            counts[image] = counts.get(image, 0)+1
            matrices += 1
        self.assertEqual(matrices, 180)
        self.assertEqual(len(counts), 15)
        self.assertEqual(set(counts.values()), {12})
        self.assertNotIn((0, 0), counts)

    def test_every_gf4_local_support_has_entrywise_upper_envelope(self):
        for state in range(1, 5):
            counts = enumerated(4, state)
            for z in (F(1, 7), F(2, 3), F(1)):
                exact_envelopes = rational_local(z, 4, state, 4)
                floating = screen.local(float(z), 4, state, 4)
                for support, marked in counts.items():
                    w = len(support)
                    exact = marked_matrix(marked, w, state, z)
                    envelope = exact_envelopes[w]
                    for i in range(2):
                        for j in range(2):
                            self.assertLessEqual(exact[i][j], envelope[i][j])
                            self.assertGreaterEqual(envelope[i][j], 0)
                    self.assertEqual(exact[0], envelope[0])
                    self.assertMatrixClose(floating[w], envelope)
                    if z == 1:
                        self.assertEqual(exact, envelope)
                        self.assertEqual([sum(row) for row in exact], [1, 1])

    def test_exact_emission_moment_and_return_probability(self):
        for state in range(1, 5):
            counts = enumerated(4, state)
            z = F(2, 3)
            for support, marked in counts.items():
                w = len(support)
                p0 = F(screen.birth_zero_count(w, state, 4), 3**w)
                unmarked = marked_matrix(marked, w, state, F(1))
                self.assertEqual(unmarked[0][0], p0)
                self.assertEqual(unmarked[1][0], (1-p0)/(4**state-1))
                exact = marked_matrix(marked, w, state, z)
                envelope = rational_local(z, 4, state, 4)[w]
                # Add back the lower return term subtracted by the envelope.
                g = envelope[1][1] + unmarked[1][0]*z**4
                self.assertEqual(g, sum(exact[1]))

    def test_gf256_combinatorial_counts(self):
        for width in (8, 16):
            for state in (1, 2, 4, 8):
                patterns = screen.zero_patterns(width, state, 256)
                self.assertEqual(sum(math.comb(width, j)*count
                                     for j, count in enumerate(patterns)), 256**state-1)
                for j, count in enumerate(patterns):
                    self.assertEqual(count, support_count(width, state, width-j, 256))
                for w in range(width+1):
                    count = screen.birth_zero_count(w, state, 256)
                    self.assertGreaterEqual(count, 0)
                    self.assertLessEqual(count, 255**w)
                    self.assertEqual(count, 1 if w == 0 else
                                     support_count(w, max(0, w-state), w, 256))

    def test_region_conditioning_against_literal_support_average(self):
        z = F(2, 3)
        original = screen.local
        for q in (4, 256):
            for state in (1, 2, 3, 4):
                with patch.object(screen, 'local',
                                  side_effect=lambda z, w, s: original(z, w, s, q)):
                    result = screen.regions(float(z), 3, 8, 4, state)
                for h in range(4):
                    reference = rational_region(z, h, 8, 4, state, q)
                    self.assertMatrixClose(result[h], reference)

    def test_coupled_normalization_matches_direct_rational_expression(self):
        # Tiny dimensions avoid underflow and allow exact rational replay of
        # the selected floating z. This checks scaling, not outward rounding.
        active, k, width, state, b, dim = 2, 16, 4, 3, 4, 2
        result = screen.coupled(active, k, width, state, b, dim, 0.1, samples=3)
        self.assertIsNotNone(result)
        self.assertFalse(result['certificate'])
        z = F.from_float(result['z'])
        rows = k//dim
        matrix = [[F(0), F(0)], [F(0), F(0)]]
        for h in range(active+1):
            matrix = add(matrix, scale(rational_region(z, h, rows, width, state),
                                       math.comb(active, h)*256**h))
        total = [[F(1), F(0)], [F(0), F(1)]]
        for _ in range(b):
            total = matmul(total, matrix)
        e = math.ceil(0.1*rows*b)-1
        bound = (math.comb(rows, active)*sum(total[0])
                 / (255*256**(active*(b-dim))*z**e))
        reference_bits = (math.log(bound.denominator)-math.log(bound.numerator))/math.log(2)
        self.assertAlmostEqual(result['margin_bits'], reference_bits, places=10)

    def test_marking_parameter_domain_is_explicit(self):
        for invalid in (0, -1, 1.01, float('nan'), float('inf'), -float('inf')):
            with self.assertRaises(ValueError):
                screen.local(invalid)

    def test_gf16_monomial_map_is_self_orthogonal(self):
        # The exact dense-return improvement uses this extra property; it is
        # not a consequence of the generic MDS interface alone.
        for i in range(4):
            for j in range(4):
                self.assertEqual(dot16(tuple(a[i] for a in GF16_A),
                                       tuple(a[j] for a in GF16_A)), 0)
        b = [[int(i == j) ^ dot16(GF16_A[i], GF16_A[j])
              for j in range(16)] for i in range(16)]
        for i in range(16):
            for j in range(16):
                self.assertEqual(dot16(b[i], tuple(b[k][j] for k in range(16))),
                                 int(i == j))

    def test_gf16_monomial_zero_patterns_by_enumeration(self):
        histogram = [0]*16
        for r in product(range(16), repeat=4):
            if any(r):
                histogram[sum(dot16(a, r) == 0 for a in GF16_A)] += 1
        exact = screen.zero_patterns(16, 4, 16)
        self.assertEqual(histogram[:4], [math.comb(16, i)*exact[i] for i in range(4)])
        self.assertEqual(histogram[4:], [0]*12)
        self.assertEqual(sum(histogram), 16**4-1)

    def test_gf16_singleton_return_polynomial(self):
        counts = [0]*17
        by_lane = []
        for j in range(16):
            lane_weights = set()
            for x in range(1, 16):
                r = tuple(GF16_PRODUCT[a][x] for a in GF16_A[j])
                y = tuple(dot16(a, r) ^ (x if ell == j else 0)
                          for ell, a in enumerate(GF16_A))
                # This is exactly the unique refreshed state causing return.
                self.assertEqual(tuple(r[i] ^ GF16_PRODUCT[GF16_A[j][i]][x]
                                       for i in range(4)), (0, 0, 0, 0))
                weight = sum(value != 0 for value in y)
                lane_weights.add(weight)
                counts[weight] += 1
            self.assertEqual(len(lane_weights), 1)
            by_lane.append(lane_weights.pop())
        self.assertEqual(counts, [0]*14+[30, 195, 15])
        self.assertEqual(by_lane[0], 15)
        self.assertEqual(by_lane[1], 16)
        exceptional = [j for j in range(2, 16) if by_lane[j] == 14]
        self.assertEqual(len(exceptional), 2)
        for j in exceptional:
            self.assertEqual(mul16(j, j) ^ j ^ 1, 0)
        from reduced_distance_certificate import local_model
        for z in (F(1, 7), F(2, 3), F(1)):
            actual = sum(count*z**weight for weight, count in enumerate(counts))
            actual /= 16*15*(16**4-1)
            expected = (13*z**15+z**16+2*z**14)/(16*(16**4-1))
            self.assertEqual(actual, expected)
            self.assertAlmostEqual(local_model(float(z), q=16, state=4)[1, 1, 0],
                                   float(expected), places=16)

    def test_dense_return_rank_bound_and_universal_row_sum(self):
        # GF4, t=4, s=2 has CA=0 and hence invertible B=I+AC. Check the
        # positive-subset rank upper bound against every exact input support.
        for z in (F(1, 7), F(2, 3), F(1)):
            f = z+(1-z)/3
            for support, marked in enumerated(4, 2).items():
                w = len(support)
                actual = marked_matrix(marked, w, 2, z)
                upper = sum(F(math.comb(4, h))*z**(4-h)*(1-z)**h
                            / 3**max(0, h-4+w) for h in range(5))/15
                self.assertLessEqual(actual[1][0], upper)
                for state in range(1, 5):
                    generic = marked_matrix(enumerated(4, state)[support], w, state, z)
                    self.assertLessEqual(sum(generic[0]), f**w)
                    self.assertLessEqual(sum(generic[1]), f**w)

    def test_directed_gf16_local_intervals_enclose_rational_envelope(self):
        from gf16_certificate import Interval, interval_local_envelope
        z = F(2, 3)
        reference = rational_local(z, 16, 4, 16)
        joint = (13*z**15+z**16+2*z**14)/(16*(16**4-1))
        reference[1][1][1] += z**16/(16**4-1)-joint
        reference[1][1][0] = joint
        for w in range(2, 17):
            dense = sum(F(math.comb(16, h), (16**4-1)*15**max(0, h-16+w))
                        * z**(16-h)*(1-z)**h for h in range(17))
            reference[w][1][0] = min(reference[w][1][0], dense)
        actual = interval_local_envelope(Interval.rational(2, 3))
        for w in range(17):
            for i in range(2):
                for j in range(2):
                    self.assertLessEqual(F(actual[w][i][j].lo), reference[w][i][j])
                    self.assertGreaterEqual(F(actual[w][i][j].hi), reference[w][i][j])

    def test_directed_source_dp_matches_independent_rational_replay(self):
        from decimal import Decimal
        from gf16_certificate import HIGH, upper_exact_source
        regional = [[[F(1, 3), F(1, 7)], [F(1, 5), F(1, 2)]],
                    [[F(1, 11), F(1, 9)], [F(1, 13), F(1, 4)]],
                    [[F(1, 17), F(1, 19)], [F(1, 23), F(1, 8)]]]
        rounded = [[HIGH.divide(Decimal(x.numerator), Decimal(x.denominator))
                    for row in matrix for x in row] for matrix in regional]
        for active in (1, 2):
            current = {((0,)*active, 0): F(1)}
            masks = list(product((0, 1), repeat=active))
            for _ in range(16):
                following = {}
                for (degrees, incoming), probability in current.items():
                    for mask in masks:
                        next_degrees = tuple(a+b for a, b in zip(degrees, mask))
                        for outgoing in range(2):
                            key = (next_degrees, outgoing)
                            following[key] = following.get(key, F(0)) + (
                                probability*regional[sum(mask)][incoming][outgoing])
                current = following
            exact = sum(probability*math.prod(support_count(16, 8, w, 16) for w in degrees)
                        for (degrees, _), probability in current.items())
            actual = F(upper_exact_source(rounded, active))
            self.assertGreaterEqual(actual, exact)
            self.assertLess(actual, exact*(1+F(1, 10**70)))


if __name__ == '__main__':
    unittest.main()
